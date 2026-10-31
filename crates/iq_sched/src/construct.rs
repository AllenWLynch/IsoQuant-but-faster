//! Block driver (design §7.1; order = 4.1 `GraphBasedModelConstructor.process`,
//! `model_construction.py:113-134`), composing U10b → U11 → U12 → U13-U19 → U20 counter.

use iq_annot::AnnotationIndex;
use iq_assign::AssignParams;
use iq_model::ends::Ends;
use iq_model::fl::{gene_id_ranks, FlAux, LiveAssigner};
use iq_model::{Ctx, DetectedKnown};
use iq_quant::TxCounter;
use iq_types::{BlockStore, GeneKey, IsoformIndex, ModelStore, Params, RefSlice, TranscriptModel};

/// Final models and counters of one block.
pub struct BlockModels {
    pub models: Vec<TranscriptModel>,
    pub counter: TxCounter,
}

/// Runs model construction and quantification for one collected block. `detected_known` is
/// threaded through the blocks of a work unit in genomic order. `reference` must cover the
/// block's reads (model boundaries ±9, FL-path introns).
pub fn construct_block(
    ann: &AnnotationIndex,
    params: &Params,
    index: &IsoformIndex,
    store: &BlockStore,
    reference: RefSlice<'_>,
    detected_known: &mut DetectedKnown,
) -> BlockModels {
    let mut counter = TxCounter::new();
    if store.is_empty() {
        return BlockModels { models: Vec::new(), counter };
    }
    // U10b: polyA (+TSS with --fl_data) predictions          parallel_workers.py:320-330
    let preds = iq_terminal::block_predictions(store, index, iq_terminal::bundled_models(), params);
    // U11 + U12                                                context.py:212-224
    let graph = iq_graph::build_graph(store, index, &preds, params);
    let paths = iq_graph::thread_paths(store, &graph, index, params);
    let assigner = LiveAssigner { params: AssignParams::new(params.matching.clone(), params.splice.clone(), false) };
    let ctx = Ctx { params, index, store, graph: &graph, paths: &paths, predictions: &preds, assigner: &assigner };
    let names: Vec<&str> = index
        .genes
        .iter()
        .map(|g| match g.key {
            GeneKey::Ref(gi) => ann.gene_id(gi),
            _ => "",
        })
        .collect();
    let gene_rank = gene_id_ranks(&names);

    let mut ms = ModelStore::new(store.len());
    iq_model::fl::construct_fl(&mut ms, detected_known, &ctx, &FlAux { reference, gene_rank: &gene_rank });
    iq_model::assign_based::construct_assignment_based(&mut ms, detected_known, &ctx);
    iq_model::filter::pre_filter(&mut ms, &ctx);
    iq_model::reassign::assign_reads_to_models(&mut ms, &ctx);
    iq_model::filter::filter(&mut ms, &ctx, &Ends::default());
    if params.splice_site_correction {
        iq_model::splice_site::correct_transcript_splice_sites(&mut ms, store, &reference);
    }
    iq_model::reassign::assign_reads_to_models(&mut ms, &ctx);
    iq_model::reassign::resolve_ownership(&mut ms);
    iq_model::joiner::join_genes(&mut ms.models, index);
    iq_model::final_assign::final_assignments(&ms, &ctx, |r, fa| {
        counter.add(&fa, store.view(r as usize).group());
    });
    BlockModels { models: ms.models, counter }
}

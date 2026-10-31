//! Serde round trip of every fixture payload type: a hand-written, Python-formatted example
//! line per `(unit, fn)` (written from SCHEMA.md v1) must re-encode byte-for-byte through
//! the typed DTOs. Also checks DTO ↔ core conversions for the shared types.

use iq_types::fixture::roundtrip::{reencode_block_line, reencode_case_line, reencode_index, reencode_params};
use iq_types::fixture::schema::*;
use iq_types::fixture::{IdMaps, TxSpace};
use iq_types::*;

const A: &str = r#"{"atype":10,"gene_atype":1,"matches":[{"gene":"gene:G1","tx":"transcript:T1","strand":"+","class":10,"penalty":0.0,"events":[{"t":11,"iso":[2147483648,2147483648],"read":[0,1],"info":0}]},{"gene":null,"tx":null,"strand":".","class":0,"penalty":1.2000000000000002,"events":[{"t":1003,"iso":[1073741823,1073741825],"read":[-1,2147483647],"info":-35}]}]}"#;
const M: &str = r#"{"id":"transcript3.8.nic","gene":"novel_gene_8_1","chr":"8","strand":"+","type":2,"exons":[[100,200],[300,400]],"intron_path":[[-21,100],[201,299],[-10,400]],"polya_confirmed":true,"tss_confirmed":false,"source":"IsoQuant","other_features":[[150,200,"CDS"]],"additional_info":[["k","v"]]}"#;
const P: &str = r#"{"delta":4,"minor_exon_extension":50,"major_exon_extension":300,"min_abs_exon_overlap":10,"min_rel_exon_overlap":0.2,"minimal_exon_overlap":5,"minimal_intron_absence_overlap":20,"micro_intron_length":50,"max_fake_terminal_exon_len":20,"max_intron_shift":30,"max_missed_exon_len":50,"max_intron_abs_diff":30,"max_intron_rel_diff":0.2,"max_suspicious_intron_abs_len":0,"max_suspicious_intron_rel_len":0.0,"resolve_ambiguous":20,"apa_delta":50,"count_exons":false,"correct_fuzzy_junctions":true,"correct_intron_shifts":false,"correct_skipped_exons":true,"correct_fake_terminal_exons":false,"correct_terminal_exons":false,"correct_microintron_retention":true}"#;
const PROF: &str = r#"{"gene_profile":[1,-1,-2,0],"read_profile":[1],"read_features":[[201,299]],"gene_profile_range":[0,2]}"#;
const ALN: &str = r#"{"qname":"r0","flag":16,"pos0":99,"ref_end0":400,"mapq":60,"is_reverse":true,"cigar":[[4,5],[0,101],[3,99],[0,101]],"seq":"ACGTN"}"#;
const B: &str = r#"{"sample":"S","slice":"sl","chr":"8","block_no":3}"#;
const G: &str = r#"{"clustered_introns":[[[201,299],5]],"intron_correction_map":[[[202,299],[201,299]]],"discarded_introns":[[500,600]],"outgoing":[[[201,299],[[-10,400]]],[[-21,100],[[201,299]]]],"incoming":[[[201,299],[[-21,100]]]],"starting_known_positions":[[[201,299],[100,101]]],"terminal_known_positions":[[[201,299],[400]]],"max_coverage":5}"#;

fn store() -> String {
    format!(
        r#"{{"models":[{M}],"transcript_read_ids":[["transcript3.8.nic",[0,2,2]]],"internal_counter":[["transcript3.8.nic",3],["transcript9.8.nnic",0]],"read_assignment_counts":[["r0",1],["r2",2]],"construction_assignment":[["r0","transcript3.8.nic"]],"read_assignment_counts_idx":[[0,1],[2,2]],"construction_assignment_idx":[[0,"transcript3.8.nic"]],"nonfl_known_reads":[["transcript:T1",[1]]],"detected_known_isoforms":["transcript:T1"],"id_counter":7}}"#
    )
}
fn read_record() -> String {
    format!(
        r#"{{"idx":0,"assignment_id":17,"name":"r0","chr":"8","region":[1,1000],"exons":[[100,200],[300,400]],"corrected_exons":[[101,200],[300,400]],"boundary_deletions":[[150,3]],"polya":[400,-1,-1,-1],"polya_found":true,"cage_found":false,"strand":"+","mapped_strand":"-","mapq":60,"multimapper":false,"read_group":["A"],"assignment":{A},"introns_match":true,"additional_attributes":{{"XX":"y\u00e9"}}}}"#
    )
}
fn cprof() -> String {
    format!(r#"{{"intron":{PROF},"exon":null,"split_exon":{PROF},"polya":[400,-1,-1,-1]}}"#)
}
fn pre(graph: &str, paths: &str) -> String {
    format!(r#""block":{B},"index":"0123abcd","graph":{graph},"paths":{paths}"#)
}

fn line(unit: &str, func: &str, input: &str, out: &str) -> String {
    let sep = if input.is_empty() { "" } else { "," };
    format!(r#"{{"unit":"{unit}","v":1,"case":"S:sl:8:3:n0","in":{{"fn":"{func}"{sep}{input}}},"out":{out}}}"#)
}

fn check(unit: &str, func: &str, input: &str, out: &str) {
    let l = line(unit, func, input, out);
    match reencode_case_line(&l) {
        Ok(s) => assert_eq!(s, l, "{unit}/{func} re-encoding differs"),
        Err(e) => panic!("{unit}/{func}: {e}\n{l}"),
    }
}

#[test]
fn every_payload_roundtrips() {
    let s = store();
    let rr = read_record();
    let cp = cprof();
    let pre_full = pre(r#""S:sl:8:3:n0""#, r#""S:sl:8:3:IntronPathStorage.fill:n0""#);
    let pre_null = pre("null", "null");
    let mut n = 0;
    let mut c = |u: &str, f: &str, i: &str, o: &str| {
        check(u, f, i, o);
        n += 1;
    };

    c("U1", "dump_annotation", r#""db":"x.db","gff3":"a.gff3","chr":"8","gene_types":"isoquant""#,
      r#"{"chr":"8","genes":[{"id":"gene:G1","start":1,"end":900,"strand":"-","source":"ensembl","attributes":[["ID",["gene:G1"]],["Name",["a","b"]]],"attrs":"Name \"a\"; ","transcripts":[{"id":"transcript:T1","type":"mRNA","start":1,"end":900,"strand":"-","source":"ensembl","attributes":[],"attrs":"","exons":[[1,10],[20,900]],"other_features":[[5,10,"CDS"]],"children_order_index":0}]}]}"#);
    c("U1", "get_gene_info_for_region", r#""chr":"8","region":[99,1000],"delta":4,"needs_reference":true"#,
      r#"{"gene_ids":["gene:G1"],"index":"0123abcd","empty":false}"#);
    c("U1", "from_models", &format!(r#""models":[{M}],"delta":5,"call_site":"assign_reads_to_models""#), r#"{"index":"0123abcd"}"#);
    c("U2", "fetch", r#""bam":"s.bam","region":"8:1001-2000","n":3"#,
      r#"{"qname":"r0","flag":0,"ref_id":0,"ref_name":"8","pos0":1000,"ref_end0":1100,"mapq":60,"cigar":[[0,100]],"seq":"ACGT","qual_len":4,"tags":{"HP":{"type":"Z","value":"A"},"XI":{"type":"i","value":2}}}"#);
    c("U3", "split_coverage_regions", r#""chr":"8","region":[99,5000],"read_count":10,"coverage":[[0,3],[1,0],[19,7]]"#,
      r#"{"splits":[[[99,512],3],[[4865,5000],7]]}"#);
    c("U3", "process_alignments_in_region", r#""chr":"8","region":[99,5000],"max_coverage":7,"chr_len":59000000,"small_chr":false,"coverage_cutoff":-1"#,
      r#"{"block_no":3,"skip_read_fraction":1,"gene_ids":["gene:G1"],"index":"0123abcd","genic":true}"#);
    c("U3", "cluster", r#""chr":"8","region":[99,5000],"read_count":10"#, r#"{"n_blocks":2,"first_block_no":3}"#);
    c("U3", "genes_in_region", r#""chr":"8","region":[99,5000]"#, r#"{"gene_ids":[]}"#);
    c("U4", "process_alignment",
      &format!(r#""read":{{"qname":"r0","aln_no":4}},"aln":{ALN},"polya_window":16,"polya_fraction":0.75,"max_fake_terminal_exon_len":20,"external_strand":null,"external_no_tail":false"#),
      r#"{"init":{"read_exons":[[100,200],[300,400]],"read_blocks":[[5,105],[106,206]],"cigar_blocks":[[1,1],[3,3]]},"polya":[-1,-1,-1,-1],"read_exons":[[100,200],[300,400]],"read_blocks":[[5,105],[106,206]],"cigar_blocks":[[1,1],[3,3]],"exons_changed":false,"polya_called":true}"#);
    c("U4", "process_alignment",
      &format!(r#""read":{{"qname":"r0","aln_no":4}},"aln":{ALN},"polya_window":16,"polya_fraction":0.75,"max_fake_terminal_exon_len":20,"external_strand":"+","external_no_tail":false"#),
      r#"{"init":{"read_exons":[],"read_blocks":[],"cigar_blocks":[]},"polya":null,"read_exons":[],"read_blocks":[],"cigar_blocks":[],"exons_changed":false,"polya_called":false}"#);
    c("U4", "extract_relevant_deletions", r#""read":{"qname":"r0","aln_no":4},"cigar":[[0,10],[2,3],[0,10]],"ref_start":100,"exons":[[101,123]],"margin":10"#,
      r#"{"deletions":[[111,3]]}"#);
    c("U5", "construct_profiles", &format!(r#""call_site":"collect","read":{{"qname":"r0","aln_no":4}},"index":"0123abcd","params":{P},"blocks":[[100,200],[300,400]],"polya":[400,-1,-1,-1]"#),
      &format!(r#"{{"intron":{PROF},"exon":null,"split_exon":{PROF}}}"#));
    c("U5", "construct_profiles", &format!(r#""call_site":"similar","read":null,"index":"0123abcd","params":{P},"blocks":[[100,200]],"polya":[-1,-1,-1,-1]"#),
      &format!(r#"{{"intron":{PROF},"exon":null,"split_exon":{PROF}}}"#));
    c("U5", "match_genomic_features", r#""call_site":"exon_corrector","read":{"read_idx":3},"known_features":[[1,5]],"delta":4,"read_features":[[2,5]]"#,
      r#"{"features":[[1,5]]}"#);
    c("U6", "assign_to_isoform", &format!(r#""call_site":"fl","read":{{"id":"transcript3.8.nic"}},"read_id":"transcript3.8.nic","index":"0123abcd","params":{P},"quick_mode":false,"profile":{cp}"#),
      &format!(r#"{{"assignment":{A}}}"#));
    c("U6", "assign_to_isoform", &format!(r#""call_site":"reassign","read":{{"read_idx":7}},"read_id":"r7","index":"0123abcd","params":{P},"quick_mode":true,"profile":{cp}"#),
      r#"{"assignment":null}"#);
    c("U7", "correct_assigned_read",
      &format!(r#""read":{{"qname":"r0","aln_no":4}},"index":"0123abcd","params":{P},"aln":{ALN},"read_exons":[[100,200]],"read_start":100,"read_end":200,"intron_profile":{PROF},"assignment":{A},"ref":[{{"start0":99,"seq":"ACGTacgt"}}]"#),
      r#"{"corrected_exons":[[100,200]]}"#);
    c("U7", "get_assignment_strand", &format!(r#""read":{{"qname":"r0","aln_no":4}},"assignment":{A},"polya":[-1,-1,-1,-1],"exons":[[1,5],[10,20]],"corrected_introns":[[6,9]],"intron_sites":[[[6,9],"GT","AG"]]"#),
      r#"{"strand":"-"}"#);
    c("U7", "check_antisense", &format!(r#""read":{{"qname":"r0","aln_no":4}},"strand":"+","assignment":{A}"#), &format!(r#"{{"assignment":{A}}}"#));
    c("U8", "collect_block", r#""chr":"8","block_no":3,"region":[99,5000],"index":"0123abcd","genic":true,"skip_read_fraction":1,"max_coverage":7"#,
      &format!(r#"{{"reads":[{rr}]}}"#));
    c("U9", "resolve", r#""strategy":2,"assignments":[{"assignment_id":1,"read_id":"r0","chr":"8","start":100,"end":400,"region":[99,5000],"multimapper":true,"polya_found":false,"atype":10,"gene_atype":10,"penalty":0.0,"genes":["gene:G1"],"isoforms":["transcript:T1","transcript:T2"]}]"#,
      r#"{"assignments":[]}"#);
    c("U10a", "detect_peaks_batch", r#""model":"polya","histograms":[[[400,3],[402,1]],null]"#,
      r#"{"records":[{"start":400,"mode":400,"var":0.75,"skew":"NaN","entropy":0.5623351446188083,"mean_height":0.17391304347826086,"peak_count":1,"histogram":[0,3,0,1,0],"raw_peaks":{"peak_location":[10],"peak_width":[2.04],"peak_left":[0],"peak_right":[1],"peak_heights":[3]}},null],"rows":[{"hist":0,"features":[0.75,"NaN",1.0,2.04,0.56,0.17,3.0,1.0],"proba":0.9123,"margin":2.3456,"label":1}],"peaks":[[[400,3,399,401]],[]]}"#);
    c("U10a", "xgb_predict", r#""model":"tss","caller":"_predict_rows","features":[[1.0,"Infinity",-1e-05,1e+16,0.0,0.0,0.0,0.0]]"#,
      r#"{"labels":[0],"proba":[0.1],"margin":[-2.1972246]}"#);
    c("U10b", "flush", &format!(r#""kind":"polya","block":{B},"transcripts":[["transcript:T1",{{"chr":"8","gene_id":"gene:G1","data":[400,400,402],"annotated":400}}]]"#),
      r#"{"predictions":[["transcript:T1",400]]}"#);
    c("U11", "IntronGraph.__init__", &format!(r#""block":{B},"index":"0123abcd","polya_predictions":[["transcript:T1",400]],"tss_predictions":null"#),
      &format!(r#"{{"graph":{G}}}"#));
    c("U12", "IntronPathStorage.fill", &format!(r#""block":{B},"graph":"S:sl:8:3:n0""#),
      r#"{"paths":[[[[-21,100],[201,299],[-10,400]],4]],"fl_paths":[[[-21,100],[201,299],[-10,400]]],"paths_to_reads":[[[[-21,100],[201,299],[-10,400]],[0,1,2,3]]]}"#);
    c("U12", "get_known_spliced_isoforms", &format!(r#""block":{B},"graph":"S:sl:8:3:n0","index":"0123abcd""#),
      r#"{"known":[[[[201,299]],"transcript:T1"]]}"#);
    let store_in = format!(r#"{pre_full},"store":{s}"#);
    let store_reads_in = format!(r#"{pre_full},"store":{s},"reads":[0,1,2]"#);
    let store_out = format!(r#"{{"store":{s}}}"#);
    for (u, f, i) in [
        ("U13", "construct_fl_isoforms", &store_in),
        ("U14", "construct_assignment_based_isoforms", &store_reads_in),
        ("U15a", "pre_filter_transcripts", &store_in),
        ("U15a", "filter_transcripts", &store_in),
        ("U15b", "_add_known_alternative_end_models", &store_in),
        ("U15b", "_drop_duplicate_alt_end_models", &store_in),
        ("U16", "correct_transcript_splice_sites", &store_in),
        ("U17", "assign_reads_to_models", &store_reads_in),
        ("U17", "resolve_read_ownership_conflicts", &store_in),
    ] {
        c(u, f, i, &store_out);
    }
    c("U15a", "detect_similar_isoforms", &format!(r#"{pre_null},"models":[{M}],"remove_unconfirmed":true"#),
      r#"{"to_substitute":[["transcript3.8.nic","transcript:T1"]]}"#);
    c("U15b", "correct_novel_transcript_ends", &format!(r#"{pre_full},"model":{M},"reads":[0,2]"#), &format!(r#"{{"model":{M}}}"#));
    c("U15b", "_mark_terminal_confirmation", &format!(r#"{pre_full},"model":{M},"reads":[0,2]"#), &format!(r#"{{"model":{M}}}"#));
    c("U15b", "derive_alternative_end_models", &format!(r#"{pre_full},"model":{M},"reads":[0,2],"id_counter":7"#),
      &format!(r#"{{"result":[[{M},[0]]],"id_counter":8}}"#));
    c("U16", "_correct_one_model_splice_sites",
      &format!(r#"{pre_null},"model":{M},"reads":[{{"idx":0,"start":100,"end":400,"boundary_deletions":[[199,4]]}}],"ref":[{{"start0":135,"seq":"ACGT"}}]"#),
      r#"{"exons":null}"#);
    c("U18", "join_transcripts", &format!(r#""block":{B},"index":"0123abcd","models":[{M}]"#), &format!(r#"{{"models":[{M}]}}"#));
    c("U19", "build_model_read_assignments", &store_reads_in, &format!(r#"{{"assignments":[{{"read_idx":0,"assignment":{A}}}]}}"#));
    c("U20", "gff_dump", &format!(r#""block":{B},"printer":"transcript_models","index":"0123abcd","models":[{M}],"printed_gene_ids":[],"exon_id_counter":0"#),
      r#"{"lines":["8\tIsoQuant\tgene\t100\t400\t.\t+\t.\tgene_id \"novel_gene_8_1\"; transcripts \"1\"; \n"],"exon_id_counter":2}"#);
    c("U20", "counter_add_read_info", &format!(r#""counter":{{"name":"OUT.discovered_transcript_counts","grouped":false,"group_index":-1}},"read":{{"read_idx":3}},"assignment":{A},"read_group":["A"]"#),
      r#"{"features":[["transcript:T1",[["NA",0.0,1.0]]]],"ambiguous":[0,1],"not_assigned":[5,5],"confirmed_added":[]}"#);
    c("U20", "counter_add_read_info", &format!(r#""counter":{{"name":"x","grouped":true,"group_index":0}},"read":null,"assignment":{A},"read_group":[]"#),
      r#"{"features":[],"ambiguous":[0,0],"not_assigned":[0,0],"confirmed_added":["transcript:T1"]}"#);
    c("U20", "counter_dump", r#""counter":{"name":"x","grouped":false,"group_index":-1},"feature_counter":[["transcript:T1",[["NA",2.0]]]],"all_features":["transcript:T1"],"confirmed_features":[],"ambiguous":1,"not_assigned":2,"not_aligned":0"#,
      r#"{"files":[["x.tsv","feature_id\tcount\ntranscript:T1\t2.00\n"]]}"#);
    // exception out
    c("U6", "assign_to_isoform", &format!(r#""call_site":"collect","read":{{"read_idx":7}},"read_id":"r7","index":"0123abcd","params":{P},"quick_mode":false,"profile":{cp}"#),
      r#"{"exception":"KeyError: 'x'"}"#);

    // every (unit, fn) of the payload table has been exercised
    assert!(n >= PAYLOAD_TABLE.len(), "{n} cases for {} payload entries", PAYLOAD_TABLE.len());
}

#[test]
fn side_tables_roundtrip() {
    let header = r#"{"kind":"block_header","v":1,"sample":"S","slice":"sl","chr":"8","block_no":3,"region":[100,900],"index":"0123abcd","gene_ids":["gene:G1"],"n_reads":1,"collect_region":[99,5000],"skip_read_fraction":1}"#;
    assert_eq!(reencode_block_line(header, true).unwrap(), header);
    let rr = read_record();
    assert_eq!(reencode_block_line(&rr, false).unwrap(), rr);

    let idx = index_json();
    assert_eq!(reencode_index(&idx).unwrap(), idx);

    let params = "{\"apa_delta\":50,\"data_type\":\"nanopore\",\"delta\":6,\"min_mapq\":null,\"polya_fraction\":0.75,\"resolve_ambiguous\":20,\"resolve_ambiguous__name\":\"monoexon_and_fsm\"}\n";
    assert_eq!(reencode_params(params).unwrap(), params);
}

fn index_json() -> String {
    r#"{"kind":"genes","chr":"8","start":1,"end":900,"delta":4,"genes":["gene:G1"],"gene_regions":[["gene:G1",[1,900]]],"isoform_exons":[["transcript:T2",[[1,10],[20,900]]],["transcript:T1",[[1,10],[30,900]]]],"isoform_introns":[["transcript:T2",[[11,19]]],["transcript:T1",[[11,29]]]],"isoform_strands":[["transcript:T1","-"],["transcript:T2","-"]],"gene_strands":[["gene:G1","-"]],"gene_id_map":[["transcript:T1","gene:G1"],["transcript:T2","gene:G1"]],"intron_features":[[11,19],[11,29]],"intron_profiles":[["transcript:T2",[1,-1]],["transcript:T1",[-1,1]]],"intron_ranges":[["transcript:T2",[0,1]],["transcript:T1",[1,2]]],"exon_features":[],"exon_profiles":[],"exon_ranges":[],"split_exon_features":[[1,10],[20,29],[30,900]],"split_exon_profiles":[["transcript:T2",[1,1,1]],["transcript:T1",[1,-1,1]]],"split_exon_ranges":[["transcript:T2",[0,3]],["transcript:T1",[0,3]]],"sources":[["gene:G1","ensembl"]],"other_features":[["transcript:T1",[[5,10,"CDS"]]]],"feature_attributes":[["gene:G1","Name \"a\"; "]]}"#.to_owned()
}

#[test]
fn index_from_dump() {
    let d: IndexDump = serde_json::from_str(&index_json()).unwrap();
    let mut ids = IdMaps::new();
    let idx = IsoformIndex::from_dump(&d, &mut ids).unwrap();
    assert_eq!(idx.n_iso(), 2);
    assert_eq!(ids.iso_name(idx.iso[0].key), "transcript:T2");
    assert_eq!(idx.iso[0].strand, Strand::Minus);
    assert_eq!(idx.iso_introns(1), &[Iv::new(11, 29)]);
    assert_eq!(idx.intron_prof.row(1), &[-1, 1]);
    assert_eq!(idx.split_prof.row(1), &[1, -1, 1]);
    assert_eq!(idx.prof_range[1][0], (1, 2));
    assert_eq!(idx.id_rank, vec![1, 0]); // T1 < T2
    assert_eq!(idx.genes.len(), 1);
    assert_eq!(idx.overlapping_isoforms(Iv::new(25, 25)), vec![0, 1]);
}

#[test]
fn read_record_core_roundtrip() {
    let f: ReadRecordFixture = serde_json::from_str(&read_record()).unwrap();
    let mut ids = IdMaps::new();
    let (core, ex) = f.to_core(&mut ids);
    assert_eq!(core.assignment.matches[0].tx, Some(IsoRef::Ref(0)));
    assert!(core.group > 0);
    // core → BlockStore → core → fixture
    let st = BlockStore::from_records(std::slice::from_ref(&core));
    let back = ReadRecordFixture::from_core(&st.record(0), &ids, &ex);
    assert_eq!(back, f);
    assert_eq!(fixture::to_py_string(&back).unwrap(), read_record());
}

#[test]
fn store_and_model_core_roundtrip() {
    let sf: StoreFixture = serde_json::from_str(&store()).unwrap();
    let mut ids = IdMaps::new();
    let names = ["r0", "r1", "r2"];
    let (st, ex) = sf.to_core(&mut ids, &names).unwrap();
    assert_eq!(st.read_assignment_counts, vec![1, 0, 2]);
    assert_eq!(st.construction_assignment[0], Some(0));
    let back = StoreFixture::from_core(&st, &ids, &ex, &names);
    assert_eq!(back.models, sf.models);
    assert!(fixture::store_diff(&back, &sf, &Default::default()).is_none());

    let g: GraphFixture = serde_json::from_str(G).unwrap();
    let core = IntronGraph::from(&g);
    assert_eq!(GraphFixture::from(&core), g);

    let a: AssignmentFixture = serde_json::from_str(A).unwrap();
    let ca = a.to_core(&mut ids, TxSpace::Model);
    assert_eq!(AssignmentFixture::from_core(&ca, &ids), a);
}

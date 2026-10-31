"""Encoders: IsoQuant objects -> JSON-ready (deep-copied) structures, per SCHEMA.md."""
import hashlib
from enum import Enum

from .sink import fval, dumps


def iv(x):
    return [int(x[0]), int(x[1])]


def ivs(xs):
    return [[int(a), int(b)] for a, b in (xs or [])]


def ev(v):
    return v.value if isinstance(v, Enum) else v


def scalar(v):
    if isinstance(v, Enum):
        return v.value
    if isinstance(v, bool) or v is None or isinstance(v, str):
        return v
    if isinstance(v, int):
        return int(v)
    if isinstance(v, float):
        return fval(v)
    return v


def vertex(v):
    return [int(v[0]), int(v[1])]


def path(p):
    return [vertex(v) for v in p]


def polya(pi):
    if pi is None:
        return None
    return [int(pi.external_polya_pos), int(pi.external_polyt_pos),
            int(pi.internal_polya_pos), int(pi.internal_polyt_pos)]


def event(e):
    return {"t": int(e.event_type.value), "iso": [int(e.isoform_region[0]), int(e.isoform_region[1])],
            "read": [int(e.read_region[0]), int(e.read_region[1])], "info": int(e.event_info)}


def match(m):
    return {"gene": m.assigned_gene, "tx": m.assigned_transcript, "strand": m.transcript_strand,
            "class": int(m.match_classification.value), "penalty": fval(m.penalty_score),
            "events": [event(e) for e in m.match_subclassifications]}


def assignment(ra):
    if ra is None:
        return None
    return {"atype": int(ra.assignment_type.value), "gene_atype": int(ra.gene_assignment_type.value),
            "matches": [match(m) for m in (ra.isoform_matches or [])]}


def _read_group(ra):
    try:
        return [str(g) for g in ra.read_group]
    except Exception:
        return []


def read_record(ra, idx):
    return {
        "idx": idx,
        "assignment_id": int(ra.assignment_id),
        "name": ra.read_id,
        "chr": ra.chr_id,
        "region": [int(ra.genomic_region[0]), int(ra.genomic_region[1])],
        "exons": ivs(ra.exons),
        "corrected_exons": ivs(ra.corrected_exons),
        "boundary_deletions": ivs(ra.boundary_deletions),
        "polya": polya(ra.polya_info),
        "polya_found": bool(ra.polyA_found),
        "cage_found": bool(ra.cage_found),
        "strand": ra.strand,
        "mapped_strand": ra.mapped_strand,
        "mapq": int(ra.mapping_quality),
        "multimapper": bool(ra.multimapper),
        "read_group": _read_group(ra),
        "assignment": assignment(ra),
        "introns_match": bool(ra.introns_match),
        "additional_attributes": {str(k): str(v) for k, v in (ra.additional_attributes or {}).items()},
    }


def basic(a):
    return {"assignment_id": int(a.assignment_id), "read_id": a.read_id, "chr": a.chr_id,
            "start": int(a.start), "end": int(a.end), "region": [int(a.genomic_region[0]), int(a.genomic_region[1])],
            "multimapper": bool(a.multimapper), "polya_found": bool(a.polyA_found),
            "atype": int(a.assignment_type.value), "gene_atype": int(a.gene_assignment_type.value),
            "penalty": fval(a.penalty_score), "genes": list(a.genes), "isoforms": list(a.isoforms)}


def prof(p):
    if p is None:
        return None
    return {"gene_profile": [int(x) for x in p.gene_profile], "read_profile": [int(x) for x in p.read_profile],
            "read_features": ivs(p.read_features), "gene_profile_range": [int(p.gene_profile_range[0]), int(p.gene_profile_range[1])]}


def cprof(cp, with_polya=True):
    d = {"intron": prof(cp.read_intron_profile), "exon": prof(cp.read_exon_profile),
         "split_exon": prof(cp.read_split_exon_profile)}
    if with_polya:
        d["polya"] = polya(cp.polya_info)
    return d


def aln(a):
    return {"qname": a.query_name, "flag": int(a.flag), "pos0": int(a.reference_start),
            "ref_end0": int(a.reference_end) if a.reference_end is not None else None,
            "mapq": int(a.mapping_quality), "is_reverse": bool(a.is_reverse),
            "cigar": [[int(o), int(l)] for o, l in (a.cigartuples or [])], "seq": a.query_sequence}


ASSIGN_PARAMS = ["delta", "minor_exon_extension", "major_exon_extension", "min_abs_exon_overlap",
                 "min_rel_exon_overlap", "minimal_exon_overlap", "minimal_intron_absence_overlap",
                 "micro_intron_length", "max_fake_terminal_exon_len", "max_intron_shift", "max_missed_exon_len",
                 "max_intron_abs_diff", "max_intron_rel_diff", "max_suspicious_intron_abs_len",
                 "max_suspicious_intron_rel_len", "resolve_ambiguous", "apa_delta", "count_exons",
                 "correct_fuzzy_junctions", "correct_intron_shifts", "correct_skipped_exons",
                 "correct_fake_terminal_exons", "correct_terminal_exons", "correct_microintron_retention"]


def assign_params(params):
    return {k: scalar(getattr(params, k, None)) for k in ASSIGN_PARAMS}


def model(m):
    return {"id": m.transcript_id, "gene": m.gene_id, "chr": m.chr_id, "strand": m.strand,
            "type": int(m.transcript_type.value), "exons": ivs(m.exon_blocks),
            "intron_path": path(m.intron_path or ()),
            "polya_confirmed": bool(getattr(m, "polya_confirmed", False)),
            "tss_confirmed": bool(getattr(m, "tss_confirmed", False)),
            "source": m.source,
            "other_features": [[int(a), int(b), str(t)] for a, b, t in (m.other_features or [])],
            "additional_info": [[str(k), str(v)] for k, v in m.additional_info.items()]}


def ridx(a):
    return getattr(a, "_iq_idx", -1)


def store(st, ctx, name_idx=None):
    from isoquant_lib.model_construction.context import ModelStore
    idd = getattr(ctx, "id_distributor", None) if ctx is not None else None
    if name_idx is None:
        from . import hooks
        name_idx = hooks.CTX.name_idx if hooks.CTX is not None else {}
    return {
        "models": [model(m) for m in st.transcript_model_storage],
        "transcript_read_ids": [[k, [ridx(a) for a in v]] for k, v in st.transcript_read_ids.items()],
        "internal_counter": [[k, int(v)] for k, v in st.internal_counter.items()],
        "read_assignment_counts": [[k, int(v)] for k, v in st.read_assignment_counts.items()],
        "construction_assignment": [[k, v] for k, v in st.construction_assignment.items()],
        "read_assignment_counts_idx": [[i, int(v)] for k, v in st.read_assignment_counts.items()
                                       for i in name_idx.get(k, [])],
        "construction_assignment_idx": [[i, v] for k, v in st.construction_assignment.items()
                                        for i in name_idx.get(k, [])],
        "nonfl_known_reads": [[k, [ridx(a) for a in v]] for k, v in st.nonfl_known_reads.items()],
        "detected_known_isoforms": sorted(ModelStore.detected_known_isoforms),
        "id_counter": int(idd.value) if idd is not None else None,
    }


def _pairs(d, f=lambda v: v):
    return [[k, f(v)] for k, v in d.items()]


def _profiles(fp):
    return ([iv(x) for x in fp.features],
            [[k, [int(x) for x in v]] for k, v in fp.profiles.items()],
            [[k, [int(v[0]), int(v[1])]] for k, v in fp.profile_ranges.items()])


def gene_info_dump(gi):
    """IsoformIndex dump (SCHEMA §2.1). Returns (bytes, sha1)."""
    if getattr(gi, "db", None) is not None and gi.gene_db_list:
        kind = "genes"
    elif hasattr(gi, "_iq_kind"):
        kind = gi._iq_kind
    else:
        kind = "region"
    genes = [g.id for g in gi.gene_db_list] if kind == "genes" else []
    gene_regions = [[g.id, [int(g.start), int(g.end)]] for g in gi.gene_db_list] if kind == "genes" else []
    i_f, i_p, i_r = _profiles(gi.intron_profiles)
    e_f, e_p, e_r = _profiles(gi.exon_profiles)
    s_f, s_p, s_r = _profiles(gi.split_exon_profiles)
    d = {
        "kind": kind,
        "chr": gi.chr_id,
        "start": int(gi.start), "end": int(gi.end), "delta": int(gi.delta),
        "genes": genes,
        "gene_regions": gene_regions,
        "isoform_exons": _pairs(gi.all_isoforms_exons, ivs),
        "isoform_introns": _pairs(gi.all_isoforms_introns, ivs),
        "isoform_strands": _pairs(gi.isoform_strands),
        "gene_strands": _pairs(gi.gene_strands),
        "gene_id_map": _pairs(gi.gene_id_map),
        "intron_features": i_f, "intron_profiles": i_p, "intron_ranges": i_r,
        "exon_features": e_f, "exon_profiles": e_p, "exon_ranges": e_r,
        "split_exon_features": s_f, "split_exon_profiles": s_p, "split_exon_ranges": s_r,
        "sources": _pairs(gi.sources),
        "other_features": [[k, [[int(a), int(b), str(t)] for a, b, t in v]] for k, v in gi.other_features.items()],
        "feature_attributes": _pairs(gi.feature_attributes, str) if kind == "genes" else [],
    }
    b = dumps(d).encode("utf-8")
    return b, hashlib.sha1(b).hexdigest()

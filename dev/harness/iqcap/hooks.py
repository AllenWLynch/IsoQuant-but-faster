"""Recording monkeypatches for IsoQuant 4.1.0 (design 07 §10.2 patch table).

All wrappers encode their inputs (a JSON-ready deep copy) BEFORE calling the original and never
change arguments, return values or IsoQuant state, except for the private `_iq_*` attributes they
attach to objects (read identity, cached index hashes).
"""
import functools
import logging
import os
import sys
import traceback
from collections import defaultdict

from . import enc
from .sink import keep, fval, BlockWriter

log = logging.getLogger("iqcap")

ASSIGN_CALL_SITES = {
    "process_genic": "collect", "construct_fl_isoforms": "fl", "detect_similar_isoforms": "similar",
    "assign_reads_to_models": "reassign", "_model_read_assignment_single": "final_single",
    "build_model_read_assignments": "final_full", "compare_models_with_known": "sqanti",
    "correct_assigned_read": "exon_corrector", "process_events": "exon_corrector",
}
CONSTRUCT_UNITS = {"U8", "U10b", "U11", "U12", "U13", "U14", "U15a", "U15b", "U16", "U17", "U18", "U19"}


class Ctx:
    def __init__(self, sinks, sample, slice_name, out_dir):
        self.S = sinks
        self.sample = sample
        self.slice = slice_name
        self.out_dir = out_dir
        self.phase = "setup"
        self.chr = None
        # collect
        self.cblock = defaultdict(int)  # chr -> next collect block no
        self.block = None
        self.block_info = {}
        self.aln_no = -1
        self.cur_aln = None
        self.cur_pos0 = None  # reference_start of the current collect alignment (U4 extract_relevant_deletions)
        self.collect_regions = {}  # (chr, block) -> dict
        # construct
        self.kblock = defaultdict(int)
        self.block_selected = False
        self.graph_key = None
        self.paths_key = None
        self.block_index = None
        self.counters = defaultdict(int)
        self.guard = 0
        self.u8_pending = None
        self.blocks_meta = []
        self.skipped_blocks = []
        self.cluster_first_block = None
        self.name_idx = {}

    # ---- keys / sampling ----
    def bkey(self, chr_id=None, block=None):
        c = chr_id if chr_id is not None else (self.chr or "-")
        b = block if block is not None else (self.block if self.block is not None else "-")
        return "%s:%s:%s:%s" % (self.sample, self.slice, c, b)

    def local(self, unit, fn, base):
        k = (unit, self.bkey(), fn, base)
        n = self.counters[k]
        self.counters[k] += 1
        return "%s:%s.%d" % (fn, base, n)

    def selected(self, unit, sample_key):
        if not self.S.wants(unit) or self.guard:
            return False
        if self.phase == "construct" or unit in CONSTRUCT_UNITS:
            if self.phase == "construct" and not self.block_selected:
                return False
            if unit in CONSTRUCT_UNITS:
                return True
        return keep(sample_key, self.S.rate(unit))

    def index(self, gi):
        if gi is None:
            return None
        sha = getattr(gi, "_iq_sha", None)
        if sha is None:
            b, sha = enc.gene_info_dump(gi)
            self.S.write_index(b, sha)
            try:
                gi._iq_sha = sha
            except Exception:
                pass
        return sha

    def block_ref(self):
        return {"sample": self.sample, "slice": self.slice, "chr": self.chr, "block_no": self.block}

    def read_ident(self, ra=None):
        if self.phase == "collect" and self.cur_aln is not None:
            return {"qname": self.cur_aln[0], "aln_no": self.cur_aln[1]}
        if ra is not None and hasattr(ra, "_iq_idx"):
            return {"read_idx": ra._iq_idx}
        return None


CTX = None


def _caller_site(depth=2, table=ASSIGN_CALL_SITES):
    f = sys._getframe(depth)
    first = f.f_code.co_name
    while f is not None:
        name = f.f_code.co_name
        if name in table:
            return table[name], f
        f = f.f_back
    return first, None


def _find_local_read(frame):
    f = frame
    for _ in range(6):
        if f is None:
            break
        for nm in ("source", "assignment", "read_assignment", "a"):
            v = f.f_locals.get(nm)
            if v is not None and hasattr(v, "_iq_idx"):
                return v
        f = f.f_back
    return None


def _write(unit, case_local, inp, out):
    CTX.S.write(unit, CTX.bkey() + ":" + case_local, inp, out)


def _exc_out(e):
    return {"exception": "%s: %s" % (type(e).__name__, e)}


def patch(obj, name, factory):
    orig = getattr(obj, name)
    w = factory(orig)
    try:
        functools.update_wrapper(w, orig)
    except Exception:
        pass
    setattr(obj, name, w)
    return orig


# ============================================================================ install
def install(ctx):
    global CTX
    CTX = ctx
    from . import ensure_isoquant_on_path
    ensure_isoquant_on_path()
    _install_identity()
    _install_u1()
    _install_u3()
    _install_u4()
    _install_u5_u6()
    _install_u7()
    _install_u8()
    _install_u9()
    _install_u10()
    _install_u11_u12()
    _install_construct()
    _install_u20()


# ============================================================================ identity / phases
def _install_identity():
    from isoquant_lib.alignment.alignment_processor import AlignmentCollector
    from isoquant_lib.assignment import assignment_loader as AL
    from isoquant_lib.gene_info import GeneInfo

    def f_init(orig):
        def w(self, chr_id, *a, **kw):
            CTX.phase = "collect"
            CTX.chr = chr_id
            CTX.block = None
            return orig(self, chr_id, *a, **kw)
        return w
    patch(AlignmentCollector, "__init__", f_init)

    def f_loader_init(orig):
        def w(self, *a, **kw):
            CTX.phase = "construct"
            CTX.block = None
            CTX.block_selected = False
            return orig(self, *a, **kw)
        return w
    patch(AL.ReadAssignmentLoader, "__init__", f_loader_init)

    def f_get_next(orig):
        def w(self):
            gene_info, storage = orig(self)
            if gene_info is None and storage is None:
                return gene_info, storage
            chr_id = gene_info.chr_id if gene_info is not None else CTX.chr
            CTX.chr = chr_id
            b = CTX.kblock[chr_id]
            CTX.kblock[chr_id] += 1
            CTX.block = b
            CTX.graph_key = None
            CTX.paths_key = None
            CTX.name_idx = {}
            for i, ra in enumerate(storage):
                ra._iq_idx = i
                CTX.name_idx.setdefault(ra.read_id, []).append(i)
            n = len(storage)
            sel = keep(CTX.bkey(), CTX.S.block_rate)
            if sel and n > CTX.S.block_cap:
                sel = False
                CTX.skipped_blocks.append({"chr": chr_id, "block_no": b, "reason": "block_cap", "n_reads": n})
            CTX.block_selected = sel and bool(CTX.S.units & (CONSTRUCT_UNITS | {"U20", "U5", "U6", "U1", "U10a"}))
            meta = {"chr": chr_id, "block_no": b, "region": [int(gene_info.start), int(gene_info.end)],
                    "n_reads": n, "captured": CTX.block_selected}
            CTX.blocks_meta.append(meta)
            if CTX.block_selected:
                sha = CTX.index(gene_info)
                CTX.block_index = sha
                cr = CTX.collect_regions.get((chr_id, b), {})
                bw = BlockWriter(CTX.out_dir, CTX.sample, chr_id, b)
                bw.write({"kind": "block_header", "v": 1, "sample": CTX.sample, "slice": CTX.slice,
                          "chr": chr_id, "block_no": b, "region": meta["region"], "index": sha,
                          "gene_ids": [g.id for g in gene_info.gene_db_list], "n_reads": n,
                          "collect_region": cr.get("region"), "skip_read_fraction": cr.get("skip_read_fraction")})
                for i, ra in enumerate(storage):
                    bw.write(enc.read_record(ra, i))
                bw.close()
            return gene_info, storage
        return w
    patch(AL.ReadAssignmentLoader, "get_next", f_get_next)

    # kind tags for index dumps
    for meth, kind in (("from_models", "models"), ("from_model", "model"), ("from_region", "region")):
        orig_cm = getattr(GeneInfo, meth).__func__

        def mk(orig_cm, kind, meth):
            def w(cls, *a, **kw):
                inp = None
                if meth == "from_models" and CTX.phase == "construct" and CTX.selected("U1", CTX.bkey()):
                    site, _ = _caller_site(2, {})
                    models = a[0] if a else kw.get("transcript_model_storage")
                    delta = a[1] if len(a) > 1 else kw.get("delta", 0)
                    inp = {"fn": "from_models", "models": [enc.model(m) for m in models], "delta": int(delta),
                           "call_site": sys._getframe(1).f_code.co_name}
                gi = orig_cm(cls, *a, **kw)
                gi._iq_kind = kind
                if inp is not None:
                    _write("U1", CTX.local("U1", "from_models", "n"), inp, {"index": CTX.index(gi)})
                return gi
            return w
        setattr(GeneInfo, meth, classmethod(mk(orig_cm, kind, meth)))


# ============================================================================ U1 (b)
def _install_u1():
    from isoquant_lib.alignment.alignment_processor import AlignmentCollector

    def f(orig):
        def w(self, current_region):
            gi = orig(self, current_region)
            if CTX.selected("U1", CTX.bkey()):
                inp = {"fn": "get_gene_info_for_region", "chr": self.chr_id,
                       "region": [int(current_region[0]), int(current_region[1])],
                       "delta": int(self.params.delta), "needs_reference": bool(self.params.needs_reference)}
                out = {"gene_ids": [g.id for g in gi.gene_db_list], "index": CTX.index(gi), "empty": bool(gi.empty())}
                _write("U1", CTX.local("U1", "get_gene_info_for_region", "n"), inp, out)
            CTX.block_info["gene_info"] = gi
            return gi
        return w
    patch(AlignmentCollector, "get_gene_info_for_region", f)


# ============================================================================ U3
def _install_u3():
    from isoquant_lib.alignment.alignment_processor import AlignmentCollector

    orig_split = AlignmentCollector.__dict__["split_coverage_regions"].__func__

    def w_split(genomic_region, alignment_storage):
        inp = None
        if CTX.S.wants("U3"):
            key_local = CTX.local("U3", "split_coverage_regions", "n")
            if keep(CTX.bkey() + ":" + key_local, CTX.S.rate("U3")):
                inp = {"fn": "split_coverage_regions", "chr": CTX.chr,
                       "region": [int(genomic_region[0]), int(genomic_region[1])],
                       "read_count": int(alignment_storage.get_read_count()),
                       "coverage": [[int(k), int(v)] for k, v in sorted(alignment_storage.coverage_dict.items())]}
        res = orig_split(genomic_region, alignment_storage)
        if inp is not None:
            _write("U3", key_local, inp, {"splits": [[[int(r[0][0]), int(r[0][1])], int(r[1])] for r in res]})
        return res
    AlignmentCollector.split_coverage_regions = staticmethod(w_split)

    def f_fwd(orig):
        def w(self, alignment_storage):
            region = alignment_storage.region
            n = alignment_storage.get_read_count()
            first = CTX.cblock[self.chr_id]
            saved_block = CTX.block
            CTX.block = None
            CTX.chr = self.chr_id
            # the cluster case is keyed outside blocks
            key_local = CTX.local("U3", "cluster", "n")
            CTX.block = saved_block
            for res in orig(self, alignment_storage):
                yield res
            if CTX.S.wants("U3") and keep(key_local, CTX.S.rate("U3")):
                CTX.S.write("U3", "%s:%s:%s:-:%s" % (CTX.sample, CTX.slice, self.chr_id, key_local),
                            {"fn": "cluster", "chr": self.chr_id, "region": [int(region[0]), int(region[1])],
                             "read_count": int(n)},
                            {"n_blocks": CTX.cblock[self.chr_id] - first, "first_block_no": first})
        return w
    patch(AlignmentCollector, "forward_alignments", f_fwd)

    def f_pair(orig):
        def w(self, current_region, alignment_storage, max_coverage):
            b = CTX.cblock[self.chr_id]
            CTX.cblock[self.chr_id] += 1
            CTX.chr = self.chr_id
            CTX.block = b
            CTX.aln_no = -1
            CTX.cur_aln = None
            CTX.block_info = {"region": [int(current_region[0]), int(current_region[1])],
                              "max_coverage": int(max_coverage)}
            res = orig(self, current_region, alignment_storage, max_coverage)
            gi = res[0]
            info = CTX.block_info
            CTX.collect_regions[(self.chr_id, b)] = {"region": info["region"],
                                                     "skip_read_fraction": info.get("skip", 1)}
            if CTX.S.wants("U3") and keep(CTX.bkey() + ":pair", CTX.S.rate("U3")):
                small = self.chr_id in self.SMALL_CHR_IDS or len(self.chr_record) < self.SMALL_CHR_LEN
                inp = {"fn": "process_alignments_in_region", "chr": self.chr_id,
                       "region": [int(current_region[0]), int(current_region[1])],
                       "max_coverage": int(max_coverage), "chr_len": len(self.chr_record), "small_chr": bool(small),
                       "coverage_cutoff": int(self.small_chr_max_coverage if small else self.usual_gene_max_coverage)}
                out = {"block_no": b, "skip_read_fraction": int(info.get("skip", 1)),
                       "gene_ids": [g.id for g in gi.gene_db_list], "index": CTX.index(gi), "genic": not gi.empty()}
                _write("U3", CTX.local("U3", "process_alignments_in_region", "n"), inp, out)
            return res
        return w
    patch(AlignmentCollector, "process_alignments_in_region", f_pair)

    for nm, genic in (("process_genic", True), ("process_intergenic", False)):
        def f_proc(orig, genic=genic):
            if genic:
                def w(self, alignment_storage, gene_info, region, skip_read_fraction=1):
                    CTX.block_info["skip"] = skip_read_fraction
                    CTX.block_info["genic"] = True
                    return orig(self, alignment_storage, gene_info, region, skip_read_fraction)
            else:
                def w(self, alignment_storage, region, skip_read_fraction=1):
                    CTX.block_info["skip"] = skip_read_fraction
                    CTX.block_info["genic"] = False
                    return orig(self, alignment_storage, region, skip_read_fraction)
            return w
        patch(AlignmentCollector, nm, f_proc)

    def f_genes(orig):
        def w(self, current_region):
            res = orig(self, current_region)
            if CTX.S.wants("U3"):
                k = "%s:%s:%s:-:" % (CTX.sample, CTX.slice, self.chr_id)
                ckey = ("U3", k, "genes_in_region")
                n = CTX.counters[ckey]
                CTX.counters[ckey] += 1
                key = k + "genes_in_region:n.%d" % n
                if keep(key, CTX.S.rate("U3")):
                    CTX.S.write("U3", key, {"fn": "genes_in_region", "chr": self.chr_id,
                                            "region": [int(current_region[0]), int(current_region[1])]},
                                {"gene_ids": sorted(res)})
            return res
        return w
    patch(AlignmentCollector, "get_genes_in_region", f_genes)


# ============================================================================ U4
def _install_u4():
    from isoquant_lib.alignment.alignment_info import AlignmentInfo
    from isoquant_lib.alignment import alignment_processor as AP

    def f_init(orig):
        def w(self, alignment):
            orig(self, alignment)
            if CTX.phase == "collect":
                CTX.aln_no += 1
                CTX.cur_aln = (alignment.query_name, CTX.aln_no)
                CTX.cur_pos0 = int(alignment.reference_start) if alignment.reference_start is not None else None
                self._iq_aln = CTX.aln_no
                if CTX.S.wants("U4") and keep(CTX.bkey() + ":a%d" % CTX.aln_no, CTX.S.rate("U4")):
                    self._iq_u4 = {"aln": enc.aln(alignment),
                                   "init": {"read_exons": enc.ivs(self.read_exons),
                                            "read_blocks": enc.ivs(self.read_blocks),
                                            "cigar_blocks": enc.ivs(self.cigar_blocks)}}
                    if not self.read_exons:
                        _write("U4", "process_alignment:a%d.0" % CTX.aln_no,
                               {"fn": "process_alignment", "read": CTX.read_ident(), "aln": self._iq_u4["aln"],
                                "polya_window": None, "polya_fraction": None, "max_fake_terminal_exon_len": None,
                                "external_strand": None, "external_no_tail": False},
                               {"init": self._iq_u4["init"], "polya": None, "read_exons": [], "read_blocks": [],
                                "cigar_blocks": [], "exons_changed": False, "polya_called": False})
        return w
    patch(AlignmentInfo, "__init__", f_init)

    def f_polya(orig):
        def w(self, polya_finder, polya_fixer, external_strand=None, external_no_tail=False):
            u4 = getattr(self, "_iq_u4", None)
            if u4 is None or CTX.guard:
                return orig(self, polya_finder, polya_fixer, external_strand=external_strand,
                            external_no_tail=external_no_tail)
            inp = {"fn": "process_alignment", "read": CTX.read_ident(), "aln": u4["aln"],
                   "polya_window": int(polya_finder.window_size), "polya_fraction": fval(polya_finder.min_polya_fraction),
                   "max_fake_terminal_exon_len": int(polya_fixer.params.max_fake_terminal_exon_len),
                   "external_strand": external_strand, "external_no_tail": bool(external_no_tail)}
            res = orig(self, polya_finder, polya_fixer, external_strand=external_strand,
                       external_no_tail=external_no_tail)
            out = {"init": u4["init"], "polya": enc.polya(self.polya_info), "read_exons": enc.ivs(self.read_exons),
                   "read_blocks": enc.ivs(self.read_blocks), "cigar_blocks": enc.ivs(self.cigar_blocks),
                   "exons_changed": bool(self.exons_changed), "polya_called": True}
            _write("U4", CTX.local("U4", "process_alignment", "a%d" % self._iq_aln), inp, out)
            del self._iq_u4
            return res
        return w
    patch(AlignmentInfo, "add_polya_info", f_polya)

    orig_erd = AP.extract_relevant_deletions

    def _erd_pos0():
        # alignment.reference_start of the alignment being processed (caller local `alignment`;
        # falls back to the value recorded by AlignmentInfo.__init__)
        try:
            al = sys._getframe(2).f_locals.get("alignment")
            if al is not None and al.reference_start is not None:
                return int(al.reference_start)
        except Exception:
            pass
        return CTX.cur_pos0

    def w_erd(cigartuples, ref_start, exons, *a, **kw):
        sel = (CTX.phase == "collect" and CTX.cur_aln is not None and CTX.S.wants("U4")
               and keep(CTX.bkey() + ":a%d" % CTX.cur_aln[1], CTX.S.rate("U4")))
        if sel:
            margin = a[0] if a else kw.get("margin")
            if margin is None:
                from isoquant_lib.model_construction.transcript_splice_site_corrector import BOUNDARY_MARGIN
                margin = BOUNDARY_MARGIN
            inp = {"fn": "extract_relevant_deletions", "read": CTX.read_ident(),
                   "cigar": [[int(o), int(l)] for o, l in (cigartuples or [])], "ref_start": int(ref_start),
                   "exons": enc.ivs(exons), "margin": int(margin), "pos0": _erd_pos0()}
        res = orig_erd(cigartuples, ref_start, exons, *a, **kw)
        if sel:
            _write("U4", CTX.local("U4", "extract_relevant_deletions", "a%d" % CTX.cur_aln[1]), inp,
                   {"deletions": enc.ivs(res)})
        return res
    AP.extract_relevant_deletions = w_erd


# ============================================================================ U5 / U6
def _read_key_and_ident(frame_read=None):
    """Sampling key, local base, read identity for per-read units."""
    if CTX.phase == "collect" and CTX.cur_aln is not None:
        base = "a%d" % CTX.cur_aln[1]
        return CTX.bkey() + ":" + base, base, {"qname": CTX.cur_aln[0], "aln_no": CTX.cur_aln[1]}
    if frame_read is not None:
        base = "r%d" % frame_read._iq_idx
        return CTX.bkey() + ":" + base, base, {"read_idx": frame_read._iq_idx}
    return None, None, None


def _install_u5_u6():
    from isoquant_lib.assignment.long_read_profiles import CombinedProfileConstructor, OverlappingFeaturesProfileConstructor
    from isoquant_lib.assignment.long_read_assigner import LongReadAssigner

    def f_cp(orig):
        def w(self, sorted_blocks, polya_info, cage_hits):
            if not CTX.S.wants("U5") or CTX.guard:
                return orig(self, sorted_blocks, polya_info, cage_hits)
            site, fr = _caller_site(2)
            rd = _find_local_read(sys._getframe(1)) if CTX.phase == "construct" else None
            skey, base, ident = _read_key_and_ident(rd)
            if ident is None:
                base = "n"
                ident = None
                skey = CTX.bkey() + ":" + site
            if not CTX.selected("U5", skey + ":" + site):
                return orig(self, sorted_blocks, polya_info, cage_hits)
            inp = {"fn": "construct_profiles", "call_site": site, "read": ident, "index": CTX.index(self.gene_info),
                   "params": enc.assign_params(self.params), "blocks": enc.ivs(sorted_blocks),
                   "polya": enc.polya(polya_info)}
            res = orig(self, sorted_blocks, polya_info, cage_hits)
            _write("U5", CTX.local("U5", "construct_profiles", base), inp, enc.cprof(res, with_polya=False))
            return res
        return w
    patch(CombinedProfileConstructor, "construct_profiles", f_cp)

    def f_mgf(orig):
        def w(self, read_features):
            if not CTX.S.wants("U5") or CTX.guard:
                return orig(self, read_features)
            skey, base, ident = _read_key_and_ident(None)
            if skey is None or not CTX.selected("U5", skey):
                return orig(self, read_features)
            site, _ = _caller_site(2)
            inp = {"fn": "match_genomic_features", "call_site": site, "read": ident,
                   "known_features": enc.ivs(self.known_features), "delta": _effective_delta(self),
                   "read_features": enc.ivs(read_features)}
            res = orig(self, read_features)
            _write("U5", CTX.local("U5", "match_genomic_features", base), inp, {"features": enc.ivs(res)})
            return res
        return w
    patch(OverlappingFeaturesProfileConstructor, "match_genomic_features", f_mgf)

    def _effective_delta(self):
        # the delta the comparator actually uses (partial(equal_ranges, delta=...)); the constructor
        # kwarg self.delta is 0 at the exon_corrector call site (see SCHEMA.md changelog)
        kw = getattr(self.comparator, "keywords", None) or {}
        return int(kw["delta"]) if "delta" in kw else int(self.delta)

    def f_assign(orig):
        def w(self, read_id, combined_read_profile):
            if not CTX.S.wants("U6") or CTX.guard:
                return orig(self, read_id, combined_read_profile)
            site, fr = _caller_site(2)
            ident = None
            if CTX.phase == "construct":
                rd = _find_local_read(sys._getframe(1))
                if site in ("fl", "similar", "sqanti") or rd is None:
                    base = "id:" + str(read_id)
                    ident = {"id": str(read_id)}
                    skey = CTX.bkey() + ":" + site + ":" + base
                else:
                    skey, base, ident = _read_key_and_ident(rd)
                    skey += ":" + site
            else:
                skey, base, ident = _read_key_and_ident(None)
                if skey is None:
                    base = "id:" + str(read_id)
                    skey = CTX.bkey() + ":" + base
            if not CTX.selected("U6", skey):
                return orig(self, read_id, combined_read_profile)
            inp = {"fn": "assign_to_isoform", "call_site": site, "read": ident, "read_id": str(read_id),
                   "index": CTX.index(self.gene_info), "params": enc.assign_params(self.params),
                   "quick_mode": bool(self.quick_mode), "profile": enc.cprof(combined_read_profile)}
            try:
                res = orig(self, read_id, combined_read_profile)
            except Exception as e:
                _write("U6", CTX.local("U6", "assign_to_isoform", base), inp, _exc_out(e))
                raise
            _write("U6", CTX.local("U6", "assign_to_isoform", base), inp, {"assignment": enc.assignment(res)})
            return res
        return w
    patch(LongReadAssigner, "assign_to_isoform", f_assign)


# ============================================================================ U7
def _refwin(chr_record, start0, end0):
    start0 = max(0, int(start0))
    end0 = min(len(chr_record), int(end0))
    return {"start0": start0, "seq": str(chr_record[start0:end0]) if end0 > start0 else ""}


U7_PAD = 32


def _exon_windows(chr_record, exons):
    """RefWins covering each read exon +-U7_PAD bp (merged when overlapping); enough for get_error_count."""
    wins = []
    for s, e in exons:
        a, b = s - 1 - U7_PAD, e + U7_PAD
        if wins and a <= wins[-1][1]:
            wins[-1][1] = max(wins[-1][1], b)
        else:
            wins.append([a, b])
    return [_refwin(chr_record, a, b) for a, b in wins]


def _install_u7():
    from isoquant_lib.assignment.exon_corrector import ExonCorrector
    from isoquant_lib.alignment.alignment_processor import AlignmentCollector
    from isoquant_lib.common import junctions_from_blocks

    def sel7():
        if CTX.phase != "collect" or CTX.cur_aln is None or not CTX.S.wants("U7") or CTX.guard:
            return False
        return keep(CTX.bkey() + ":a%d" % CTX.cur_aln[1], CTX.S.rate("U7"))

    def f_corr(orig):
        def w(self, alignment_info, read_assignment):
            if not sel7():
                return orig(self, alignment_info, read_assignment)
            a = alignment_info.alignment
            cp = alignment_info.combined_profile
            inp = {"fn": "correct_assigned_read", "read": CTX.read_ident(), "index": CTX.index(self.gene_info),
                   "params": enc.assign_params(self.params), "aln": enc.aln(a),
                   "read_exons": enc.ivs(alignment_info.read_exons),
                   "read_start": int(alignment_info.read_start), "read_end": int(alignment_info.read_end),
                   "intron_profile": enc.prof(cp.read_intron_profile) if cp is not None else None,
                   "assignment": enc.assignment(read_assignment),
                   "ref": _exon_windows(self.chr_record, alignment_info.read_exons) if self.chr_record is not None else []}
            res = orig(self, alignment_info, read_assignment)
            _write("U7", CTX.local("U7", "correct_assigned_read", "a%d" % CTX.cur_aln[1]), inp,
                   {"corrected_exons": enc.ivs(res)})
            return res
        return w
    patch(ExonCorrector, "correct_assigned_read", f_corr)

    def f_strand(orig):
        def w(self, read_assignment):
            if not sel7():
                return orig(self, read_assignment)
            cr = self.strand_detector.chr_record
            introns = read_assignment.corrected_introns or []
            sites = []
            for s, e in introns:
                sites.append([[int(s), int(e)], str(cr[s - 1:s + 1]) if cr is not None else "",
                              str(cr[e - 2:e]) if cr is not None else ""])
            inp = {"fn": "get_assignment_strand", "read": CTX.read_ident(), "assignment": enc.assignment(read_assignment),
                   "polya": enc.polya(read_assignment.polya_info), "exons": enc.ivs(read_assignment.exons),
                   "corrected_introns": enc.ivs(introns), "intron_sites": sites}
            res = orig(self, read_assignment)
            _write("U7", CTX.local("U7", "get_assignment_strand", "a%d" % CTX.cur_aln[1]), inp, {"strand": res})
            return res
        return w
    patch(AlignmentCollector, "get_assignment_strand", f_strand)

    orig_anti = AlignmentCollector.__dict__["check_antisense"].__func__

    def w_anti(read_assignment):
        if not sel7():
            return orig_anti(read_assignment)
        inp = {"fn": "check_antisense", "read": CTX.read_ident(), "strand": read_assignment.strand,
               "assignment": enc.assignment(read_assignment)}
        res = orig_anti(read_assignment)
        _write("U7", CTX.local("U7", "check_antisense", "a%d" % CTX.cur_aln[1]), inp,
               {"assignment": enc.assignment(read_assignment)})
        return res
    AlignmentCollector.check_antisense = staticmethod(w_anti)


# ============================================================================ U8
def _flush_u8():
    p = CTX.u8_pending
    CTX.u8_pending = None
    if p is None or p.get("skip"):
        return
    CTX.S.write("U8", p["key"], p["in"], {"reads": p["reads"]})


def _install_u8():
    from isoquant_lib.assignment.assignment_io import TmpFileAssignmentPrinter

    def f_gi(orig):
        def w(self, gene_info):
            _flush_u8()
            res = orig(self, gene_info)
            if CTX.S.wants("U8") and CTX.phase == "collect" and keep(CTX.bkey(), CTX.S.block_rate):
                info = CTX.block_info
                CTX.u8_pending = {
                    "key": CTX.bkey() + ":collect_block:n.0",
                    "in": {"fn": "collect_block", "chr": CTX.chr, "block_no": CTX.block, "region": info.get("region"),
                           "index": CTX.index(gene_info), "genic": bool(info.get("genic", not gene_info.empty())),
                           "skip_read_fraction": int(info.get("skip", 1)), "max_coverage": info.get("max_coverage")},
                    "reads": []}
            return res
        return w
    patch(TmpFileAssignmentPrinter, "add_gene_info", f_gi)

    def f_ri(orig):
        def w(self, read_assignment):
            p = CTX.u8_pending
            if p is not None and not p.get("skip") and read_assignment.assignment_type is not None \
                    and read_assignment.isoform_matches is not None:
                if len(p["reads"]) >= CTX.S.block_cap:
                    p["skip"] = True
                    p["reads"] = []
                    CTX.skipped_blocks.append({"chr": CTX.chr, "block_no": CTX.block, "reason": "U8 block_cap"})
                else:
                    p["reads"].append(enc.read_record(read_assignment, len(p["reads"])))
            return orig(self, read_assignment)
        return w
    patch(TmpFileAssignmentPrinter, "add_read_info", f_ri)

    def f_close(orig):
        def w(self):
            if not self.dumper.closed:
                _flush_u8()
            return orig(self)
        return w
    patch(TmpFileAssignmentPrinter, "close", f_close)


# ============================================================================ U9
def _install_u9():
    from isoquant_lib.assignment.multimap_resolver import MultimapResolver

    def f(orig):
        def w(self, assignment_list):
            if not CTX.S.wants("U9") or not assignment_list:
                return orig(self, assignment_list)
            a0 = assignment_list[0]
            k = "%s:%s:%s:-:" % (CTX.sample, CTX.slice, a0.chr_id)
            ckey = ("U9", k, a0.read_id)
            n = CTX.counters[ckey]
            CTX.counters[ckey] += 1
            key = k + "resolve:%s.%d" % (a0.read_id, n)
            if not keep(key, CTX.S.rate("U9")):
                return orig(self, assignment_list)
            inp = {"fn": "resolve", "strategy": int(self.strategy.value), "assignments": [enc.basic(a) for a in assignment_list]}
            res = orig(self, assignment_list)
            CTX.S.write("U9", key, inp, {"assignments": [enc.basic(a) for a in (res or [])]})
            return res
        return w
    patch(MultimapResolver, "resolve", f)


# ============================================================================ U10
def _model_name(model):
    return getattr(model, "_iq_name", None)


def _install_u10():
    from xgboost import XGBClassifier
    from isoquant_lib.terminal_prediction import terminal_peaks as TP
    from isoquant_lib.terminal_prediction.terminal_counter import TerminalCounter
    import numpy as np
    import pandas as pd

    def f_load(orig):
        def w(self, fname, *a, **kw):
            r = orig(self, fname, *a, **kw)
            b = os.path.basename(str(fname)).lower()
            self._iq_name = "polya" if "polya" in b else ("tss" if "tss" in b else b)
            return r
        return w
    patch(XGBClassifier, "load_model", f_load)

    def f_pred(orig):
        def w(self, X, *a, **kw):
            if CTX.guard or a or kw or not CTX.S.wants("U10a"):
                return orig(self, X, *a, **kw)
            caller = sys._getframe(1).f_code.co_name
            skey = CTX.bkey() + ":xgb:" + caller
            if not CTX.selected("U10a", skey):
                return orig(self, X, *a, **kw)
            arr = np.asarray(X, dtype=float)
            res = orig(self, X)
            CTX.guard += 1
            try:
                proba = self.predict_proba(X)[:, 1]
                margin = orig(self, X, output_margin=True)
            finally:
                CTX.guard -= 1
            inp = {"fn": "xgb_predict", "model": _model_name(self), "caller": caller,
                   "features": [[fval(v) for v in row] for row in arr.tolist()]}
            out = {"labels": [int(v) for v in np.asarray(res).tolist()], "proba": [fval(v) for v in proba.tolist()],
                   "margin": [fval(v) for v in np.asarray(margin).ravel().tolist()]}
            _write("U10a", CTX.local("U10a", "xgb_predict", "n"), inp, out)
            return res
        return w
    patch(XGBClassifier, "predict", f_pred)

    orig_dpb = TP.detect_peaks_batch

    def w_dpb(histograms, model):
        if CTX.guard or not CTX.S.wants("U10a") or not CTX.selected("U10a", CTX.bkey() + ":dpb"):
            return orig_dpb(histograms, model)
        hs = [None if not h else [[int(p), int(c)] for p, c in h.items()] for h in histograms]
        CTX.guard += 1
        try:
            records, rows = [], []
            for i, counts in enumerate(histograms):
                if not counts:
                    records.append(None)
                    continue
                r = TP._histogram_record(dict(counts))
                rr = {"start": int(r["start"]), "mode": int(r["mode"]), "var": fval(r["var"]), "skew": fval(r["skew"]),
                      "entropy": fval(r["entropy"]), "mean_height": fval(r["mean_height"]),
                      "peak_count": int(r["peak_count"]), "histogram": [int(x) for x in r["histogram"]],
                      "raw_peaks": None}
                if r["peak_count"]:
                    rp = r["raw_peaks"]
                    rr["raw_peaks"] = {"peak_location": [int(x) for x in rp["peak_location"]],
                                       "peak_width": [fval(x) for x in rp["peak_width"]],
                                       "peak_left": [int(x) for x in rp["peak_left"]],
                                       "peak_right": [int(x) for x in rp["peak_right"]],
                                       "peak_heights": [int(x) for x in rp["peak_heights"]]}
                    for pk in TP._rank_peaks(r):
                        rows.append((i, [r["var"], r["skew"], r["peak_count"], pk["peak_width"], r["entropy"],
                                         r["mean_height"], pk["peak_heights"], pk["relative_height"]]))
                records.append(rr)
            row_out = []
            if rows:
                feats = pd.DataFrame([x[1] for x in rows], columns=TP.FEATURE_COLUMNS).astype(float)
                lab = model.predict(feats)
                proba = model.predict_proba(feats)[:, 1]
                margin = model.predict(feats, output_margin=True)
                for j, (i, fr) in enumerate(rows):
                    row_out.append({"hist": i, "features": [fval(v) for v in feats.iloc[j].tolist()],
                                    "proba": fval(proba[j]), "margin": fval(np.asarray(margin).ravel()[j]),
                                    "label": int(lab[j])})
        finally:
            CTX.guard -= 1
        res = orig_dpb(histograms, model)
        out = {"records": records, "rows": row_out,
               "peaks": [[[int(p.position), int(p.count), int(p.left), int(p.right)] for p in pl] for pl in res]}
        _write("U10a", CTX.local("U10a", "detect_peaks_batch", "n"),
               {"fn": "detect_peaks_batch", "model": _model_name(model), "histograms": hs}, out)
        return res
    TP.detect_peaks_batch = w_dpb

    def f_flush(orig):
        def w(self):
            if (not self.ignore_read_groups or self._collecting_training or not CTX.selected("U10b", CTX.bkey())
                    or CTX.phase != "construct"):
                return orig(self)
            kind = "polya" if "polya" in type(self).__name__.lower() else "tss"
            tr = [[t, {"chr": e.get("chr"), "gene_id": e.get("gene_id"), "data": [int(x) for x in e["data"]],
                       "annotated": int(e["annotated"]) if e.get("annotated") is not None else None}]
                  for t, e in self.transcripts.items()]
            inp = {"fn": "flush", "kind": kind, "block": CTX.block_ref(), "transcripts": tr}
            res = orig(self)
            _write("U10b", CTX.local("U10b", "flush", kind), inp,
                   {"predictions": [[t, int(p)] for t, p in self.last_gene_predictions]})
            return res
        return w
    patch(TerminalCounter, "flush", f_flush)


# ============================================================================ U11 / U12
def graph_dump(g):
    ic = g.intron_collector
    return {
        "clustered_introns": [[enc.iv(k), int(v)] for k, v in ic.clustered_introns.items()],
        "intron_correction_map": [[enc.iv(k), enc.iv(v)] for k, v in ic.intron_correction_map.items()],
        "discarded_introns": [enc.iv(x) for x in sorted(ic.discarded_introns)],
        "outgoing": [[enc.vertex(k), [enc.vertex(x) for x in sorted(v)]] for k, v in g.outgoing_edges.items()],
        "incoming": [[enc.vertex(k), [enc.vertex(x) for x in sorted(v)]] for k, v in g.incoming_edges.items()],
        "starting_known_positions": [[enc.iv(k), [int(x) for x in v]] for k, v in g.starting_known_positions.items()],
        "terminal_known_positions": [[enc.iv(k), [int(x) for x in v]] for k, v in g.terminal_known_positions.items()],
        "max_coverage": int(g.max_coverage),
    }


def _preds(p):
    return None if p is None else [[t, int(x)] for t, x in p]


def _install_u11_u12():
    from isoquant_lib.model_construction.intron_graph import IntronGraph
    from isoquant_lib.model_construction.intron_path import IntronPathStorage
    from isoquant_lib.model_construction.context import ModelConstructionContext

    def f_graph(orig):
        def w(self, params, gene_info, read_assignments, polya_predictions=None, tss_predictions=None):
            sel = CTX.selected("U11", CTX.bkey())
            if sel:
                inp = {"fn": "IntronGraph.__init__", "block": CTX.block_ref(), "index": CTX.index(gene_info),
                       "polya_predictions": _preds(polya_predictions), "tss_predictions": _preds(tss_predictions)}
            orig(self, params, gene_info, read_assignments, polya_predictions=polya_predictions,
                 tss_predictions=tss_predictions)
            if sel:
                loc = CTX.local("U11", "IntronGraph.__init__", "n")
                CTX.graph_key = CTX.bkey() + ":" + loc
                _write("U11", loc, inp, {"graph": graph_dump(self)})
        return w
    patch(IntronGraph, "__init__", f_graph)

    def f_fill(orig):
        def w(self, read_assignments):
            res = orig(self, read_assignments)
            if CTX.selected("U12", CTX.bkey()):
                loc = CTX.local("U12", "IntronPathStorage.fill", "n")
                CTX.paths_key = CTX.bkey() + ":" + loc
                _write("U12", loc, {"fn": "IntronPathStorage.fill", "block": CTX.block_ref(), "graph": CTX.graph_key},
                       {"paths": [[enc.path(p), int(c)] for p, c in self.paths.items()],
                        "fl_paths": [enc.path(p) for p in sorted(self.fl_paths)],
                        "paths_to_reads": [[enc.path(p), [enc.ridx(a) for a in v]] for p, v in self.paths_to_reads.items()]})
            return res
        return w
    patch(IntronPathStorage, "fill", f_fill)

    def f_known(orig):
        def w(self, gene_info, s="known"):
            res = orig(self, gene_info, s)
            if CTX.selected("U12", CTX.bkey()):
                _write("U12", CTX.local("U12", "get_known_spliced_isoforms", "n"),
                       {"fn": "get_known_spliced_isoforms", "block": CTX.block_ref(), "graph": CTX.graph_key,
                        "index": CTX.index(gene_info)},
                       {"known": [[enc.path(p), t] for p, t in res.items()]})
            return res
        return w
    patch(ModelConstructionContext, "get_known_spliced_isoforms", f_known)


# ============================================================================ U13 - U19
def _base_in(fn, obj):
    return {"fn": fn, "block": CTX.block_ref(), "index": CTX.index(obj.ctx.gene_info),
            "graph": CTX.graph_key, "paths": CTX.paths_key}


def _store_method(cls, name, unit, with_reads_arg=False, extra=None):
    def f(orig):
        if with_reads_arg:
            def w(self, read_assignments):
                if not CTX.selected(unit, CTX.bkey()):
                    return orig(self, read_assignments)
                inp = _base_in(name, self)
                inp["store"] = enc.store(self.store, self.ctx)
                inp["reads"] = [enc.ridx(a) for a in read_assignments]
                if extra is not None:
                    inp.update(extra(self))
                res = orig(self, read_assignments)
                _write(unit, CTX.local(unit, name, "n"), inp, {"store": enc.store(self.store, self.ctx)})
                return res
        else:
            def w(self):
                if not CTX.selected(unit, CTX.bkey()):
                    return orig(self)
                inp = _base_in(name, self)
                inp["store"] = enc.store(self.store, self.ctx)
                if extra is not None:
                    inp.update(extra(self))
                res = orig(self)
                _write(unit, CTX.local(unit, name, "n"), inp, {"store": enc.store(self.store, self.ctx)})
                return res
        return w
    patch(cls, name, f)


def _u13_intron_sites(self):
    """U13 `intron_sites`: for every vertex of every FL path's intron_path (path[1:-1]), the bases
    get_intron_strand reads: (vertex, chr[a-1:a+1], chr[b-2:b]) verbatim; sorted, unique."""
    cr = self.strand_detector.chr_record
    vs = set()
    for path in self.ctx.path_storage.fl_paths:
        for v in path[1:-1]:
            vs.add((int(v[0]), int(v[1])))
    sites = []
    for a, b in sorted(vs):
        sites.append([[a, b], str(cr[a - 1:a + 1]) if cr is not None else "",
                      str(cr[b - 2:b]) if cr is not None else ""])
    return {"intron_sites": sites}


U16_PAD_L, U16_PAD_R = 65, 64


def _u16_store_ref(self):
    """U16 store-level `ref`: windows [b-65, b+64) around every exon boundary b of every non-known
    model with >= 2 exons in the store before the call, merged when overlapping (sorted by start)."""
    from isoquant_lib.gene_info import TranscriptModelType
    if self.chr_record is None:
        return {"ref": []}
    iv = []
    for m in self.store.transcript_model_storage:
        if m.transcript_type == TranscriptModelType.known or len(m.exon_blocks) < 2:
            continue
        for e in m.exon_blocks:
            for b in e:
                iv.append([b - U16_PAD_L, b + U16_PAD_R])
    iv.sort()
    wins = []
    for a, b in iv:
        if wins and a <= wins[-1][1]:
            wins[-1][1] = max(wins[-1][1], b)
        else:
            wins.append([a, b])
    return {"ref": [_refwin(self.chr_record, a, b) for a, b in wins]}


def _install_construct():
    from isoquant_lib.model_construction.fl_graph_constructor import GraphBasedConstructor
    from isoquant_lib.model_construction.assignment_based_constructor import AssignmentBasedConstructor
    from isoquant_lib.model_construction.model_filter import ModelFilter
    from isoquant_lib.model_construction.end_processor import TranscriptEndProcessor
    from isoquant_lib.model_construction.read_assigner import ConstructionReadAssigner
    from isoquant_lib.model_construction.transcript_to_gene_joiner import TranscriptToGeneJoiner
    from isoquant_lib.model_construction.model_read_counter import ModelReadCounter

    _store_method(GraphBasedConstructor, "construct_fl_isoforms", "U13", extra=_u13_intron_sites)
    _store_method(AssignmentBasedConstructor, "construct_assignment_based_isoforms", "U14", True)
    _store_method(ModelFilter, "pre_filter_transcripts", "U15a")
    _store_method(ModelFilter, "filter_transcripts", "U15a")
    _store_method(TranscriptEndProcessor, "_add_known_alternative_end_models", "U15b")
    _store_method(TranscriptEndProcessor, "_drop_duplicate_alt_end_models", "U15b")
    _store_method(ModelFilter, "correct_transcript_splice_sites", "U16", extra=_u16_store_ref)
    _store_method(ConstructionReadAssigner, "assign_reads_to_models", "U17", True)
    _store_method(ConstructionReadAssigner, "resolve_read_ownership_conflicts", "U17")

    def f_sim(orig):
        def w(self, model_storage, remove_unconfirmed=False):
            if not CTX.selected("U15a", CTX.bkey()):
                return orig(self, model_storage, remove_unconfirmed)
            inp = _base_in("detect_similar_isoforms", self)
            inp["models"] = [enc.model(m) for m in model_storage]
            inp["remove_unconfirmed"] = bool(remove_unconfirmed)
            res = orig(self, model_storage, remove_unconfirmed)
            _write("U15a", CTX.local("U15a", "detect_similar_isoforms", "n"), inp,
                   {"to_substitute": [[k, v] for k, v in res.items()]})
            return res
        return w
    patch(ModelFilter, "detect_similar_isoforms", f_sim)

    def f_ends(orig):
        def w(self, transcript_model, assigned_reads):
            if not CTX.selected("U15b", CTX.bkey()):
                return orig(self, transcript_model, assigned_reads)
            inp = _base_in("correct_novel_transcript_ends", self)
            inp["model"] = enc.model(transcript_model)
            inp["reads"] = [enc.ridx(a) for a in assigned_reads]
            res = orig(self, transcript_model, assigned_reads)
            _write("U15b", CTX.local("U15b", "correct_novel_transcript_ends", transcript_model.transcript_id), inp,
                   {"model": enc.model(transcript_model)})
            return res
        return w
    patch(TranscriptEndProcessor, "correct_novel_transcript_ends", f_ends)

    def f_mark(orig):
        def w(self, model):
            if not CTX.selected("U15b", CTX.bkey()):
                return orig(self, model)
            inp = _base_in("_mark_terminal_confirmation", self)
            inp["model"] = enc.model(model)
            inp["reads"] = [enc.ridx(a) for a in self.store.transcript_read_ids.get(model.transcript_id, [])]
            res = orig(self, model)
            _write("U15b", CTX.local("U15b", "_mark_terminal_confirmation", model.transcript_id), inp,
                   {"model": enc.model(model)})
            return res
        return w
    patch(TranscriptEndProcessor, "_mark_terminal_confirmation", f_mark)

    def f_derive(orig):
        def w(self, source_model, assigned_reads):
            if not CTX.selected("U15b", CTX.bkey()):
                return orig(self, source_model, assigned_reads)
            inp = _base_in("derive_alternative_end_models", self)
            inp["model"] = enc.model(source_model)
            inp["reads"] = [enc.ridx(a) for a in assigned_reads]
            inp["id_counter"] = int(self.ctx.id_distributor.value)
            res = orig(self, source_model, assigned_reads)
            _write("U15b", CTX.local("U15b", "derive_alternative_end_models", source_model.transcript_id), inp,
                   {"result": [[enc.model(m), [enc.ridx(a) for a in rs]] for m, rs in res],
                    "id_counter": int(self.ctx.id_distributor.value)})
            return res
        return w
    patch(TranscriptEndProcessor, "derive_alternative_end_models", f_derive)

    def f_ss(orig):
        def w(self, model, reads):
            if not CTX.selected("U16", CTX.bkey()):
                return orig(self, model, reads)
            inp = _base_in("_correct_one_model_splice_sites", self)
            inp["model"] = enc.model(model)
            inp["reads"] = [{"idx": enc.ridx(r), "start": int(r.exons[0][0]) if r.exons else None,
                             "end": int(r.exons[-1][1]) if r.exons else None,
                             "boundary_deletions": enc.ivs(r.boundary_deletions)} for r in reads]
            wins = []
            if self.chr_record is not None:
                for e in model.exon_blocks:
                    for b in e:
                        wins.append(_refwin(self.chr_record, b - 65, b + 64))
            inp["ref"] = wins
            res = orig(self, model, reads)
            _write("U16", CTX.local("U16", "_correct_one_model_splice_sites", model.transcript_id), inp,
                   {"exons": enc.ivs(res) if res is not None else None})
            return res
        return w
    patch(ModelFilter, "_correct_one_model_splice_sites", f_ss)

    def f_join(orig):
        def w(self):
            if not CTX.selected("U18", CTX.bkey()):
                return orig(self)
            inp = {"fn": "join_transcripts", "block": CTX.block_ref(), "index": CTX.index(self.gene_info),
                   "models": [enc.model(m) for m in self.transcipt_model_storage]}
            res = orig(self)
            _write("U18", CTX.local("U18", "join_transcripts", "n"), inp, {"models": [enc.model(m) for m in res]})
            return res
        return w
    patch(TranscriptToGeneJoiner, "join_transcripts", f_join)

    def f_final(orig):
        def w(self, read_assignments):
            if not CTX.selected("U19", CTX.bkey()):
                return orig(self, read_assignments)
            inp = _base_in("build_model_read_assignments", self)
            inp["store"] = enc.store(self.store, self.ctx)
            inp["reads"] = [enc.ridx(a) for a in read_assignments]
            res = orig(self, read_assignments)
            out = [{"read_idx": enc.ridx(src), "assignment": enc.assignment(ra)}
                   for src, ra in zip(read_assignments, self.model_read_assignments)]
            _write("U19", CTX.local("U19", "build_model_read_assignments", "n"), inp, {"assignments": out})
            return res
        return w
    patch(ModelReadCounter, "build_model_read_assignments", f_final)


# ============================================================================ U20
class _Tee:
    def __init__(self, fh):
        self.fh = fh
        self.buf = []

    def write(self, s):
        self.buf.append(s)
        return self.fh.write(s)

    def __getattr__(self, k):
        return getattr(self.fh, k)


def _counter_desc(c):
    return {"name": os.path.basename(c.output_counts_file_name), "grouped": not c.ignore_read_groups,
            "group_index": int(getattr(c, "group_index", 0))}


def _is_model_counter(c):
    n = os.path.basename(getattr(c, "output_counts_file_name", "") or "")
    return ".discovered_" in n


def _install_u20():
    from isoquant_lib.model_construction.transcript_printer import GFFPrinter
    from isoquant_lib.quantification.long_read_counter import AssignedFeatureCounter

    def f_dump(orig):
        def w(self, gene_info, transcript_model_storage):
            kind = "transcript_models" if self.model_fname.endswith(".transcript_models.gtf") else "extended_annotation"
            if (kind != "transcript_models" or CTX.phase != "construct" or not CTX.S.wants("U20")
                    or not CTX.block_selected):
                return orig(self, gene_info, transcript_model_storage)
            inp = {"fn": "gff_dump", "block": CTX.block_ref(), "printer": kind, "index": CTX.index(gene_info),
                   "models": [enc.model(m) for m in (transcript_model_storage or [])],
                   "printed_gene_ids": sorted(self.printed_gene_ids),
                   "exon_id_counter": int(self.exon_id_storage.id_distributor.value)}
            real = self.out_gff
            tee = _Tee(real)
            self.out_gff = tee
            try:
                res = orig(self, gene_info, transcript_model_storage)
            finally:
                self.out_gff = real
            _write("U20", CTX.local("U20", "gff_dump", "n"), inp,
                   {"lines": tee.buf, "exon_id_counter": int(self.exon_id_storage.id_distributor.value)})
            return res
        return w
    patch(GFFPrinter, "dump", f_dump)

    def gname(c, g):
        try:
            return str(c._get_group_name(g))
        except Exception:
            return str(g)

    def f_add(orig):
        def w(self, read_assignment=None):
            if (read_assignment is None or not _is_model_counter(self) or CTX.phase != "construct"
                    or not CTX.S.wants("U20") or not CTX.block_selected):
                return orig(self, read_assignment)
            src = _find_local_read(sys._getframe(1))
            base = "r%d" % src._iq_idx if src is not None else "n"
            skey = CTX.bkey() + ":" + base + ":" + os.path.basename(self.output_counts_file_name)
            if not keep(skey, CTX.S.rate("U20")):
                return orig(self, read_assignment)
            try:
                feats = list(self.assignment_extractor.get_features(read_assignment) or [])
            except Exception:
                feats = []
            before = {f: dict(self.feature_counter[f].data) if f in self.feature_counter else {} for f in feats}
            amb0, na0 = self.ambiguous_reads, self.not_assigned_reads
            conf0 = set(self.confirmed_features)
            inp = {"fn": "counter_add_read_info", "counter": _counter_desc(self),
                   "read": {"read_idx": src._iq_idx} if src is not None else None,
                   "assignment": enc.assignment(read_assignment),
                   "read_group": [str(g) for g in (read_assignment.read_group or [])]}
            res = orig(self, read_assignment)
            fo = []
            for f in feats:
                after = dict(self.feature_counter[f].data) if f in self.feature_counter else {}
                groups = list(after.keys()) + [g for g in before[f] if g not in after]
                fo.append([f, [[gname(self, g), fval(before[f].get(g, 0.0)), fval(after.get(g, 0.0))] for g in groups]])
            out = {"features": fo, "ambiguous": [int(amb0), int(self.ambiguous_reads)],
                   "not_assigned": [int(na0), int(self.not_assigned_reads)],
                   "confirmed_added": sorted(set(self.confirmed_features) - conf0)}
            _write("U20", CTX.local("U20", "counter_add_read_info", base), inp, out)
            return res
        return w
    patch(AssignedFeatureCounter, "add_read_info", f_add)

    def f_cdump(orig):
        def w(self):
            if (not _is_model_counter(self) or not CTX.S.wants("U20") or CTX.phase != "construct"
                    or (not self.feature_counter and not self.all_features and not self.not_assigned_reads
                        and not self.ambiguous_reads)):
                return orig(self)
            inp = {"fn": "counter_dump", "counter": _counter_desc(self),
                   "feature_counter": [[f, [[gname(self, g), fval(v)] for g, v in d.data.items()]]
                                       for f, d in self.feature_counter.items()],
                   "all_features": sorted(x for x in self.all_features if x is not None),
                   "confirmed_features": sorted(self.confirmed_features),
                   "ambiguous": int(self.ambiguous_reads), "not_assigned": int(self.not_assigned_reads),
                   "not_aligned": int(self.not_aligned_reads)}
            res = orig(self)
            files = []
            for p in (self.output_counts_file_name, getattr(self, "output_stats_file_name", None),
                      getattr(self, "usable_file_name", None)):
                if p and os.path.exists(p):
                    with open(p) as fh:
                        files.append([os.path.basename(p), fh.read()])
            CTX.S.write("U20", "%s:%s:%s:-:counter_dump:%s" % (CTX.sample, CTX.slice, CTX.chr, inp["counter"]["name"]),
                        inp, {"files": files})
            return res
        return w
    patch(AssignedFeatureCounter, "dump", f_cdump)

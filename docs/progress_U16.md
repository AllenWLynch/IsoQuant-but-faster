# U16 progress (splice-site correction) — DONE

- `crates/iq_model/src/splice_site.rs`: `correct_exon_sites(exons, strand, reads, &RefSlice)` (core),
  `correct_model_sites(&TranscriptModel, &[SiteRead], &RefSlice)`, `valid_exon_chain`,
  `correct_transcript_splice_sites(&mut ModelStore, &BlockStore, &RefSlice)` (caller checks
  `params.splice_site_correction`). Registered in iq_model lib.rs.
- Fixtures: `_correct_one_model_splice_sites` 172/172, `correct_transcript_splice_sites` 67/67 exact, both builds,
  both modes. Only 1 fixture case actually changes exons, so a randomized differential test against the real 4.1
  Python functions (scratchpad gen2.py → U16-format jsonl, 40k cases, ~280 corrections, lowercase ref, 1-bp exons,
  '.' strand) was run: 40000/40000 exact.
- Order: mode must be unique (ties → -1, never accepted); per-site dict is keyed by position and corrected
  locations become a set → no dict/set/Counter order reaches output. Nothing to report to main.
- Runner notes: per-boundary `ref` windows stitched into one RefSlice (gaps = 'N'). Store-level runner reads
  `data/ref/genome.fa` (copy of U13's private `read_region`/`genome_path`; could be moved to fixture_util).
- Shared target/ dir races with other agents' builds (link errors); used CARGO_TARGET_DIR in scratchpad.

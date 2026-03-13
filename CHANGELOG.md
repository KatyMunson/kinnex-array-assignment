# Changelog

## [Unreleased] — 2026-03-11

### Fixed
- Corrected 10 bugs found during code review, ranging from crash-on-valid-input
  (ZeroDivisionError in split_skera_by_library.py, division-by-zero in
  visualize_posteriors.py) to silently wrong output (or-None swallowing valid
  zero counts in aggregate_pipeline_qc.py, stale baseline weights in training
  scripts, split() vs split("\t") inconsistency in assign_kinnex.py).
- Training scripts (optimize_barcode_weights.py, optimize_thresholds.py,
  optimize_thresholds_v2.py) now seed their baseline/comparison values from
  the assignment file provenance header instead of hardcoded defaults, so
  output correctly reflects the weights and thresholds actually used.

### Changed
- Lima pass 1 thresholds (min_score, min_score_lead, min_ref_span) moved from
  hardcoded Snakefile values to config.yaml (lima_pass1 section), replacing
  a dead lima: peek_guess: true entry that was never read.

## [Unreleased] — 2026-03-13

### Fixed
- `test_script/` Snakefile rules (`Snakefile_subsample` and `Snakefile_qc`)
  now capture per-rule logs: `log:` directives and Python `logging.FileHandler`
  wrappers added to all `run:` blocks so errors are written to per-rule log
  files rather than lost in Snakemake's stderr stream.
- Fixed crash in `Snakefile_subsample` `bam_path()`: `.loc[]` on a non-unique
  manifest index returned a Series, causing `not path` to raise
  `ValueError: The truth value of a Series is ambiguous`. Added uniqueness
  validation at manifest load time and switched to `.at[]` for scalar access.
- Fixed silent wrong output in `Snakefile_subsample` proportional subsampling:
  `sc_value` from the pool TSV was read as a string; `int * str` produced
  string repetition instead of multiplication. Fixed by casting to `float()`.
- Fixed all-zero aggregate QC tables in `Snakefile_qc` `aggregate_qc` rule:
  pool/knbc were extracted from the result directory name by splitting on
  `_KN_`, but pool names ending with `_KN` caused the split to hit the wrong
  delimiter, producing values that matched nothing in `unique_pools`/
  `unique_knbcs`. Fixed by iterating `zip(POOLS, KNBCS, input.summaries)`
  instead of parsing the path.

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

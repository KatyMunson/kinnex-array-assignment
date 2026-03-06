# Refactor Notes

Collected during testing of the initial refactored pipeline (Feb 2026).
Items to address before or during the next major revision.

---

## Rule granularity — Snakefile_subsample

**Issue:** `merge_and_finalize` is a monolithic rule that does three distinct things:
1. Merges subsampled BAMs
2. Creates the read lookup file (`*_reads.txt`)
3. Creates the arrays.txt file

These are bundled together, so regenerating any one output requires rerunning
the entire rule — including the expensive BAM merge. For example, fixing the
whitespace bug in arrays.txt generation required a full rerun even though the
merged BAMs and lookup files were unchanged and correct.

**Desired outcome:** Each output should be in its own rule (or at minimum,
arrays.txt should be a separate rule with no BAM dependency). This would also
open the door to a standalone `generate_arrays.py` script that could be run
directly from pool.tsv without needing any BAM data at all — since arrays.txt
only depends on the bc_set column, not on the reads.

**Related:** Do a broader audit of both Snakefile_subsample and Snakefile_qc
for other places where unrelated outputs are bundled into a single rule.
Principle: if two outputs don't share a slow/expensive dependency, they should
be in separate rules.

---

## Barcode terminology rename

**Issue:** Column headers and constants still use "informative" / "uninformative"
terminology, which is a misnomer — shared barcodes still contribute signal.

**Proposed rename:**
- Informative → Specific (library-unique barcodes)
- Uninformative → Shared (barcodes appearing in ≤N libraries, still weighted)
- Extraneous → Discordant (barcodes from a different library in the pool)
- `INF_WEIGHT` → `SPECIFIC_WEIGHT`
- `MAX_UNINF_WEIGHT` → `MAX_SHARED_WEIGHT`

**Affected files:** assign_kinnex.py (constants, output column headers),
Snakefile_qc (column references), optimize_barcode_weights.py,
optimize_thresholds.py, optimize_thresholds_v2.py, README.md.

---

## Utility module (utils.py)

**Issue:** Several scripts duplicate the same logic: loading pool.tsv/manifest,
the bc_set → IsoSeqX barcode name translation, ZMW name normalisation
(stripping /ccs suffix), and data loading patterns.

**Desired outcome:** Extract shared logic into `scripts/utils.py` and import
from there. Reduces maintenance surface and ensures consistency (e.g. the
/ccs stripping fix would only need to be made in one place).

---

## MIN_INF gate (minimum informative barcodes)

**Issue:** The classifier gates on MIN_OBS (total segments) but not on a
minimum count of specific/informative barcodes. A ZMW with many shared-only
observations could pass MIN_OBS but have low discriminating power.

**Desired outcome:** Add a MIN_INF_HIGH_CONF / MIN_INF_LOW_CONF parameter
to assign_kinnex.py, and use optimize_thresholds_v2.py to find the optimal
value empirically. See earlier design discussion for context.

---

---

## Code commenting pass

**Issue:** Comments throughout the codebase are currently written at a
developer level — they describe *what* the code does but often not *why*,
and assume familiarity with Python idioms and bioinformatics conventions.

**Desired outcome:** A thorough commenting pass aimed at a reader who is
comfortable with the biology but relatively new to coding. Specifically:

- Every function should have a docstring explaining its purpose, its inputs,
  and what it returns, in plain English
- Non-obvious logic (Bayesian scoring, Pareto frontier computation, ZMW name
  normalisation, wildcard constraint syntax, etc.) should have a plain-English
  explanation of what it's doing and why, not just a restatement of the code
- Section headers in long files (especially assign_kinnex.py and the
  Snakefiles) should give enough context that someone can navigate to the
  relevant part without reading everything
- Any place where a magic number or threshold appears, there should be a
  comment explaining where it came from and what changing it would do

**Priority files:** assign_kinnex.py (core algorithm), Snakefile (production
workflow), Snakefile_subsample (most complex logic), optimize_thresholds_v2.py
(Pareto/grid search is hard to follow without explanation).

---

## Conda environment directives

**Issue:** The Snakefiles don't currently use `conda:` directives per rule, so
the right tools need to be on PATH manually. The envs/ folder exists but isn't
wired into Snakemake's environment management.

**Desired outcome:** Add `conda:` directives to each rule that calls an
external tool, and run with `snakemake --use-conda`. Mapping would be:
- skera rules → envs/skera.yaml
- lima rules → envs/lima.yaml
- isoseq refine rules → envs/isoseq.yaml
- Python script rules → envs/python.yaml

This would make the pipeline fully self-contained and reproducible without
any manual environment setup.

---

## Document commit_configs.sh in README

**Issue:** `commit_configs.sh` exists at the repo root but isn't mentioned
anywhere in the README. A new user (or future-you) would have no idea it's
there or why they need it.

**Desired outcome:** Add a short section to README.md under a "Repository
Utilities" or "Git Workflow" heading explaining:
- Why config files are protected by default (they contain local paths)
- The three commands and when to use each
- The typical workflow for committing a template config change

---

## Pool design rules — theoretical justification

**Context:** Analysis of why TRAIN_3way pools produce zero HIGH_CONF reads led
to a theoretical derivation of maximum allowable barcode overlap given pool
size and scorer weights.

**Finding:** With 6 barcodes per library, the maximum shared barcodes that
still allow HIGH_CONF (posterior ≥ 0.840) depends on pool size:

| Pool size | Max shared BCs | Min unique BCs | Max achievable posterior |
|-----------|---------------|----------------|------------------------|
| 2 libs    | 4             | 2              | 0.900                  |
| 3 libs    | 3             | 3              | 0.931                  |
| 4 libs    | 3             | 3              | 0.900                  |
| 5 libs    | 3             | 3              | 0.871                  |
| 6 libs    | 3             | 3              | 0.844                  |

The production rule of max 2 shared barcodes is conservative but well-justified
— it provides comfortable headroom above the theoretical floor of 3 unique BCs,
and accounts for real-world noise that pushes observed posteriors below the
theoretical ceiling.

**Note:** This model assumes shared barcodes appear in exactly 2 libraries.
When a barcode is shared across 3+ libraries, the weight is
min(1/N_sharing, MAX_UNINF_WEIGHT) — capped at 0.2 regardless, so the ceiling
holds as long as MAX_UNINF_WEIGHT is the binding constraint.

**Desired outcome:** Add this analysis and the design rule derivation to the
README as a "Pool Design" section, providing scientific justification for
barcode overlap limits rather than just empirical rules.

---

## Standardize assignment file parsing via utils.py

**Context:** Adding provenance `#` comment headers to `assign_kinnex.py` output
required patching 8 separate `read_csv` / `open()` calls across 6 files to
handle comment lines. `split_skera_by_library.py` was missed on the first pass
and caused a runtime crash.

**Desired outcome:** Create `prod_script/scripts/utils.py` (or a shared module)
with a single `load_assignments(path)` function that:
- Skips `#` comment lines
- Skips the column header
- Returns a consistent data structure (pandas DataFrame or dict)
- Is imported by all scripts that read `.txt` assignment files:
  - `prod_script/scripts/split_skera_by_library.py`
  - `prod_script/scripts/aggregate_pipeline_qc.py`
  - `test_script/Snakefile_qc`
  - `test_script/scripts/visualize_posteriors.py`
  - `train_script/optimize_barcode_weights.py`
  - `train_script/optimize_thresholds.py`
  - `train_script/optimize_thresholds_v2.py`

Any future change to the assignment file format (new columns, different comment
syntax, etc.) would then require only one code change.

---

## Training scripts: read scoring parameters from assigned file headers

**Context:** `assign_kinnex.py` now writes scoring parameters as `#` comment
headers in every assigned output file (INF_WEIGHT, MAX_UNINF_WEIGHT,
EXTRANEOUS_PENALTY, POSTERIOR_HIGH_CONF, etc.). Currently the training scripts
(`optimize_barcode_weights.py`, `optimize_thresholds.py`,
`optimize_thresholds_v2.py`) ignore these headers entirely and assume the
parameters they were trained on match the current hardcoded constants in
`assign_kinnex.py`.

**Desired outcome:** When loading an assigned file, parse the `#` header lines
to extract the parameters used to generate that file. Use these to:
1. Warn if training files were generated with different parameters than each
   other (mixing parameter sets invalidates training)
2. Warn if training files were generated with different parameters than the
   current `assign_kinnex.py` constants (results may not transfer)
3. Optionally surface the git hash per file so provenance is fully tracked
   through the training run

A shared `parse_assignment_header(path)` helper (in the future `utils.py`)
would make this straightforward across all three scripts.

---

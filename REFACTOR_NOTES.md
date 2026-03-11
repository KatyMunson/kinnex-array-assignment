# Refactor Notes

Collected during development and testing of the pipeline.
Items to address before or during the next major revision.

---
##Sankey plot broken

**Issue:** `generate_sankey.py` / rule sankey_plot is not producing a
correct output. Needs investigation.
**To do:** Run generate_sankey.py manually against a completed results
directory and check stderr for errors. Likely a data format mismatch between
`qc_summary.tsv` (which feeds it) and what the script expects — possibly
predates the single/multi-library unification and assumes all samples are
multi-library.


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

## ~~Barcode terminology rename~~ — COMPLETED (Mar 2026)

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

**Resolution:** All terminology updated across the codebase. The old names
(`informative`, `uninformative`, `extraneous`) no longer appear in any file.

---

## ~~MIN_INF gate (minimum informative barcodes)~~ — COMPLETED (Mar 2026)

**Issue:** The classifier gates on MIN_OBS (total segments) but not on a
minimum count of specific/informative barcodes. A ZMW with many shared-only
observations could pass MIN_OBS but have low discriminating power.

**Desired outcome:** Add a MIN_INF_HIGH_CONF / MIN_INF_LOW_CONF parameter
to assign_kinnex.py, and use optimize_thresholds_v2.py to find the optimal
value empirically. See earlier design discussion for context.

**Resolution:** Added `MIN_SPECIFIC_HIGH_CONF` and `MIN_SPECIFIC_LOW_CONF`
constants to `assign_kinnex.py` (defaults: 1 and 0 respectively). Both are
applied in the classification logic and written to the assignment file header
for provenance. `optimize_thresholds_v2.py` was updated to grid-search over
these parameters.

---

## ~~Bug fix pass~~ — COMPLETED (Mar 2026)

**Issue:** Code review identified 10 bugs: 3 crash-on-valid-input, 5 silently
wrong output, 2 documentation mismatches. See CHANGELOG.md for the full list.

**Resolution:** All 10 fixed. Key fixes: ZeroDivisionError in
split_skera_by_library.py; or-None dropping valid zero FLNC counts in
aggregate_pipeline_qc.py; stale hardcoded baseline weights in training scripts
(now read from assignment file header); split() vs split("\t") inconsistency
between assign_kinnex.py and utils.py.

---

## ~~Lima pass 1 flags hardcoded in Snakefile~~ — COMPLETED (Mar 2026)

**Issue:** min_score, min_score_lead, min_ref_span were hardcoded in the
lima_isoseq rule shell block. Users troubleshooting low pass rates had to
edit the Snakefile directly. A dead lima: peek_guess: true key in config.yaml
gave the false impression that Lima options were already configurable.

**Resolution:** All three params moved to a lima_pass1: section in config.yaml
with notes on PacBio defaults and when to adjust them. The dead config key
was removed. Snakefile rule now reads from config[lima_pass1].

---

## Sentinel BAMs missing for zero-assignment libraries

**Issue:** split_skera_by_library.py creates sentinel (header-only) BAMs for
zero-read confidence tiers to prevent Snakemake MissingInputException. However,
it only creates them for libraries that appear in the assignment file. A library
present in the arrays file that received zero ZMW assignments will have no
sentinel BAM, which can cause a downstream MissingInputException.

**To do:** Pass the arrays file path as an additional input to the split script
so it can create sentinel BAMs for all expected libraries, not just those with
at least one assignment.

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

**Note:** The new Snakefile and utils.py written for the single/multi-library
unification (Mar 2026) follow the desired commenting standard and can serve
as a reference.

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

## ~~Extend utils.py to test_script and train_script~~ — COMPLETED (Mar 2026)

**Context:** utils.py was created in prod_script/scripts/ (Mar 2026) and
provides `load_assignments()`, `load_assignments_df()`, and
`parse_assignment_header()`. The following files still use their own
inline assignment file parsing and do not yet import from utils.py:

- `test_script/Snakefile_qc`
- `test_script/scripts/visualize_posteriors.py`
- `train_script/optimize_barcode_weights.py`
- `train_script/optimize_thresholds.py`
- `train_script/optimize_thresholds_v2.py`

**Desired outcome:** Make utils.py importable from these locations (e.g. by
moving it to a top-level `scripts/` directory or adding it to the Python
path in each script) and update all five files to use
`utils.load_assignments_df()` and `utils.parse_assignment_header()`.

**Note:** Training scripts should also use `parse_assignment_header()` to
warn when assignment files being used for training were generated with
different scoring parameters than each other, or than the current constants
in assign_kinnex.py.

**Resolution:** All five files now import utils.py via `sys.path.insert` to
`prod_script/scripts/` and use `utils.load_assignments_df()` and
`utils.parse_assignment_header()` for loading assignment data and headers.

**Follow-up fix (Mar 2026):** `visualize_posteriors.py` was additionally
corrected to use the header-parsed thresholds for all plot annotations.
Previously, threshold lines in the posterior distribution plots were hardcoded
to exploratory values (0.89, 0.93, 0.6) rather than the actual
`POSTERIOR_HIGH_CONF` / `POSTERIOR_LOW_CONF` / `MIN_OBS_HIGH_CONF` values
from the assignment file header. The plot now reads and displays the real
thresholds used during classification.

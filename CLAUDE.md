# CLAUDE.md — Kinnex Array Assignment Pipeline

## Project Overview

This is a **PacBio Kinnex Iso-Seq bioinformatics pipeline** for demultiplexing full-length RNA reads from multiplexed Kinnex arrays. The core problem it solves: when multiple libraries share a Kinnex barcode (KN_BC), each individual read (ZMW) must be assigned to its originating library using Bayesian classification before downstream processing.

The repository has three sub-pipelines:
- **`prod_script/`** — Production Snakemake workflow (Skera → Lima → Assign → Split → Lima → Refine → QC)
- **`test_script/`** — Synthetic pool generation and accuracy validation against ground truth
- **`train_script/`** — Parameter optimization scripts for thresholds and barcode weights

---

## Repository Structure

```
kinnex-array-assignment/
├── README.md                        # Comprehensive user-facing documentation
├── REFACTOR_NOTES.md                # Known issues and planned improvements
├── commit_configs.sh                # Git helper for config files with local paths
├── prod_script/
│   ├── Snakefile                    # Main workflow (719 lines)
│   ├── config.yaml                  # Cluster resources, script paths
│   ├── manifest.tab.example         # Example manifest input
│   ├── scripts/
│   │   ├── assign_kinnex.py         # Bayesian ZMW classifier (core logic)
│   │   ├── split_skera_by_library.py# Splits BAM by library + confidence tier
│   │   ├── aggregate_flnc_qc.py     # FLNC polyA count aggregation
│   │   ├── aggregate_pipeline_qc.py # Multi-stage QC report generator
│   │   ├── generate_sankey.py       # Interactive read-flow diagram (currently broken)
│   │   ├── plot_posteriors.py       # Posterior distribution PNG (no ground truth)
│   │   └── utils.py                 # Shared utilities (manifest, arrays, assignments)
│   └── envs/
│       ├── python.yaml              # pandas, numpy, pysam, sklearn, matplotlib, plotly
│       ├── skera.yaml               # PacBio Skera tool
│       ├── lima.yaml                # PacBio Lima tool
│       └── isoseq.yaml              # PacBio IsoSeq tool
├── test_script/
│   ├── Snakefile_subsample          # Creates synthetic multiplexed pools with ground truth
│   ├── Snakefile_qc                 # Validates assignment accuracy vs ground truth
│   ├── config_subsample.yaml        # Subsampling parameters
│   ├── config_qc.yaml               # QC validation paths
│   └── scripts/
│       └── visualize_posteriors.py  # Posterior plots with ground truth comparison
└── train_script/
    ├── optimize_barcode_weights.py  # ML-based weight optimization
    ├── optimize_thresholds.py       # 1D posterior threshold optimization
    └── optimize_thresholds_v2.py    # Multi-dimensional grid search (824 lines)
```

---

## Key Concepts

### Processing Modes

The pipeline automatically detects the processing mode per KN_BC:
- **Single-library**: Only 1 library maps to a KN_BC → Skera BAM passed directly as HIGH_CONF (no Bayesian assignment needed)
- **Multi-library**: 2+ libraries share a KN_BC → Full Bayesian assignment pipeline runs

Both modes can coexist in the same run. Mode can be manually overridden in the manifest.

### Input File Formats

**Manifest (`manifest.tab`)** — Tab-separated, 4–5 columns:
```
sample    KN_BC    bam_path    arrays_path    [mode]
```
- `mode` is optional; auto-detected if omitted

**Arrays file** — Tab-separated, no header:
```
library_name    KN_BC    barcode1    barcode2    ...
```
- Multiple barcodes per library are listed as additional columns
- Barcodes are shared strings (e.g., `bc1019`, `bc1020`)

### Assignment Logic (`assign_kinnex.py`)

The Bayesian scorer computes a score per library for each ZMW:

```
score = Σ(specific_barcodes * SPECIFIC_WEIGHT)
      + Σ(shared_barcodes * min(1.0/count, MAX_SHARED_WEIGHT))
      + Σ(discordant_barcodes * DISCORDANT_PENALTY)
```

Posteriors are computed via softmax over library scores.

**Default thresholds** (hardcoded constants, optimizable via `train_script/`):
```python
SPECIFIC_WEIGHT = 1.0
MAX_SHARED_WEIGHT = 0.2
DISCORDANT_PENALTY = -0.10
POSTERIOR_HIGH_CONF = 0.840
POSTERIOR_LOW_CONF = 0.50
MIN_OBS_HIGH_CONF = 3
MIN_OBS_LOW_CONF = 2
MIN_SPECIFIC_HIGH_CONF = 1
MIN_SPECIFIC_LOW_CONF  = 0
```

**Classification tiers**:
- `HIGH_CONF`: posterior ≥ 0.840 AND total observations ≥ 3 AND specific barcodes ≥ 1
- `LOW_CONF`: posterior ≥ 0.50 AND total observations ≥ 2
- `UNASSIGNED`: below thresholds

### Assignment Output Schema

The assignment TSV has a comment header block (`#`-prefixed) with provenance (git hash, run date, parameters), followed by tab-separated data:
```
ZMW | Assigned_Array | Classification | Top_Posterior | N_Observations |
Informative_Barcodes | Uninformative_Barcodes | Extraneous_Barcodes | Array_Kinnex | All_Barcodes
```

Use `utils.load_assignments()` or `utils.load_assignments_df()` to load — these skip the comment header automatically.

---

## Development Workflows

### Running the Production Pipeline

```bash
cd prod_script/
# Dry run
snakemake --configfile config.yaml -n

# Local execution
snakemake --configfile config.yaml --cores 8 --use-conda

# SGE cluster submission
snakemake --configfile config.yaml --use-conda --cluster "qsub ..." --jobs 50
```

### Running the Test Pipeline

```bash
# Step 1: Generate synthetic ground truth pools
cd test_script/
snakemake -s Snakefile_subsample --configfile config_subsample.yaml --cores 4 --use-conda

# Step 2: Run production pipeline on synthetic pools
cd ../prod_script/
snakemake --configfile config.yaml --cores 8 --use-conda

# Step 3: Validate accuracy against ground truth
cd ../test_script/
snakemake -s Snakefile_qc --configfile config_qc.yaml --cores 4 --use-conda
```

### Parameter Optimization

```bash
cd train_script/

# Optimize barcode weights
python optimize_barcode_weights.py --assignments <file> --lookup <dir>

# Optimize thresholds (multi-dimensional grid search)
python optimize_thresholds_v2.py --assignments <file> --lookup <dir>
```

### Config File Management

Config files contain local absolute paths and are hidden from git by default using `--assume-unchanged`. Use the helper script to manage them:

```bash
./commit_configs.sh --protect   # Default: hide configs from git (local paths safe)
./commit_configs.sh --upload    # Expose configs for committing template changes
./commit_configs.sh --restore   # Pull config defaults from git, discard local changes
```

**Important**: Always run `--protect` before committing to avoid leaking local paths.

---

## Shared Utilities (`utils.py`)

All scripts that consume manifests, arrays files, or assignments should import from `utils.py` in `prod_script/scripts/`. Key functions:

| Function | Purpose |
|---|---|
| `parse_manifest(path)` | Parse manifest.tab, resolve mode per KN_BC |
| `resolve_mode(mode, libraries)` | Validate/auto-detect single vs multi mode |
| `read_arrays(path, kn_bc)` | Parse arrays file, filter to a KN_BC |
| `libraries_for_knbc(arrays, kn_bc)` | Get library names for a given KN_BC |
| `load_assignments(path)` | Load assignment file as dict (skips headers) |
| `load_assignments_df(path)` | Load assignment file as DataFrame |
| `parse_assignment_header(path)` | Extract provenance/parameters from header |

**Note**: `test_script/` and `train_script/` do not yet import `utils.py` — this is a planned improvement (tracked in `REFACTOR_NOTES.md`).

---

## Code Conventions

### Python Style

- Standard library imports first, then third-party (pandas, numpy, pysam), then local
- Functions are generally short and single-purpose
- Argparse used for all CLI scripts
- Logging via `print()` to stderr (not Python `logging` module)
- ZMW name normalization: take only the first two slash-separated fields (`movie/zmw`)

### Snakemake Conventions

- Each tool (Skera, Lima, IsoSeq) runs in its own conda environment (`envs/`)
- Resources specified in `config.yaml` and referenced via `config["rule_name"]["threads"]`
- Checkpoints used to gate downstream rules on dynamically determined outputs (after split by library)
- Sentinel BAMs (header-only) created for zero-read confidence tiers to prevent missing file errors
- Shell commands in rules use `log:` directive and redirect stderr to log file

### Reproducibility

- Git hash recorded in assignment output headers automatically
- Run date and all scoring parameters written to assignment file comments
- Random seed fixed at 42 for subsampling
- Parameter provenance parseable by `parse_assignment_header()` for downstream scripts

### Barcode Terminology

The codebase uses two naming conventions (a rename is planned):
- Current code: `informative` / `uninformative` barcodes
- Intended new names: `specific` / `shared` barcodes

When reading or modifying code, be aware both terms refer to the same concepts.

---

## Known Issues & Planned Work

From `REFACTOR_NOTES.md`:

1. **Sankey plot broken** — `generate_sankey.py` has a data format mismatch; do not rely on it
2. **Rule granularity** — `merge_and_finalize` in `Snakefile_subsample` bundles 3 separate outputs into one rule; should be split
3. **Terminology rename** — `informative`→`specific`, `uninformative`→`shared` throughout scripts
4. **MIN_INF gate** — Add minimum specific-barcode observation requirement to assignment classifier
5. **Comment pass** — Many functions and non-obvious logic blocks need inline documentation
6. **Utils extension** — `test_script/` and `train_script/` should import from `utils.py` instead of duplicating logic
7. **commit_configs.sh undocumented** — Should be added to README

When working on these, keep changes focused and create separate commits per item.

---

## Output Directory Structure

Production pipeline outputs (relative to `prod_script/`):
```
results/
├── skera/{sample}/           # Skera S-reads BAMs
├── lima_pass1/{sample}/      # Lima pass 1 reports (barcode observation counts)
├── assigned/{sample}/        # Assignment TSV files (with provenance headers)
├── split/{sample}/           # Per-library HIGH_CONF, LOW_CONF, unassigned BAMs
├── lima_pass2/{sample}/{lib}/# Fully demultiplexed BAMs
├── refine/{sample}/{lib}/    # FLNC BAMs with polyA filter
└── qc/                       # TSV summaries, Excel workbook, Sankey HTML, posterior PNGs
```

---

## Dependencies

All tools run in isolated conda environments. No global installation needed except Snakemake.

| Tool | Environment | Purpose |
|---|---|---|
| Snakemake | (base) | Workflow orchestration |
| Python 3.x + pandas/numpy/pysam/sklearn | `envs/python.yaml` | Custom scripts |
| Skera | `envs/skera.yaml` | Array segmentation |
| Lima | `envs/lima.yaml` | Barcode detection/demux |
| IsoSeq | `envs/isoseq.yaml` | FLNC generation |
| matplotlib/seaborn/plotly | `envs/python.yaml` | Visualization |

---

## What AI Assistants Should Know

- **No unit tests exist** — validation is integration-based via synthetic ground truth
- **BAM file handling** uses `pysam`; always open with appropriate mode (`rb` for read, `wb` for write)
- **Assignment files** have comment-prefixed headers — always use `utils.load_assignments*()` to load them safely
- **Snakemake checkpoints** are critical for the split→lima_pass2 dependency chain; changes there require understanding checkpoint mechanics
- **The Sankey plot is broken** — do not attempt to use or fix it without understanding the full data flow
- **Config files** intentionally contain local paths and are git-ignored via `--assume-unchanged`; never commit them with real paths
- **Single-library mode** completely bypasses `assign_kinnex.py` — do not assume it always runs
- **Sentinel BAMs** (header-only) exist by design for zero-read tiers to prevent Snakemake missing-file failures

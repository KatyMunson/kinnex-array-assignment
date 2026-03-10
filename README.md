# Kinnex assignment Iso-Seq Pipeline - Determine Library from Data

Pipeline for processing PacBio Kinnex full-length RNA data. Designed to support data-based library assignment of pools of Kinnex libraries sharing a Kinnex barcode. Single-library and multi-library samples can be processed in the same run, with per-sample mode detection from the arrays file.

---

## Overview

Three separate pipelines work together:

| Pipeline | Snakefile | Folder | Purpose |
|---|---|---|---|
| **Production** | `Snakefile` | `prod_script/` | Full processing: Skera → (assignment) → Lima → Refine → QC |
| **Subsample** | `Snakefile_subsample` | `test_script/` | Creates synthetic test pools from pre-Skera BAMs with known ground truth |
| **QC** | `Snakefile_qc` | `test_script/` | Validates assignment accuracy against ground truth from subsampling |

The subsample and QC pipelines are used together to validate the classifier
and live in the same directory. The production pipeline runs separately.
Optimizer scripts for retraining the classifier live in `train_script/`.

---

## Dependencies

The following tools must be available on `PATH` (or via conda/envmodules —
see [Running](#running)):

- `skera` — PacBio Kinnex array segmentation
- `lima` — IsoSeq barcode demultiplexing
- `isoseq` — IsoSeq refine (FLNC generation)
- `samtools` — BAM operations

Python packages:

```
pandas
numpy
matplotlib
seaborn
scikit-learn
openpyxl
pysam
```

---

## Processing Modes

The production pipeline handles two processing modes per sample, determined
automatically from the arrays file.

### Single-library

One library maps to the KN_BC for this sample. No Lima pass 1, no Bayesian
assignment, no split step. The Skera output BAM is linked directly as
`{library}.highconf.bam`. HIGH_CONF is definitionally correct — there is
nothing to disambiguate.

### Multi-library

Two or more libraries share a KN_BC. Full Bayesian assignment pipeline:
Lima pass 1 (report only) → Bayesian ZMW classifier → split by library →
Lima pass 2 → IsoSeq refine.

From Lima pass 2 onward, single-library and multi-library samples are
processed identically. Both produce `highconf` outputs; multi-library
samples additionally produce `lowconf` outputs for ZMWs that did not meet
the high-confidence posterior threshold.

---

## Input Files

### Arrays file

Tab-separated, no header, one library per line:

```
{library_name}    {KN_BC}    {barcode1}    {barcode2}    ...
```

Example:
```
1KG_Pool076_KN    bcM0003    IsoSeqX_bc01_5p    IsoSeqX_bc02_5p    IsoSeqX_bc03_5p
1KG_Pool057_KN    bcM0003    IsoSeqX_bc04_5p    IsoSeqX_bc05_5p    IsoSeqX_bc06_5p
1KG_Pool038_KN    bcM0004    IsoSeqX_bc01_5p    IsoSeqX_bc02_5p    IsoSeqX_bc03_5p
```

A single master arrays file covering all KN_BCs on a SMRT Cell is
recommended. The pipeline subsets by KN_BC per sample automatically.

### Production manifest (`config.yaml` → `manifest`)

Tab-separated, no header, `#` lines ignored. Four required columns plus an
optional fifth:

```
{sample}    {KN_BC}    {path/to/skera_input.bam}    {path/to/arrays.txt}    [{mode}]
```

| Column | Description |
|---|---|
| `sample` | Unique identifier for this pool / KN_BC combination |
| `KN_BC` | Kinnex barcode (e.g. `bcM0003`) |
| `bam_path` | Path to pre-Skera CCS/HiFi BAM |
| `arrays_path` | Path to arrays file (may be a master file for the full SMRT Cell) |
| `mode` | **Optional.** `single` or `multi`. Overrides auto-detection (see below) |

#### Mode auto-detection

When the `mode` column is absent, the pipeline counts libraries sharing this
sample's KN_BC in the arrays file at startup:

- 1 library → **single-library** mode
- 2+ libraries → **multi-library** mode

The detected mode is printed to the log for every sample before any jobs
are submitted.

#### Mode override cross-validation

When a `mode` column is present:

| Declared | Libraries found | Behaviour |
|---|---|---|
| `single` | 1 | Proceeds as single-library |
| `single` | >1 | **Error** — library name is ambiguous; cannot proceed |
| `multi` | >1 | Proceeds as multi-library |
| `multi` | 1 | **Warning**, then proceeds through the assignment path as declared |

#### Mixed-mode example

Single-library and multi-library samples can coexist in one manifest:

```
# Two KN_BCs with 2 libraries each → auto-detected MULTI
1_C01_M0003    bcM0003    /path/to/Cell1.ccs.bam    arrays/Cell1_master.txt
1_C01_M0004    bcM0004    /path/to/Cell1.ccs.bam    arrays/Cell1_master.txt
# Two KN_BCs with 1 library each → auto-detected SINGLE
1_C01_M0001    bcM0001    /path/to/Cell1.ccs.bam    arrays/Cell1_master.txt
1_C01_M0002    bcM0002    /path/to/Cell1.ccs.bam    arrays/Cell1_master.txt
```

---

## Production Pipeline

### Setup

```
prod_script/
├── Snakefile
├── config.yaml
├── manifest.tab              ← your run manifest (not committed — see manifest.tab.example)
├── primers/
│   └── primers.fasta
├── adapters/
│   └── adapters.fasta
├── envs/
│   └── *.yaml                ← conda environment definitions
└── scripts/
    ├── utils.py              ← shared manifest parsing and mode resolution
    ├── assign_kinnex.py
    ├── split_skera_by_library.py
    ├── aggregate_flnc_qc.py
    ├── aggregate_pipeline_qc.py
    └── generate_sankey.py
```

### Configuration (`config.yaml`)

```yaml
manifest: manifest.tab
primers_fasta:  primers/primers.fasta
adapters_fasta: adapters/adapters.fasta

assign_script:       scripts/assign_kinnex.py
split_skera_script:  scripts/split_skera_by_library.py
aggregate_script:    scripts/aggregate_flnc_qc.py
pipeline_qc_script:  scripts/aggregate_pipeline_qc.py
sankey_script:       scripts/generate_sankey.py

refine:
  require_polya:    true
  min_polya_length: 20

resources:
  skera:     { threads: 8, mem: 1,  hrs: 6  }
  lima:      { threads: 8, mem: 1,  hrs: 12 }
  refine:    { threads: 8, mem: 1,  hrs: 8  }
  assign:    { threads: 2, mem: 6,  hrs: 8  }
  split:     { threads: 2, mem: 6,  hrs: 8  }
  aggregate: { threads: 1, mem: 4,  hrs: 1  }
  qc:        { threads: 1, mem: 8,  hrs: 2  }
```

`mem` is GB **per slot**. The cluster multiplies `mem × threads` for the
total memory request.

### Running

```bash
# Dry run
snakemake -s Snakefile --configfile config.yaml -n

# Local
snakemake -s Snakefile --configfile config.yaml --cores 8 --use-conda

# SGE cluster (--retries 3 for memory scaling)
snakemake -s Snakefile --configfile config.yaml \
    --cluster "qsub -pe serial {threads} \
               -l h_rt={resources.hrs}:00:00 \
               -l mfree={resources.mem}G \
               -o log/ -e log/ -cwd" \
    --jobs 48 --retries 3
```

### Pipeline steps

| Step | Rule | Samples | Notes |
|---|---|---|---|
| 1 | `skera_split` | All | Segments Kinnex arrays into S-reads |
| 2 | `lima_isoseq` | Multi-library only | Barcode detection, `--no-output`; report feeds classifier |
| 3 | `assign_arrays` | Multi-library only | Bayesian ZMW → library classification |
| 4a | `split_skera_by_library` | Multi-library only | Splits Skera BAM into per-library highconf/lowconf BAMs |
| 4b | `passthrough_single_library` | Single-library only | Hard-links Skera BAM directly as `{library}.highconf.bam` |
| 4c | `split_skera_collect` | All | Checkpoint: gates Lima pass 2 until all step 4 jobs finish |
| 5 | `lima_pass2` | All | Full demux per library/confidence, `--peek-guess` (checkpoint) |
| 6 | `isoseq_refine` | All | FLNC reads with polyA filtering |
| 7 | `aggregate_flnc_qc` | All | Summarises FLNC counts per sample/library/confidence |
| 8 | `pipeline_qc` | All | Excel QC report and intermediate TSVs |
| 9 | `sankey_plot` | All | Interactive HTML read-flow diagram |

### Output structure

```
results/
├── skera/
│   └── {sample}.skera.bam
├── lima/                                      ← multi-library only
│   └── {sample}/
│       └── *.lima.report
├── assigned/                                  ← multi-library only
│   └── {sample}.txt
├── split_skera/
│   └── {sample}/
│       ├── {library}.highconf.bam             ← all samples
│       ├── {library}.lowconf.bam              ← multi-library only
│       ├── unassigned.bam                     ← multi-library only
│       └── split_skera_qc.tsv
├── lima2/
│   └── {sample}/{library}/{highconf|lowconf}/
│       └── *.fl.*.bam
├── refine/
│   └── {sample}/{library}/{highconf|lowconf}/
│       └── *.flnc.bam
└── qc/
    ├── flnc_qc_summary.tsv
    ├── pipeline_qc.xlsx
    └── pipeline_read_flow.html
```

### Assignment classification thresholds (multi-library only)

| Classification | Posterior threshold | Min observations |
|---|---|---|
| `HIGH_CONF` | ≥ 0.840 | ≥ 3 segments |
| `LOW_CONF` | ≥ 0.500 | ≥ 2 segments |
| `UNASSIGNED` | below thresholds | — |

### QC report (`pipeline_qc.xlsx`)

The Excel workbook contains one tab per pipeline stage. A
`processing_mode` column (`SINGLE-LIBRARY` / `MULTI-LIBRARY`) appears
on the Summary, Skera, Split_Skera, Lima_Pass2, Lima2_Barcodes, and
Refine tabs, colour-coded green (single) and blue (multi). The Lima_Pass1,
Assign, and Assign_Library tabs contain only multi-library sample rows.

---

## Subsampling Pipeline

Creates synthetic test pools by subsampling pre-Skera (CCS/HiFi) BAMs from
individual libraries and merging them. Because the source of each read is
known, these pools have exact ground truth for QC validation.

### Setup

```
test_script/
├── Snakefile_subsample
├── Snakefile_qc
├── config_subsample.yaml
├── config_qc.yaml
├── pool.tsv                   ← pool metadata (not committed — see pool.tsv.example)
├── manifest_preskera.tsv      ← pre-Skera BAM paths per library (not committed)
└── scripts/
    └── visualize_posteriors.py
```

### Pre-Skera manifest (`config_subsample.yaml` → `manifest`)

Tab-separated with header:

```
Kinnex_Lib          PreSkera_BAM
1KG_Pool076_KN      /path/to/1KG_Pool076_KN.ccs.bam
1KG_Pool057_KN      /path/to/1KG_Pool057_KN.ccs.bam
```

These should be CCS/HiFi BAMs **before** Skera segmentation — one read per ZMW.

Pool/KN_BC combinations with only one library are automatically skipped
(no assignment test needed for single-library pools).

### `pool.tsv`

The central metadata file. Tab-separated with header:

| Column | Description |
|---|---|
| `Pool` | Pool identifier (e.g. `1`, `Pool076`) |
| `KN_BC` | Kinnex barcode ID (e.g. `bcM0003`) |
| `Kinnex_Lib` | Library name (e.g. `1KG_Pool076_KN`) |
| `bc_set` | Comma-separated IsoSeq barcodes for this library (e.g. `bc01,bc02,bc03`) |
| `SC` | Relative cell count for proportional subsampling |

Multiple libraries sharing the same `Pool` + `KN_BC` form one multiplexed pool.

### Configuration (`config_subsample.yaml`)

```yaml
pool_tsv:  "pool.tsv"
manifest:  "manifest_preskera.tsv"

# "proportional" uses SC values from pool.tsv to weight library contributions
# "even" splits reads equally among all libraries in the pool
pooling_mode: "proportional"

# Total reads per pool (split across constituent libraries)
target_total_reads: 1000000

random_seed: 42
```

### Running

```bash
snakemake -s Snakefile_subsample --configfile config_subsample.yaml --cores 4
```

### Output structure

```
results/
├── subsampled/          ← individual per-library BAMs (deleted after merge)
├── merged/
│   └── {pool}_KN_{knbc}.bam          ← synthetic pool BAM
├── lookup/
│   └── {pool}_KN_{knbc}_reads.txt    ← ground truth: read → library mapping
└── arrays/
    └── {pool}_KN_{knbc}.txt          ← arrays file for the production pipeline
```


### Next steps after subsampling

The subsampling pipeline produces everything the production pipeline needs
directly — no manual steps required:

- `results/merged/{pool}_KN_{knbc}.bam` — synthetic pre-Skera BAM
- `results/arrays/{pool}_KN_{knbc}.txt` — arrays file
- `results/manifest.tab` — production manifest, ready to use

`results/manifest.tab` is generated automatically with the correct
`bam_path` and `arrays_path` columns pointing at the merged BAMs and
arrays files. Copy or symlink it into `prod_script/` and run the
production pipeline as normal:

```bash
cp test_script/results/manifest.tab prod_script/manifest.tab
cd prod_script
snakemake -s Snakefile --configfile config.yaml -n
```

Because all pools from subsampling are multi-library by construction,
mode will auto-detect as `MULTI-LIBRARY` for every sample.

---

## QC Pipeline

Compares production pipeline assignment results against the ground truth
generated by the subsampling step. Runs only on multi-library
Pool/KN_BC combinations (single-library pools are skipped automatically).

### Configuration (`config_qc.yaml`)

```yaml
pool_tsv: "pool.tsv"

# Directory of assignment .txt files from the production pipeline
assignment_dir: "path/to/production/results/assigned"

# Directory of ground truth lookup files from Snakefile_subsample
lookup_dir: "results/lookup"

# Path to visualization script
visualize_script: "scripts/visualize_posteriors.py"

# Optional: enables FLNC-level analysis
# refine_dir: "path/to/production/results/refine"
```

### Running

```bash
snakemake -s Snakefile_qc --configfile config_qc.yaml --cores 4
```

### Output structure

```
results/qc/
├── per_pool/
│   └── {pool}_KN_{knbc}/
│       ├── summary.txt
│       ├── errors.txt
│       ├── threshold_analysis.csv
│       ├── distributions.png
│       ├── posterior_threshold_analysis.csv
│       ├── flnc_summary.txt             ← if refine_dir set
│       └── flnc_threshold_analysis.csv  ← if refine_dir set
└── aggregate/
    ├── HIGH_CONF_correct_table.tsv
    ├── HIGH_CONF_incorrect_table.tsv
    ├── LOW_CONF_correct_table.tsv
    ├── LOW_CONF_incorrect_table.tsv
    ├── UNASSIGNED_table.tsv
    └── TOTAL_FRACTION_table.tsv
```

The aggregate tables have pools as rows and KN_BCs as columns, showing the
fraction of ZMWs in each classification category.

### FLNC analysis

When `refine_dir` is configured, the QC pipeline adds FLNC-level analysis:
of the ZMWs misassigned at HIGH_CONF, how many actually produced FLNC reads
that would appear in downstream analysis? The `FLNC_Incorrect_HC` column in
`flnc_threshold_analysis.csv` is the key metric.

---

## Optimizer Scripts

Standalone scripts in `train_script/` for tuning classifier parameters on
validated data. Run manually after generating ground truth with the
subsampling pipeline.

### `optimize_barcode_weights.py`

Finds optimal `SPECIFIC_WEIGHT`, `MAX_SHARED_WEIGHT`, and
`DISCORDANT_PENALTY` values for the Bayesian scorer.

```bash
python train_script/optimize_barcode_weights.py \
    --arrays test_script/results/arrays/Pool1_KN_bcM0003.txt \
    --assignments prod_script/results/assigned/*.txt \
    --lookups test_script/results/lookup/*_reads.txt \
    --output optimized_weights.json
```

### `optimize_thresholds.py`

1D posterior threshold optimisation using ROC and precision-recall analysis.

```bash
python train_script/optimize_thresholds.py \
    --assignments prod_script/results/assigned/*.txt \
    --lookups test_script/results/lookup/*_reads.txt \
    --output threshold_recommendations.json
```

### `optimize_thresholds_v2.py`

Multi-dimensional grid search over posterior threshold, minimum observations,
and minimum informative barcodes. Computes the Pareto frontier.

```bash
python train_script/optimize_thresholds_v2.py \
    --assignments prod_script/results/assigned/*.txt \
    --lookups test_script/results/lookup/*_reads.txt \
    --output threshold_recommendations_v2.json
```

After running, update `POSTERIOR_HIGH_CONF`, `POSTERIOR_LOW_CONF`,
`SPECIFIC_WEIGHT`, `MAX_SHARED_WEIGHT`, and `DISCORDANT_PENALTY` in
`prod_script/scripts/assign_kinnex.py` accordingly.

---

## Pool Design Guidelines

The maximum allowable barcode overlap between libraries depends on pool size.
With 6 barcodes per library and the current scorer weights, the theoretical
ceiling on achievable posterior (and therefore whether HIGH_CONF calls are
possible at all) is:

| Pool size | Max shared BCs | Min unique BCs | Max achievable posterior |
|---|---|---|---|
| 2 libraries | 4 | 2 | 0.900 |
| 3 libraries | 3 | 3 | 0.931 |
| 4 libraries | 3 | 3 | 0.900 |
| 5 libraries | 3 | 3 | 0.871 |
| 6 libraries | 3 | 3 | 0.844 |

The production rule of **max 2 shared barcodes** is conservative but
well-justified — it provides comfortable headroom above the theoretical
floor and accounts for real-world noise that pushes observed posteriors
below the theoretical ceiling.

---

## Shared Utilities (`scripts/utils.py`)

`utils.py` provides shared logic imported by `Snakefile`,
`aggregate_pipeline_qc.py`, and the training scripts:

- `parse_manifest()` — parses the 4- or 5-column manifest and resolves
  processing mode for every sample
- `resolve_mode()` — auto-detects or validates declared single/multi mode
  against the arrays file, with appropriate warnings and errors
- `read_arrays()` / `libraries_for_knbc()` — arrays file parsing
- `load_assignments()` / `load_assignments_df()` — loads assignment files,
  correctly skipping provenance `#` comment headers
- `parse_assignment_header()` — extracts scoring parameters from assignment
  file headers for training script validation

---

## Typical Workflow

```
1. Prepare pool.tsv with library and barcode metadata

2. [Recommended for new pool designs — validate before running on real data]
   a. Prepare pre-Skera manifest and run Snakefile_subsample
   b. Run Skera + production pipeline on merged BAMs
   c. Run Snakefile_qc to validate accuracy
   d. Adjust thresholds or weights if needed using train_script/

3. Prepare production manifest (single master arrays file recommended)
   and run Snakefile on real data

4. Review results/qc/pipeline_qc.xlsx and pipeline_read_flow.html
```

---

## Acknowledgements

Pipeline development was assisted by Claude (Anthropic). All scientific
decisions, validation, and code review were performed by the authors.

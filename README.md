# Kinnex Array Assignment Pipeline

Pipeline for demultiplexing PacBio Kinnex (IsoSeq) arrays from multiplexed single-cell RNA-seq experiments. Assigns ZMWs to their library of origin using a Bayesian classifier, then processes through lima and isoseq refine to produce FLNC reads per library.

---

## Overview

Three separate pipelines work together:

| Pipeline | Snakefile | Folder | Purpose |
|---|---|---|---|
| **Production** | `Snakefile` | `prod_script/` | Full processing: Skera → Lima → Assignment → Split → Lima2 → Refine → QC |
| **Subsample** | `Snakefile_subsample` | `test_script/` | Creates synthetic test pools from pre-Skera BAMs with known ground truth |
| **QC** | `Snakefile_qc` | `test_script/` | Validates assignment accuracy against ground truth from subsampling |

The subsample and QC pipelines are used together to validate the classifier and live in the same directory. The production pipeline runs separately in its own directory. Optimizer scripts for retraining the classifier live in `train_script/`.

---

## Dependencies

The following tools must be available on `PATH`:

- `skera` — PacBio Kinnex array segmentation
- `lima` — IsoSeq barcode demultiplexing  
- `isoseq` — IsoSeq refine (FLNC generation)
- `samtools` — BAM indexing and merging
- `pysam` — Python BAM handling (for split scripts)

Python packages:

```
pandas
numpy
matplotlib
seaborn
scikit-learn
openpyxl
```

---

## Input Files

### `pool.tsv`

The central metadata file used by all three pipelines. Tab-separated with the following required columns:

| Column | Description |
|---|---|
| `Pool` | Pool identifier (e.g. `1`, `Pool076`) |
| `KN_BC` | Kinnex barcode ID (e.g. `bcM0003`) |
| `Kinnex_Lib` | Library name (e.g. `1KG_Pool076_KN`) |
| `bc_set` | Comma-separated IsoSeq barcodes for this library (e.g. `bc01,bc02,bc03`) |
| `SC` | Relative cell count used for proportional subsampling |

Multiple libraries sharing the same `Pool` + `KN_BC` combination form one multiplexed array pool.

### Production manifest (`config.yaml` → `manifest`)

Tab-separated, no header, four columns:

```
{sample}    {KN_BC}    {path/to/skera_input.bam}    {path/to/arrays.txt}
```

### Arrays file (`arrays.txt`)

One library per line, fully tab-separated:

```
{library_name}    {KN_BC}    {barcode1}    {barcode2}    ...
```

Example:
```
1KG_Pool076_KN    bcM0003    IsoSeqX_bc01_5p    IsoSeqX_bc02_5p    IsoSeqX_bc03_5p
1KG_Pool057_KN    bcM0003    IsoSeqX_bc04_5p    IsoSeqX_bc05_5p    IsoSeqX_bc06_5p
```

This file is produced automatically by `Snakefile_subsample`. For production runs, create it manually or adapt from your pool metadata.

---

## Production Pipeline

Processes real multiplexed Kinnex data end-to-end.

### Setup

The production pipeline lives in `prod_script/` in the repository:

```
prod_script/
├── Snakefile
├── config.yaml
├── pool.tsv               ← your pool metadata (not committed — see pool.tsv.example)
├── manifest.tab           ← your run manifest (not committed — see manifest.tab.example)
├── primers/
│   └── primers.fasta
├── adapters/
│   └── adapters.fasta
├── envs/
│   └── *.yaml             ← conda environment definitions
└── scripts/
    ├── assign_kinnex.py
    ├── split_skera_by_library.py
    ├── aggregate_flnc_qc.py
    ├── aggregate_pipeline_qc.py
    └── generate_sankey.py
```

### Configuration (`config.yaml`)

```yaml
manifest: manifest.tab
primers_fasta: primers/primers.fasta
adapters_fasta: adapters/adapters.fasta

assign_script: scripts/assign_kinnex.py
split_skera_script: scripts/split_skera_by_library.py
aggregate_script: scripts/aggregate_flnc_qc.py
pipeline_qc_script: scripts/aggregate_pipeline_qc.py
sankey_script: scripts/generate_sankey.py

lima:
  peek_guess: false

refine:
  require_polya: true
  min_polya_length: 20

resources:
  skera:    { threads: 8, mem: 6,  hrs: 6  }
  lima:     { threads: 8, mem: 8,  hrs: 12 }
  refine:   { threads: 8, mem: 6,  hrs: 8  }
  assign:   { threads: 2, mem: 12, hrs: 8  }
  split:    { threads: 2, mem: 12, hrs: 8  }
  aggregate:{ threads: 1, mem: 4,  hrs: 1  }
  qc:       { threads: 1, mem: 8,  hrs: 2  }
```

### Running

```bash
# Dry run
snakemake -s Snakefile --configfile config.yaml -n

# Run locally
snakemake -s Snakefile --configfile config.yaml --cores 8

# Run on SGE cluster
snakemake -s Snakefile --configfile config.yaml \
    --cluster "qsub -pe smp {threads} -l h_vmem={resources.mem}G -l h_rt={resources.hrs}:00:00" \
    --jobs 50
```

### Output Structure

```
results/
├── skera/
│   └── {sample}.skera.bam
├── lima/
│   └── {sample}/                          ← lima pass 1 (report only)
├── assigned/
│   └── {sample}.txt                       ← ZMW → library assignments
├── split_skera/
│   └── {sample}/
│       └── {library}.{highconf|lowconf}.bam
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

### Pipeline Steps

1. **skera split** — Segments concatenated Kinnex arrays into individual reads
2. **lima (pass 1)** — Barcode detection only (`--no-output`); produces the report used for assignment
3. **assign_kinnex** — Classifies each ZMW using a Bayesian scorer. For each ZMW, scores all libraries sharing the same Kinnex barcode using specific (library-unique) and shared barcode observations weighted by specificity. Outputs `HIGH_CONF`, `LOW_CONF`, or `UNASSIGNED`
4. **split_skera_by_library** — Splits the Skera BAM into per-library highconf/lowconf BAMs
5. **lima (pass 2)** — Full demultiplexing per library with `--peek-guess`
6. **isoseq refine** — Generates FLNC reads with polyA filtering
7. **aggregate_flnc_qc** — Summarizes FLNC counts across libraries
8. **pipeline_qc** — Produces Excel QC report and Sankey diagram of read flow

### Assignment Classification Thresholds

| Classification | Posterior threshold | Min observations |
|---|---|---|
| `HIGH_CONF` | ≥ 0.840 | ≥ 3 segments |
| `LOW_CONF` | ≥ 0.500 | ≥ 2 segments |
| `UNASSIGNED` | below thresholds | — |

---

## Subsampling Pipeline

Creates synthetic test pools by subsampling pre-Skera (CCS/HiFi) BAMs from individual libraries and merging them together. Because the source of each read is known, these pools have exact ground truth for QC validation.

### Setup

The subsampling and QC pipelines share the `test_script/` directory in the repository:

```
test_script/
├── Snakefile_subsample
├── Snakefile_qc
├── config_subsample.yaml
├── config_qc.yaml
├── pool.tsv                   ← same pool.tsv as production (not committed)
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

### Configuration (`config_subsample.yaml`)

```yaml
pool_tsv: "pool.tsv"
manifest: "manifest_preskera.tsv"

# "proportional" uses SC values from pool.tsv to weight library contributions
# "even" splits reads equally among all libraries in the pool
pooling_mode: "proportional"

# Total reads per pool (split across constituent libraries)
target_total_reads: 1000000

random_seed: 42
```

For quick testing, set `target_total_reads: 10000`.

### Running

```bash
snakemake -s Snakefile_subsample --configfile config_subsample.yaml --cores 4
```

### Output Structure

```
results/
├── subsampled/          ← individual per-library BAMs (deleted after merge)
├── merged/
│   └── {pool}_KN_{knbc}.bam     ← synthetic pool BAM to run through production
├── lookup/
│   └── {pool}_KN_{knbc}_reads.txt   ← ground truth: read → library mapping
└── arrays/
    └── {pool}_KN_{knbc}.txt         ← arrays file for production pipeline
```

### Next Steps After Subsampling

The merged BAMs must be processed through Skera and the production pipeline before QC can be run. This is done manually — run Skera on each merged BAM, then run the production pipeline using the generated `results/arrays/*.txt` files as the arrays input.

```bash
# Run Skera on merged BAMs (adjust paths and adapter file as needed)
for bam in results/merged/*.bam; do
    name=$(basename $bam .bam)
    skera split $bam adapters/adapters.fasta results/skera/${name}.skera.bam
done

# Then run the production pipeline on the skera outputs,
# pointing the manifest at results/arrays/ for the arrays files
```

---

## QC Pipeline

Compares production pipeline assignment results against the ground truth generated by the subsampling step.

### Setup

Run from `test_script/` alongside the subsampling pipeline (same directory):

```
test_script/
├── Snakefile_qc
├── config_qc.yaml
├── pool.tsv
└── scripts/
    └── visualize_posteriors.py
```

### Configuration (`config_qc.yaml`)

```yaml
pool_tsv: "pool.tsv"

# Directory of assignment .txt files from the production pipeline
assignment_dir: "path/to/production/results/assigned"

# Directory of ground truth lookup files from Snakefile_subsample
lookup_dir: "results/lookup"

# Path to visualization script
visualize_script: "scripts/visualize_posteriors.py"

# Optional: enables FLNC-level analysis showing how many misassigned ZMWs
# actually produced transcripts. Requires refine output from production run.
# refine_dir: "path/to/production/results/refine"
```

### Running

```bash
snakemake -s Snakefile_qc --configfile config_qc.yaml --cores 4
```

### Output Structure

```
results/qc/
├── per_pool/
│   └── {pool}_KN_{knbc}/
│       ├── summary.txt                       ← ZMW and segment-level accuracy (uses ground truth)
│       ├── errors.txt                        ← Detailed list of misassigned ZMWs
│       ├── threshold_analysis.csv            ← Accuracy at alternative thresholds (uses ground truth)
│       ├── distributions.png                 ← Posterior distribution plots
│       ├── posterior_threshold_analysis.csv  ← Threshold analysis from posterior distributions alone
│       ├── flnc_summary.txt                  ← (if refine_dir set) FLNC outcome counts
│       └── flnc_threshold_analysis.csv       ← (if refine_dir set) FLNC-aware thresholds
└── aggregate/
    ├── HIGH_CONF_correct_table.tsv
    ├── HIGH_CONF_incorrect_table.tsv
    ├── LOW_CONF_correct_table.tsv
    ├── LOW_CONF_incorrect_table.tsv
    ├── UNASSIGNED_table.tsv
    └── TOTAL_FRACTION_table.tsv
```

The aggregate tables have pools as rows and Kinnex barcodes as columns, showing the fraction of ZMWs in each assignment category. These are the primary summary for comparing accuracy across pools.

### FLNC Analysis

When `refine_dir` is configured, the QC pipeline adds FLNC-level analysis. This answers the question: of the ZMWs that were misassigned at HIGH_CONF, how many actually produced FLNC reads (i.e., full-length transcripts that would appear in downstream analysis)?

`flnc_threshold_analysis.csv` tests alternative posterior and observation thresholds and reports both ZMW-level accuracy and the number of misassigned FLNC segments that would result from each threshold. The `FLNC_Incorrect_HC` column is the key metric — it shows actual misassigned transcripts in the high-quality output.

---

## Optimizer Scripts

Two standalone scripts (plus an extended version) live in `train_script/` and can be used to tune the classifier parameters on validated data. Run manually, outside of Snakemake, after generating ground truth data with the subsampling pipeline.

### `optimize_barcode_weights.py`

Uses machine learning to find the optimal informative/uninformative/extraneous weights for the Bayesian scorer, training on pools with known ground truth.

```bash
python train_script/optimize_barcode_weights.py \
    --arrays test_script/results/arrays/Pool1_KN_bcM0003.txt \
    --assignments prod_script/results/assigned/*.txt \
    --lookups test_script/results/lookup/*_reads.txt \
    --output optimized_weights.json
```

### `optimize_thresholds.py`

1D posterior threshold optimisation. Uses ROC and precision-recall analysis to recommend `HIGH_CONF` and `LOW_CONF` posterior thresholds given a target accuracy or yield.

```bash
python train_script/optimize_thresholds.py \
    --assignments prod_script/results/assigned/*.txt \
    --lookups test_script/results/lookup/*_reads.txt \
    --output threshold_recommendations.json
```

### `optimize_thresholds_v2.py`

Multi-dimensional grid search over posterior threshold, minimum observations, and minimum informative barcodes. Computes the Pareto frontier and produces heatmap and numerical table outputs.

```bash
python train_script/optimize_thresholds_v2.py \
    --assignments prod_script/results/assigned/*.txt \
    --lookups test_script/results/lookup/*_reads.txt \
    --output threshold_recommendations_v2.json
```

After running these, update `INF_WEIGHT`, `MAX_UNINF_WEIGHT`, `EXTRANEOUS_PENALTY`, `POSTERIOR_HIGH_CONF`, and `POSTERIOR_LOW_CONF` in `prod_script/scripts/assign_kinnex.py` accordingly.

---

## Typical Workflow

```
1. Prepare pool.tsv with library and barcode metadata

2. [If validating first — recommended for new pool designs]
   a. Prepare pre-Skera manifest and run Snakefile_subsample
   b. Run Skera + production pipeline on merged BAMs
   c. Run Snakefile_qc to validate accuracy
   d. Adjust thresholds or weights if needed

3. Prepare production manifest and run Snakefile on real data

4. Review results/qc/pipeline_qc.xlsx and pipeline_read_flow.html
```

---

## Acknowledgements

Pipeline development was assisted by Claude (Anthropic). All scientific
decisions, validation, and code review were performed by the authors.

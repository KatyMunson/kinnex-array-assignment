# Conda Environments

This folder contains conda environment definitions for the pipeline tools.
Each PacBio tool has its own environment to avoid dependency conflicts.

| File | Used by | Tools |
|------|---------|-------|
| `skera.yaml` | Snakefile (skera rules) | skera |
| `lima.yaml` | Snakefile (lima pass 1 & 2 rules) | lima |
| `isoseq.yaml` | Snakefile (refine rules) | isoseq3 |
| `python.yaml` | All scripts | pandas, pysam, numpy, scikit-learn, openpyxl, matplotlib, seaborn, plotly |

## Creating environments

```bash
conda env create -f envs/skera.yaml -n kinnex-skera
conda env create -f envs/lima.yaml -n kinnex-lima
conda env create -f envs/isoseq.yaml -n kinnex-isoseq
conda env create -f envs/python.yaml -n kinnex-python
```

## Notes

- The PacBio tools (skera, lima, isoseq) are available via bioconda and may
  also be pre-installed on your cluster — check with `module avail` or
  `which skera` before creating new environments.
- `python.yaml` covers all the pipeline scripts (assign_kinnex.py,
  aggregate_pipeline_qc.py, split_skera_by_library.py, optimizer scripts, etc.)
- These environments are not currently declared as `conda:` directives in the
  Snakefiles — see REFACTOR_NOTES.md for the planned improvement.

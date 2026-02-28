#!/usr/bin/env python3
"""
aggregate_flnc_qc.py

Aggregates FLNC poly-A read counts from isoseq refine JSON filter reports
into a summary TSV with one row per sample/library/confidence tier.

Usage: python aggregate_flnc_qc.py [refine_dir] [output_tsv]
  refine_dir:  path to results/refine/  (default: results/refine)
  output_tsv:  output file path         (default: results/qc/flnc_qc_summary.tsv)

Author:  KM
Created: 2026-02
"""
import sys
import json
import glob
import re
from pathlib import Path
from collections import defaultdict

refine_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/refine")
output_tsv = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results/flnc_qc_summary.tsv")

BARCODES = [f"bc{i:02d}" for i in range(1, 13)]

# counts[sample][library][confidence][barcode] = num_reads_flnc_polya
counts = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))

json_files = glob.glob(str(refine_dir / "*" / "*" / "*" / "*.flnc.filter_summary.report.json"))
if not json_files:
    print(f"WARNING: No JSON files found under {refine_dir} — writing empty summary", file=sys.stderr)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, "w") as out:
        out.write("sample\tlibrary\tconfidence\t" + "\t".join(BARCODES) + "\ttotal_flnc\n")
    print(f"Written to {output_tsv} (empty — no refine outputs found)", file=sys.stderr)
    sys.exit(0)

print(f"Found {len(json_files)} JSON files", file=sys.stderr)

for json_path in sorted(json_files):
    p = Path(json_path)
    confidence = p.parts[-2]
    library    = p.parts[-3]
    sample     = p.parts[-4]

    match = re.search(r'_(bc\d{2})_', p.name)
    if not match:
        print(f"WARNING: Could not extract barcode from {p.name}, skipping", file=sys.stderr)
        continue
    barcode = match.group(1)

    with open(p) as f:
        data = json.load(f)

    for attr in data.get("attributes", []):
        if attr.get("id") == "num_reads_flnc_polya":
            counts[sample][library][confidence][barcode] += attr["value"]
            break
    else:
        print(f"WARNING: num_reads_flnc_polya not found in {p.name}, skipping", file=sys.stderr)

# Write output
output_tsv.parent.mkdir(parents=True, exist_ok=True)
with open(output_tsv, "w") as out:
    out.write("sample\tlibrary\tconfidence\t" + "\t".join(BARCODES) + "\ttotal_flnc\n")
    for sample in sorted(counts):
        for library in sorted(counts[sample]):
            for confidence in sorted(counts[sample][library]):
                bc_counts = counts[sample][library][confidence]
                values = [bc_counts.get(bc, 0) for bc in BARCODES]
                total  = sum(values)
                out.write(f"{sample}\t{library}\t{confidence}\t" +
                          "\t".join(str(v) for v in values) +
                          f"\t{total}\n")

print(f"Written to {output_tsv}", file=sys.stderr)

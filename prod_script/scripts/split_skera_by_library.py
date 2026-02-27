#!/usr/bin/env python3
"""
split_skera_by_library.py

Splits a Skera BAM by library assignment. Reads the assign_kinnex.py output
and writes one BAM per library (highconf and lowconf separately), plus an
unassigned BAM for ZMWs that did not receive a HIGH_CONF or LOW_CONF call.
Also writes a per-library read count QC TSV.

Usage:
    python split_skera_by_library.py <skera.bam> <assign.txt> <out_dir>

Outputs:
    <out_dir>/<library>.highconf.bam  — one per library, HIGH_CONF calls only
    <out_dir>/<library>.lowconf.bam   — one per library, LOW_CONF calls only
    <out_dir>/unassigned.bam          — UNASSIGNED ZMWs
    <out_dir>/split_skera_qc.tsv      — read counts per library and tier

Author:  KM
Created: 2026-02
"""

import sys
import pysam
from collections import defaultdict
from pathlib import Path

skera_bam  = sys.argv[1]
assign_file = sys.argv[2]
out_dir    = Path(sys.argv[3])
out_dir.mkdir(parents=True, exist_ok=True)

# ------------------
# Load assignments
# ------------------
zmw2high = {}
zmw2low  = {}

with open(assign_file) as f:
    next(f)  # skip header
    for line in f:
        parts = line.rstrip().split("\t")
        zmw, lib, cls = parts[0], parts[1], parts[2]

        if cls == "HIGH_CONF":
            zmw2high[zmw] = lib
        elif cls == "LOW_CONF":
            zmw2low[zmw] = lib

# ------------------
# Stream skera BAM and write per-library BAMs
# ------------------
writers = {}
counts = defaultdict(int)

def get_writer(lib, template):
    if lib not in writers:
        outbam = out_dir / f"{lib}.highconf.bam"
        writers[lib] = pysam.AlignmentFile(str(outbam), "wb", template=template)
    return writers[lib]
    
def get_low_writer(lib, template):
    key = f"{lib}.lowconf"
    if key not in writers:
        outbam = out_dir / f"{lib}.lowconf.bam"
        writers[key] = pysam.AlignmentFile(str(outbam), "wb", template=template)
    return writers[key]

with pysam.AlignmentFile(skera_bam, "rb", check_sq=False) as bam:
    # Pre-open unassigned writer
    unassigned_writer = pysam.AlignmentFile(
        str(out_dir / "unassigned.bam"), "wb", template=bam
    )
    for read in bam:
        zmw = "/".join(read.query_name.split("/")[:2])
        if zmw in zmw2high:
            lib = zmw2high[zmw]
            get_writer(lib, bam).write(read)
            counts[f"{lib}.HIGH_CONF"] += 1

        elif zmw in zmw2low:
            lib = zmw2low[zmw]
            get_low_writer(lib, bam).write(read)
            counts[f"{lib}.LOW_CONF"] += 1

        else:
            unassigned_writer.write(read)
            counts["unassigned"] += 1

for w in writers.values():
    w.close()
unassigned_writer.close()

# ------------------
# QC
# ------------------
qc_file = out_dir / "split_skera_qc.tsv"
with open(qc_file, "w") as out:
    out.write("library\tread_count\n")
    for lib, n in sorted(counts.items()):
        out.write(f"{lib}\t{n}\n")

total = sum(counts.values())
assigned = total - counts.get("unassigned", 0)
print(f"Total reads:    {total}")
print(f"Assigned:       {assigned} ({100*assigned/total:.2f}%)")
print(f"Unassigned:     {counts.get('unassigned', 0)} ({100*counts.get('unassigned',0)/total:.2f}%)")
print(f"Libraries:      {sorted(k for k in counts if k != 'unassigned')}")
print(f"QC written to:  {qc_file}")

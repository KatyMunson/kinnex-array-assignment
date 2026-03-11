#!/usr/bin/env python3
"""
assign_kinnex.py

Single-tier Bayesian assignment of ZMWs to Kinnex array libraries.
Classifies each ZMW by scoring it against all libraries sharing its expected
Kinnex barcode, using specific (library-unique) and shared barcode
observations weighted by specificity.

Scoring weights were optimised by ML training on ground-truth pools (Feb 2026).
Run optimize_barcode_weights.py to retrain on new data.

Author:  KM
Created: 2026-02
"""

import pandas as pd
from collections import defaultdict
import math
import argparse
from multiprocessing import Pool
import os
import subprocess
import datetime

# =======================
# PROVENANCE
# =======================
def get_git_provenance():
    """
    Return a short git hash for the current commit, with '-dirty' appended if
    there are uncommitted changes in the working tree. Returns 'unknown' if
    the script is not inside a git repository or git is unavailable.

    Examples:
        'a3f2c1b'        — clean commit, fully reproducible
        'a3f2c1b-dirty'  — uncommitted changes present; parameters are
                           recorded but exact code state is not pinned
        'unknown'        — not in a git repo or git not on PATH
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        git_hash = subprocess.run(
            ['git', 'rev-parse', '--short', 'HEAD'],
            capture_output=True, text=True, cwd=script_dir, timeout=5
        ).stdout.strip()
        if not git_hash:
            return 'unknown'
        dirty = subprocess.run(
            ['git', 'status', '--porcelain'],
            capture_output=True, text=True, cwd=script_dir, timeout=5
        ).stdout.strip()
        return f"{git_hash}-dirty" if dirty else git_hash
    except Exception:
        return 'unknown'

# =======================
# DEFAULT PARAMETERS
# =======================
SPECIFIC_WEIGHT = 1.0
MAX_SHARED_WEIGHT = 0.2
DISCORDANT_PENALTY = -0.10
POSTERIOR_HIGH_CONF = 0.840
POSTERIOR_LOW_CONF = 0.50
MIN_OBS_HIGH_CONF = 3
MIN_OBS_LOW_CONF = 2
MIN_SPECIFIC_HIGH_CONF = 1   # require at least 1 specific barcode for HIGH_CONF
MIN_SPECIFIC_LOW_CONF  = 0   # no specific-barcode floor for LOW_CONF

# =======================
# INPUT FUNCTIONS
# =======================
def read_arrays(arrays_file):
    arrays = {}
    with open(arrays_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            arrays[parts[0]] = {'kinnex': parts[1], 'barcodes': set(parts[2:])}
    return arrays

def read_lima(lima_file, limit=None):
    # limit applies to rows, not ZMWs; useful for quick testing on large reports
    zmw_dict = defaultdict(list)
    df = pd.read_csv(lima_file, sep='\t', dtype=str, nrows=limit)
    for _, row in df.iterrows():
        bc = row.get('IdxFirstNamed', '')
        if bc.endswith('_5p'):
            zmw_dict[row['ZMW']].append(bc)
        else:
            bc2 = row.get('IdxCombinedNamed', '')
            if bc2.endswith('_5p'):
                zmw_dict[row['ZMW']].append(bc2)
    return zmw_dict

# =======================
# CLASSIFICATION LOGIC
# =======================
def classify_single_tier(barcodes, arrays_subset, all_arrays):
    if not arrays_subset:
        return None, "UNASSIGNED", 0.0, {}
    
    array_names = list(arrays_subset.keys())
    barcode_array_map = defaultdict(list)
    for arr_name, arr_info in all_arrays.items():
        for bc in arr_info['barcodes']:
            barcode_array_map[bc].append(arr_name)

    array_scores = {}
    possible_arrays = set(array_names)
    n_obs = len(barcodes)
    
    for arr_name in array_names:
        score, specific_count, shared_count, discordant_count = 0.0, 0, 0, 0
        arr_barcodes = arrays_subset[arr_name]['barcodes']
        for bc in barcodes:
            if bc in arr_barcodes:
                arrays_with_bc = set(barcode_array_map[bc]) & possible_arrays
                if len(arrays_with_bc) == 1:
                    score += SPECIFIC_WEIGHT
                    specific_count += 1
                else:
                    weight = min(1.0 / len(arrays_with_bc), MAX_SHARED_WEIGHT)
                    score += weight
                    shared_count += 1
            else:
                discordant_count += 1
                score += DISCORDANT_PENALTY

        array_scores[arr_name] = {'score': score, 'specific': specific_count, 'shared': shared_count,
                                  'discordant': discordant_count, 'n_obs': n_obs, 'kinnex': arrays_subset[arr_name]['kinnex']}

    scores = [v['score'] for v in array_scores.values()]
    max_score = max(scores)
    exp_scores = [math.exp(s - max_score) for s in scores]
    total = sum(exp_scores)
    
    for arr_name, exp_s in zip(array_names, exp_scores):
        array_scores[arr_name]['posterior'] = exp_s / total

    best_array = max(array_scores, key=lambda k: array_scores[k]['posterior'])
    res = array_scores[best_array]
    posterior = res['posterior']
    
    if (n_obs >= MIN_OBS_HIGH_CONF and posterior >= POSTERIOR_HIGH_CONF
            and res['specific'] >= MIN_SPECIFIC_HIGH_CONF):
        classification = "HIGH_CONF"
    elif (n_obs >= MIN_OBS_LOW_CONF and posterior >= POSTERIOR_LOW_CONF
            and res['specific'] >= MIN_SPECIFIC_LOW_CONF):
        classification = "LOW_CONF"
    else:
        classification, best_array, posterior = "UNASSIGNED", None, 0.0

    return best_array, classification, posterior, array_scores

# =======================
# PARALLEL HELPERS
# =======================
_worker_arrays = None
_worker_assigned_kinnex = None

def _worker_init(arrays, assigned_kinnex):
    global _worker_arrays, _worker_assigned_kinnex
    _worker_arrays = arrays
    _worker_assigned_kinnex = assigned_kinnex

def _worker_classify(zmw_item):
    zmw, barcodes = zmw_item
    kinnex_arrays = {name: info for name, info in _worker_arrays.items()
                     if info['kinnex'] == _worker_assigned_kinnex}
    best_array, classification, posterior, scores = \
        classify_single_tier(barcodes, kinnex_arrays, _worker_arrays)

    summary = {'specific': scores[best_array]['specific'], 'shared': scores[best_array]['shared'],
               'discordant': scores[best_array]['discordant'], 'kinnex': scores[best_array]['kinnex']} \
              if best_array else \
              {'specific': 0, 'shared': 0, 'discordant': 0, 'kinnex': 'None'}

    return (zmw, barcodes, best_array, classification, posterior, summary)

# =======================
# MAIN
# =======================
def main():
    parser = argparse.ArgumentParser(description="Single-tier ZMW assignment")
    parser.add_argument("arrays_file")
    parser.add_argument("lima_file")
    parser.add_argument("kinnex_barcode")
    parser.add_argument("out_file")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    arrays = read_arrays(args.arrays_file)
    zmw_dict = read_lima(args.lima_file, limit=args.limit)
    zmw_items = list(zmw_dict.items())
    n_threads = args.threads if args.threads > 0 else os.cpu_count()

    if n_threads > 1:
        with Pool(processes=n_threads, initializer=_worker_init,
                  initargs=(arrays, args.kinnex_barcode)) as pool:
            results = pool.map(_worker_classify, zmw_items, chunksize=100)
    else:
        _worker_init(arrays, args.kinnex_barcode)
        results = [_worker_classify(item) for item in zmw_items]

    with open(args.out_file, "w") as out_f:
        # Write provenance header so parameters and code version are recorded
        # alongside the data. Lines starting with # are skipped by downstream
        # scripts that use pandas read_csv with comment='#'.
        run_date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        git = get_git_provenance()
        out_f.write(f"# assign_kinnex.py | {run_date} | git: {git}\n")
        out_f.write(f"# SPECIFIC_WEIGHT={SPECIFIC_WEIGHT}  MAX_SHARED_WEIGHT={MAX_SHARED_WEIGHT}  DISCORDANT_PENALTY={DISCORDANT_PENALTY}\n")
        out_f.write(f"# POSTERIOR_HIGH_CONF={POSTERIOR_HIGH_CONF}  POSTERIOR_LOW_CONF={POSTERIOR_LOW_CONF}\n")
        out_f.write(f"# MIN_OBS_HIGH_CONF={MIN_OBS_HIGH_CONF}  MIN_OBS_LOW_CONF={MIN_OBS_LOW_CONF}\n")
        out_f.write(f"# MIN_SPECIFIC_HIGH_CONF={MIN_SPECIFIC_HIGH_CONF}  MIN_SPECIFIC_LOW_CONF={MIN_SPECIFIC_LOW_CONF}\n")
        out_f.write("ZMW\tAssigned_Array\tClassification\tTop_Posterior\tN_Observations\tSpecific_Barcodes\tShared_Barcodes\tDiscordant_Barcodes\tArray_Kinnex\tAll_Barcodes\n")
        for zmw, barcodes, best_array, classification, posterior, summary in results:
            out_f.write(f"{zmw}\t{best_array}\t{classification}\t{posterior:.3f}\t{len(barcodes)}\t{summary['specific']}\t{summary['shared']}\t{summary['discordant']}\t{summary['kinnex']}\t{','.join(barcodes)}\n")

if __name__ == "__main__":
    main()

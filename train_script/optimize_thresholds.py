#!/usr/bin/env python3
"""
optimize_thresholds.py

1D posterior threshold optimisation for HIGH_CONF and LOW_CONF classification.
Sweeps posterior thresholds and recommends the operating point that maximises
yield subject to a target accuracy. Produces plots and a JSON recommendations
file.

See optimize_thresholds_v2.py for multi-dimensional search over posterior,
minimum observations, and minimum specific barcode counts.

Usage:
    python optimize_thresholds.py \\
        --assignments results/assigned/*.txt \\
        --lookups results/lookup/*_reads.txt \\
        --output threshold_recommendations.json

Author:  KM
Created: 2026-02
"""

import argparse
import json
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../prod_script/scripts'))
import utils

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


def load_data(assignment_files, lookup_files):
    """Load assignment results and ground truth."""
    # Build lookup: ZMW -> ground truth
    lookup_map = {}
    for lookup_file in lookup_files:
        df = pd.read_csv(lookup_file, sep='\t', dtype=str)
        df['Read'] = df['Read'].str.replace('/ccs', '', regex=False)
        df = df.rename(columns={'Read': 'ZMW'})
        lookup_map.update(dict(zip(df['ZMW'], df['KinnexLib'])))
    
    print(f"Loaded {len(lookup_map):,} ground truth ZMWs")
    
    # Load assignments — warn if files were produced with different parameters
    seen_params = {}
    all_data = []
    for assign_file in assignment_files:
        print(f"Processing {assign_file}...")
        header = utils.parse_assignment_header(assign_file)
        param_sig = {k: v for k, v in header.items() if k not in ('run_date', 'git')}
        if seen_params and param_sig != seen_params:
            print(f"WARNING: {Path(assign_file).name} has different scoring parameters "
                  f"than previous files — mixing may produce inconsistent results.", flush=True)
        seen_params = param_sig
        df = utils.load_assignments_df(assign_file)
        df['Top_Posterior'] = pd.to_numeric(df['Top_Posterior'], errors='coerce')
        
        for _, row in df.iterrows():
            zmw = row['ZMW']
            if zmw not in lookup_map:
                continue
            
            true_array = lookup_map[zmw]
            assigned_array = row['Assigned_Array']
            
            if assigned_array is None or pd.isna(assigned_array):
                continue
            
            all_data.append({
                'zmw': zmw,
                'posterior': row['Top_Posterior'],
                'classification': row['Classification'],
                'correct': (assigned_array == true_array),
                'assigned': assigned_array,
                'true': true_array
            })
    
    df = pd.DataFrame(all_data)
    print(f"\nTotal assignments: {len(df):,}")
    print(f"  Correct: {df['correct'].sum():,} ({df['correct'].mean()*100:.2f}%)")
    print(f"  Incorrect: {(~df['correct']).sum():,} ({(~df['correct']).mean()*100:.2f}%)")

    return df, seen_params


def analyze_posterior_distribution(df):
    """Analyze how posterior relates to correctness."""
    print("\n=== POSTERIOR DISTRIBUTION ANALYSIS ===")
    
    # Bin posteriors
    bins = [0, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99, 1.0]
    df['posterior_bin'] = pd.cut(df['posterior'], bins=bins, include_lowest=True)
    
    # Calculate accuracy per bin
    accuracy_by_bin = df.groupby('posterior_bin', observed=True).agg({
        'correct': ['sum', 'count', 'mean']
    }).round(4)
    
    print("\nAccuracy by Posterior Range:")
    print(accuracy_by_bin.to_string())
    
    return accuracy_by_bin


def find_optimal_thresholds(df, target_accuracy=0.999):
    """
    Find optimal thresholds for HIGH_CONF and LOW_CONF.
    
    HIGH_CONF threshold: Maximize data retention while achieving target accuracy
    LOW_CONF threshold: Accept lower accuracy, capture remaining useful data
    """
    print(f"\n=== THRESHOLD OPTIMIZATION ===")
    print(f"Target HIGH_CONF accuracy: {target_accuracy*100:.2f}%")
    
    # Test different thresholds
    thresholds = np.arange(0.50, 1.001, 0.01)
    results = []
    
    for thresh in thresholds:
        mask = df['posterior'] >= thresh
        filtered = df[mask]
        
        if len(filtered) == 0:
            continue
        
        accuracy = filtered['correct'].mean()
        n_kept = len(filtered)
        pct_kept = len(filtered) / len(df) * 100
        n_errors = (~filtered['correct']).sum()
        
        results.append({
            'threshold': thresh,
            'accuracy': accuracy,
            'n_kept': n_kept,
            'pct_kept': pct_kept,
            'n_errors': n_errors
        })
    
    results_df = pd.DataFrame(results)
    
    # Find HIGH_CONF threshold: highest threshold that still keeps good data
    # while achieving target accuracy
    high_conf_candidates = results_df[results_df['accuracy'] >= target_accuracy]
    
    if len(high_conf_candidates) == 0:
        print(f"WARNING: Cannot achieve {target_accuracy*100:.2f}% accuracy at any threshold!")
        print(f"Best accuracy: {results_df['accuracy'].max()*100:.2f}% at threshold {results_df.loc[results_df['accuracy'].idxmax(), 'threshold']:.3f}")
        high_conf_thresh = results_df.loc[results_df['accuracy'].idxmax(), 'threshold']
    else:
        # Choose the LOWEST threshold that achieves target (keeps most data)
        high_conf_thresh = high_conf_candidates['threshold'].min()
    
    high_conf_stats = results_df[results_df['threshold'] == high_conf_thresh].iloc[0]
    
    print(f"\nHIGH_CONF Recommendation: {high_conf_thresh:.3f}")
    print(f"  Accuracy: {high_conf_stats['accuracy']*100:.4f}%")
    print(f"  Keeps: {high_conf_stats['n_kept']:,} assignments ({high_conf_stats['pct_kept']:.1f}%)")
    print(f"  Errors: {int(high_conf_stats['n_errors'])}")
    
    # Find LOW_CONF threshold: much lower, to capture remaining useful data
    # Typically we want LOW_CONF to have >95% accuracy
    low_conf_target = 0.95
    low_conf_candidates = results_df[
        (results_df['threshold'] < high_conf_thresh) & 
        (results_df['accuracy'] >= low_conf_target)
    ]
    
    if len(low_conf_candidates) == 0:
        # Fallback: find threshold that gives 90% accuracy
        low_conf_candidates = results_df[
            (results_df['threshold'] < high_conf_thresh) &
            (results_df['accuracy'] >= 0.90)
        ]
    
    if len(low_conf_candidates) > 0:
        low_conf_thresh = low_conf_candidates['threshold'].min()
        low_conf_stats = results_df[results_df['threshold'] == low_conf_thresh].iloc[0]
        
        print(f"\nLOW_CONF Recommendation: {low_conf_thresh:.3f}")
        print(f"  Accuracy: {low_conf_stats['accuracy']*100:.4f}%")
        print(f"  Keeps: {low_conf_stats['n_kept']:,} assignments ({low_conf_stats['pct_kept']:.1f}%)")
        print(f"  Errors: {int(low_conf_stats['n_errors'])}")
    else:
        low_conf_thresh = 0.5
        print(f"\nLOW_CONF: Using default 0.500 (no optimal threshold found)")
    
    # Calculate what gets lost
    unassigned = len(df[df['posterior'] < low_conf_thresh])
    print(f"\nUNASSIGNED: {unassigned:,} ({unassigned/len(df)*100:.1f}%)")
    
    return {
        'HIGH_CONF': {
            'threshold': float(high_conf_thresh),
            'accuracy': float(high_conf_stats['accuracy']),
            'n_kept': int(high_conf_stats['n_kept']),
            'pct_kept': float(high_conf_stats['pct_kept']),
            'errors': int(high_conf_stats['n_errors'])
        },
        'LOW_CONF': {
            'threshold': float(low_conf_thresh),
            'accuracy': float(low_conf_stats['accuracy']) if len(low_conf_candidates) > 0 else 0.0,
            'n_kept': int(low_conf_stats['n_kept']) if len(low_conf_candidates) > 0 else 0,
            'pct_kept': float(low_conf_stats['pct_kept']) if len(low_conf_candidates) > 0 else 0.0,
            'errors': int(low_conf_stats['n_errors']) if len(low_conf_candidates) > 0 else 0
        }
    }, results_df


def plot_threshold_analysis(results_df, recommendations, output_prefix):
    """Create visualizations of threshold optimization."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Accuracy vs Threshold
    ax = axes[0, 0]
    ax.plot(results_df['threshold'], results_df['accuracy'] * 100, 'b-', linewidth=2)
    ax.axhline(y=99.9, color='r', linestyle='--', alpha=0.5, label='99.9% target')
    ax.axvline(x=recommendations['HIGH_CONF']['threshold'], color='g', linestyle='--', 
               alpha=0.7, label=f"HIGH_CONF: {recommendations['HIGH_CONF']['threshold']:.3f}")
    ax.axvline(x=recommendations['LOW_CONF']['threshold'], color='orange', linestyle='--',
               alpha=0.7, label=f"LOW_CONF: {recommendations['LOW_CONF']['threshold']:.3f}")
    ax.set_xlabel('Posterior Threshold')
    ax.set_ylabel('Accuracy (%)')
    ax.set_title('Accuracy vs Posterior Threshold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([95, 100.1])
    
    # Plot 2: Data Retention vs Threshold
    ax = axes[0, 1]
    ax.plot(results_df['threshold'], results_df['pct_kept'], 'b-', linewidth=2)
    ax.axvline(x=recommendations['HIGH_CONF']['threshold'], color='g', linestyle='--', alpha=0.7)
    ax.axvline(x=recommendations['LOW_CONF']['threshold'], color='orange', linestyle='--', alpha=0.7)
    ax.set_xlabel('Posterior Threshold')
    ax.set_ylabel('Data Retained (%)')
    ax.set_title('Data Retention vs Posterior Threshold')
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Error Count vs Threshold
    ax = axes[1, 0]
    ax.semilogy(results_df['threshold'], results_df['n_errors'], 'r-', linewidth=2)
    ax.axvline(x=recommendations['HIGH_CONF']['threshold'], color='g', linestyle='--', alpha=0.7)
    ax.axvline(x=recommendations['LOW_CONF']['threshold'], color='orange', linestyle='--', alpha=0.7)
    ax.set_xlabel('Posterior Threshold')
    ax.set_ylabel('Number of Errors (log scale)')
    ax.set_title('Errors vs Posterior Threshold')
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Accuracy vs Data Retention (Pareto curve)
    ax = axes[1, 1]
    ax.plot(results_df['pct_kept'], results_df['accuracy'] * 100, 'b-', linewidth=2)
    
    # Mark recommended thresholds
    hc = recommendations['HIGH_CONF']
    lc = recommendations['LOW_CONF']
    ax.plot(hc['pct_kept'], hc['accuracy']*100, 'go', markersize=10, label='HIGH_CONF')
    ax.plot(lc['pct_kept'], lc['accuracy']*100, 'o', color='orange', markersize=10, label='LOW_CONF')
    
    ax.set_xlabel('Data Retained (%)')
    ax.set_ylabel('Accuracy (%)')
    ax.set_title('Accuracy vs Data Retention Trade-off')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim([95, 100.1])
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_threshold_analysis.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved plot: {output_prefix}_threshold_analysis.png")
    plt.close()


def compare_current_vs_recommended(df, current_thresholds, recommended_thresholds):
    """Compare current thresholds vs recommended."""
    print("\n=== CURRENT vs RECOMMENDED COMPARISON ===")
    
    # Current thresholds
    current_hc = df[df['posterior'] >= current_thresholds['HIGH_CONF']]
    current_lc = df[(df['posterior'] >= current_thresholds['LOW_CONF']) & 
                    (df['posterior'] < current_thresholds['HIGH_CONF'])]
    
    # Recommended thresholds
    rec_hc = df[df['posterior'] >= recommended_thresholds['HIGH_CONF']['threshold']]
    rec_lc = df[(df['posterior'] >= recommended_thresholds['LOW_CONF']['threshold']) &
                (df['posterior'] < recommended_thresholds['HIGH_CONF']['threshold'])]
    
    print("\nCURRENT Thresholds (HIGH_CONF ≥ {:.2f}, LOW_CONF ≥ {:.2f}):".format(
        current_thresholds['HIGH_CONF'], current_thresholds['LOW_CONF']))
    print(f"  HIGH_CONF: {len(current_hc):,} assignments, {current_hc['correct'].mean()*100:.4f}% accuracy")
    print(f"  LOW_CONF:  {len(current_lc):,} assignments, {current_lc['correct'].mean()*100:.4f}% accuracy")
    
    print("\nRECOMMENDED Thresholds (HIGH_CONF ≥ {:.3f}, LOW_CONF ≥ {:.3f}):".format(
        recommended_thresholds['HIGH_CONF']['threshold'],
        recommended_thresholds['LOW_CONF']['threshold']))
    print(f"  HIGH_CONF: {len(rec_hc):,} assignments, {rec_hc['correct'].mean()*100:.4f}% accuracy")
    print(f"  LOW_CONF:  {len(rec_lc):,} assignments, {rec_lc['correct'].mean()*100:.4f}% accuracy")
    
    # Show difference
    print("\nIMPACT of changing thresholds:")
    print(f"  HIGH_CONF: {len(rec_hc) - len(current_hc):+,} assignments ({(len(rec_hc) - len(current_hc))/len(df)*100:+.1f}%)")
    print(f"  LOW_CONF:  {len(rec_lc) - len(current_lc):+,} assignments ({(len(rec_lc) - len(current_lc))/len(df)*100:+.1f}%)")


def main():
    parser = argparse.ArgumentParser(description="Optimize posterior thresholds after weight tuning")
    parser.add_argument('--assignments', nargs='+', required=True, help='Assignment result files')
    parser.add_argument('--lookups', nargs='+', required=True, help='Ground truth lookup files')
    parser.add_argument('--output', default='threshold_recommendations.json', help='Output file')
    parser.add_argument('--current-high-conf', type=float, default=None,
                        help='Current HIGH_CONF threshold for comparison (default: read from assignment file header)')
    parser.add_argument('--current-low-conf', type=float, default=None,
                        help='Current LOW_CONF threshold for comparison (default: read from assignment file header)')
    parser.add_argument('--target-accuracy', type=float, default=0.9999,
                        help='Target accuracy for HIGH_CONF (default: 0.9999 = 99.99%%)')
    
    args = parser.parse_args()
    
    # Expand wildcards
    assignment_files = []
    for pattern in args.assignments:
        assignment_files.extend(Path('.').glob(pattern))
    
    lookup_files = []
    for pattern in args.lookups:
        lookup_files.extend(Path('.').glob(pattern))
    
    if not assignment_files or not lookup_files:
        print("ERROR: No files found")
        sys.exit(1)
    
    print(f"Found {len(assignment_files)} assignment files")
    print(f"Found {len(lookup_files)} lookup files")
    
    # Load data
    df, seen_params = load_data(assignment_files, lookup_files)

    # Seed --current-* from assignment file header when not explicitly provided
    if args.current_high_conf is None:
        args.current_high_conf = float(seen_params.get('POSTERIOR_HIGH_CONF', 0.840))
    if args.current_low_conf is None:
        args.current_low_conf = float(seen_params.get('POSTERIOR_LOW_CONF', 0.50))

    # Analyze distribution
    analyze_posterior_distribution(df)
    
    # Find optimal thresholds
    recommendations, results_df = find_optimal_thresholds(df, target_accuracy=args.target_accuracy)
    
    # Plot
    output_prefix = args.output.replace('.json', '')
    plot_threshold_analysis(results_df, recommendations, output_prefix)
    
    # Compare with current
    current_thresholds = {
        'HIGH_CONF': args.current_high_conf,
        'LOW_CONF': args.current_low_conf
    }
    compare_current_vs_recommended(df, current_thresholds, recommendations)
    
    # Save recommendations
    output_data = {
        'recommended_thresholds': recommendations,
        'current_thresholds': current_thresholds,
        'target_accuracy': args.target_accuracy,
        'total_assignments': len(df),
        'overall_accuracy': float(df['correct'].mean())
    }
    
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved to {args.output}")
    print("\nTo use these thresholds, update the DEFAULT PARAMETERS section of assign_kinnex.py:")
    print(f"  POSTERIOR_HIGH_CONF = {recommendations['HIGH_CONF']['threshold']:.3f}")
    print(f"  POSTERIOR_LOW_CONF = {recommendations['LOW_CONF']['threshold']:.3f}")


if __name__ == '__main__':
    main()

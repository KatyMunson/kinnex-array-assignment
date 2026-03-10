#!/usr/bin/env python3
"""
visualize_posteriors.py

Visualises assignment posterior distributions to inform threshold selection.
Compares distributions for correct vs incorrect assignments, broken down by
confidence tier (HIGH_CONF / LOW_CONF). Requires ground truth lookup from
Snakefile_subsample to label correctness.

Called by Snakefile_qc (rule visualize_posteriors), or run standalone:
    python visualize_posteriors.py <assignments.txt> <lookup.txt> <output_prefix>

Author:  KM
Created: 2026-02
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../prod_script/scripts'))
import utils

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

def load_assignment_data(assignment_file, lookup_file):
    """Load and merge assignment results with ground truth."""
    # Load data
    truth_df = pd.read_csv(lookup_file, sep="\t", dtype=str)
    assign_df = utils.load_assignments_df(assignment_file)
    
    # Strip /ccs suffix from truth
    truth_df['Read'] = truth_df['Read'].str.replace('/ccs', '', regex=False)
    truth_df = truth_df.rename(columns={'Read': 'ZMW'})
    
    # Merge
    merged = pd.merge(assign_df, truth_df, on="ZMW", how="left")
    
    # Convert numeric columns
    merged['Top_Posterior'] = pd.to_numeric(merged['Top_Posterior'])
    merged['Specific_Barcodes'] = pd.to_numeric(merged['Specific_Barcodes'])
    merged['Shared_Barcodes']   = pd.to_numeric(merged['Shared_Barcodes'])
    
    # Calculate segments per ZMW from All_Barcodes column
    merged['Segments_per_ZMW'] = merged['All_Barcodes'].str.count(',') + 1
    
    # Add correctness flag
    merged['Correct'] = merged['Assigned_Array'] == merged['KinnexLib']
    
    # Add total barcodes and specific fraction
    merged['Total_Barcodes']    = merged['Specific_Barcodes'] + merged['Shared_Barcodes']
    merged['Specific_Fraction'] = merged['Specific_Barcodes'] / merged['Total_Barcodes']
    
    return merged

def plot_posterior_distributions(df, plot_path):
    """Create comprehensive posterior distribution plots."""
    
    # Set style
    sns.set_style("whitegrid")
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Overall posterior distribution by correctness
    ax1 = plt.subplot(2, 3, 1)
    
    # Filter out posterior = 0 for better visualization
    df_filtered = df[df['Top_Posterior'] > 0].copy()
    
    for correct in [True, False]:
        data = df_filtered[df_filtered['Correct'] == correct]['Top_Posterior']
        label = 'Correct' if correct else 'Incorrect'
        color = 'green' if correct else 'red'
        # Use bins that focus on 0.5-1.0 where most data is
        ax1.hist(data, bins=50, range=(0.5, 1.0), alpha=0.6, label=label, color=color, density=True)
    ax1.set_xlim(0.5, 1.0)
    ax1.set_xlabel('Top Posterior')
    ax1.set_ylabel('Density')
    ax1.set_title('Overall Posterior Distribution (excluding 0)')
    ax1.legend()
    ax1.axvline(x=0.89, color='blue', linestyle='--', alpha=0.5, linewidth=1)
    ax1.axvline(x=0.93, color='purple', linestyle='--', alpha=0.5, linewidth=1)
    ax1.axvline(x=0.6, color='green', linestyle='--', alpha=0.5, linewidth=1)
    ax1.text(0.89, ax1.get_ylim()[1]*0.95, 'std', fontsize=8, ha='center')
    ax1.text(0.93, ax1.get_ylim()[1]*0.85, 'relax', fontsize=8, ha='center')
    ax1.text(0.6, ax1.get_ylim()[1]*0.75, 'lowconf', fontsize=8, ha='center')
    
    # 2. HIGH_CONF only
    ax2 = plt.subplot(2, 3, 2)
    high_conf = df[df['Classification'] == 'HIGH_CONF']
    high_conf_filtered = high_conf[high_conf['Top_Posterior'] > 0].copy()
    
    for correct in [True, False]:
        data = high_conf_filtered[high_conf_filtered['Correct'] == correct]['Top_Posterior']
        label = f'Correct (n={len(data)})'  if correct else f'Incorrect (n={len(data)})'
        color = 'green' if correct else 'red'
        ax2.hist(data, range=(0.5, 1.0), bins=50, alpha=0.6, label=label, color=color, density=True)
    ax2.set_xlim(0.5, 1.0)
    ax2.set_xlabel('Top Posterior')
    ax2.set_ylabel('Density')
    ax2.set_title('HIGH_CONF Only (excluding 0)')
    ax2.legend()
    
    # 3. LOW_CONF only
    ax3 = plt.subplot(2, 3, 3)
    low_conf = df[df['Classification'] == 'LOW_CONF']
    low_conf_filtered = low_conf[low_conf['Top_Posterior'] > 0].copy()
    
    for correct in [True, False]:
        data = low_conf_filtered[low_conf_filtered['Correct'] == correct]['Top_Posterior']
        label = f'Correct (n={len(data)})' if correct else f'Incorrect (n={len(data)})'
        color = 'green' if correct else 'red'
        ax3.hist(data, bins=50, range=(0.5, 1.0), alpha=0.6, label=label, color=color, density=True)
    ax3.set_xlim(0.5, 1.0)
    ax3.set_xlabel('Top Posterior')
    ax3.set_ylabel('Density')
    ax3.set_title('LOW_CONF Only (excluding 0)')
    ax3.legend()
    
    # 4. Posterior vs Specific Barcodes (correct assignments)
    ax4 = plt.subplot(2, 3, 4)
    correct_df = df[df['Correct'] == True]
    scatter = ax4.scatter(correct_df['Specific_Barcodes'], correct_df['Top_Posterior'],
                         alpha=0.3, s=1, c='green')
    ax4.set_xlabel('Specific Barcodes')
    ax4.set_ylabel('Top Posterior')
    ax4.set_title('Correct Assignments')
    ax4.axhline(y=0.89, color='blue', linestyle='--', alpha=0.5)
    ax4.axhline(y=0.93, color='purple', linestyle='--', alpha=0.5)
    
    # 5. Posterior vs Specific Barcodes (incorrect assignments)
    ax5 = plt.subplot(2, 3, 5)
    incorrect_df = df[df['Correct'] == False]
    scatter = ax5.scatter(incorrect_df['Specific_Barcodes'], incorrect_df['Top_Posterior'],
                         alpha=0.5, s=2, c='red')
    ax5.set_xlabel('Specific Barcodes')
    ax5.set_ylabel('Top Posterior')
    ax5.set_title(f'Incorrect Assignments (n={len(incorrect_df)})')
    ax5.axhline(y=0.89, color='blue', linestyle='--', alpha=0.5)
    ax5.axhline(y=0.93, color='purple', linestyle='--', alpha=0.5)
    
    # 6. Segments per ZMW distribution
    ax6 = plt.subplot(2, 3, 6)
    for correct in [True, False]:
        data = df[df['Correct'] == correct]['Segments_per_ZMW']
        label = f'Correct (n={len(data)})' if correct else f'Incorrect (n={len(data)})'
        color = 'green' if correct else 'red'
        ax6.hist(data, bins=range(1, 10), alpha=0.6, label=label, color=color, density=True, align='left')
    ax6.set_xlabel('Segments per ZMW')
    ax6.set_ylabel('Density')
    ax6.set_title('Segments per ZMW Distribution')
    ax6.set_xticks(range(1, 9))
    ax6.legend()
    ax6.axvline(x=3, color='blue', linestyle='--', alpha=0.5, label='3 segments')
    
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Saved plot: {plot_path}")
    plt.close()

def generate_threshold_table(df):
    """Generate a table showing error rates at different thresholds for HIGH_CONF only."""
    
    # Focus on HIGH_CONF assignments only
    high_conf_df = df[df['Classification'] == 'HIGH_CONF'].copy()
    
    thresholds = [
        {'name': 'Current HIGH_CONF',       'posterior': 0.0,  'min_specific_bc': 0, 'specific_frac': 0.0, 'min_segments': 0},
        {'name': 'Posterior ≥ 0.90',        'posterior': 0.90, 'min_specific_bc': 0, 'specific_frac': 0.0, 'min_segments': 0},
        {'name': 'Posterior ≥ 0.95',        'posterior': 0.95, 'min_specific_bc': 0, 'specific_frac': 0.0, 'min_segments': 0},
        {'name': 'Posterior ≥ 0.99',        'posterior': 0.99, 'min_specific_bc': 0, 'specific_frac': 0.0, 'min_segments': 0},
        {'name': 'SpecBC ≥ 2',              'posterior': 0.0,  'min_specific_bc': 2, 'specific_frac': 0.0, 'min_segments': 0},
        {'name': 'SpecBC ≥ 3',              'posterior': 0.0,  'min_specific_bc': 3, 'specific_frac': 0.0, 'min_segments': 0},
        {'name': 'Segments ≥ 3',            'posterior': 0.0,  'min_specific_bc': 0, 'specific_frac': 0.0, 'min_segments': 3},
        {'name': 'Segments ≥ 4',            'posterior': 0.0,  'min_specific_bc': 0, 'specific_frac': 0.0, 'min_segments': 4},
        {'name': 'SpecFrac ≥ 0.5',          'posterior': 0.0,  'min_specific_bc': 0, 'specific_frac': 0.5, 'min_segments': 0},
        {'name': 'Post≥0.95 + SpecBC≥3',   'posterior': 0.95, 'min_specific_bc': 3, 'specific_frac': 0.0, 'min_segments': 0},
        {'name': 'Post≥0.95 + Seg≥3',      'posterior': 0.95, 'min_specific_bc': 0, 'specific_frac': 0.0, 'min_segments': 3},
        {'name': 'Post≥0.99 + Seg≥3',      'posterior': 0.99, 'min_specific_bc': 0, 'specific_frac': 0.0, 'min_segments': 3},
        {'name': 'Post≥0.95 + SpecFrac≥0.5','posterior': 0.95, 'min_specific_bc': 0, 'specific_frac': 0.5, 'min_segments': 0},
        {'name': 'Conservative',            'posterior': 0.95, 'min_specific_bc': 3, 'specific_frac': 0.5, 'min_segments': 0},
        {'name': 'Ultra-Conservative',      'posterior': 0.99, 'min_specific_bc': 3, 'specific_frac': 0.5, 'min_segments': 3},
    ]
    
    results = []
    
    for thresh in thresholds:
        # Apply filters to HIGH_CONF only
        mask = (high_conf_df['Top_Posterior']    >= thresh['posterior']) & \
               (high_conf_df['Specific_Barcodes'] >= thresh['min_specific_bc']) & \
               (high_conf_df['Specific_Fraction'] >= thresh['specific_frac']) & \
               (high_conf_df['Segments_per_ZMW']  >= thresh['min_segments'])
        
        filtered = high_conf_df[mask]
        
        if len(filtered) == 0:
            continue
        
        total = len(filtered)
        correct = filtered['Correct'].sum()
        incorrect = total - correct
        accuracy = correct / total * 100 if total > 0 else 0
        kept_pct_of_highconf = total / len(high_conf_df) * 100
        kept_pct_of_all = total / len(df) * 100
        
        results.append({
            'Threshold': thresh['name'],
            'HIGH_CONF_Kept': total,
            'Correct': correct,
            'Incorrect': incorrect,
            'Accuracy_%': f"{accuracy:.2f}",
            'Kept_%_of_HIGH_CONF': f"{kept_pct_of_highconf:.1f}",
            'Kept_%_of_All_ZMWs': f"{kept_pct_of_all:.1f}"
        })
    
    results_df = pd.DataFrame(results)
    return results_df

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python visualize_posteriors.py <assignment_file> <lookup_file> <plot_png> <threshold_csv>")
        print("Example: python visualize_posteriors.py assigned/Pool1_KN_bcM0001.txt lookup/Pool1_KN_bcM0001_reads.txt results/qc/per_pool/Pool1_KN_bcM0001/distributions.png results/qc/per_pool/Pool1_KN_bcM0001/posterior_threshold_analysis.csv")
        sys.exit(1)

    assignment_file = sys.argv[1]
    lookup_file     = sys.argv[2]
    plot_path       = sys.argv[3]
    threshold_path  = sys.argv[4]

    print(f"Loading data...")
    df = load_assignment_data(assignment_file, lookup_file)

    print(f"\nDataset summary:")
    print(f"  Total reads: {len(df):,}")
    print(f"  HIGH_CONF: {(df['Classification']=='HIGH_CONF').sum():,}")
    print(f"    - Correct: {df[(df['Classification']=='HIGH_CONF') & (df['Correct']==True)].shape[0]:,}")
    print(f"    - Incorrect: {df[(df['Classification']=='HIGH_CONF') & (df['Correct']==False)].shape[0]:,}")
    print(f"  LOW_CONF: {(df['Classification']=='LOW_CONF').sum():,}")
    print(f"    - Correct: {df[(df['Classification']=='LOW_CONF') & (df['Correct']==True)].shape[0]:,}")
    print(f"    - Incorrect: {df[(df['Classification']=='LOW_CONF') & (df['Correct']==False)].shape[0]:,}")
    print(f"  UNASSIGNED: {(df['Classification']=='UNASSIGNED').sum():,}")
    print(f"  Overall accuracy: {df['Correct'].sum()/len(df)*100:.2f}%")

    print(f"\nGenerating plots...")
    plot_posterior_distributions(df, plot_path)

    print(f"\nThreshold analysis (HIGH_CONF only):")
    threshold_table = generate_threshold_table(df)
    print(threshold_table.to_string(index=False))

    # Save threshold table to explicit path declared by Snakemake.
    # Named distinctly from check_assignments' threshold_analysis.csv, which
    # includes ground-truth correctness labels. This version uses posterior
    # distributions alone and is useful when ground truth is unavailable.
    threshold_table.to_csv(threshold_path, index=False)
    print(f"\nSaved threshold analysis: {threshold_path}")

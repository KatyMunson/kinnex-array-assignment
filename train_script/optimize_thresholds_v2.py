#!/usr/bin/env python3
"""
optimize_thresholds_v2.py

Multi-dimensional threshold optimisation for HIGH_CONF and LOW_CONF
classification. Extends optimize_thresholds.py with a grid search over:
  - posterior threshold
  - minimum total observations (segments) per ZMW
  - minimum informative barcode observations per ZMW (optional)

For each combination, computes accuracy and yield for HIGH_CONF, surfaces the
Pareto frontier, and recommends the operating point that maximises yield
subject to a target accuracy.

Usage:
    python optimize_thresholds_v2.py \\
        --assignments results/assigned/*.txt \\
        --lookups results/lookup/*_reads.txt \\
        --output threshold_recommendations_v2.json

    # Tighter accuracy target, skip informative sweep for speed:
    python optimize_thresholds_v2.py \\
        --assignments results/assigned/*.txt \\
        --lookups results/lookup/*_reads.txt \\
        --target-accuracy 0.9999 \\
        --no-informative-sweep \\
        --output threshold_recommendations_v2.json

Author:  KM
Created: 2026-02
"""

import argparse
import json
import sys
from pathlib import Path
from itertools import product

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')


# =======================
# DATA LOADING
# =======================

def load_data(assignment_files, lookup_files):
    """Load assignment results merged with ground truth."""
    lookup_map = {}
    for lookup_file in lookup_files:
        df = pd.read_csv(lookup_file, sep='\t', dtype=str)
        df['Read'] = df['Read'].str.replace('/ccs', '', regex=False)
        df = df.rename(columns={'Read': 'ZMW'})
        lookup_map.update(dict(zip(df['ZMW'], df['KinnexLib'])))

    print(f"Loaded {len(lookup_map):,} ground truth ZMWs from {len(lookup_files)} file(s)")

    all_rows = []
    for assign_file in assignment_files:
        print(f"  Loading {Path(assign_file).name}...")
        df = pd.read_csv(assign_file, sep='\t', dtype=str)

        for col in ('Top_Posterior', 'N_Observations', 'Informative_Barcodes',
                    'Uninformative_Barcodes', 'Extraneous_Barcodes'):
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        for _, row in df.iterrows():
            zmw = row['ZMW']
            if zmw not in lookup_map:
                continue
            assigned = row.get('Assigned_Array')
            if pd.isna(assigned):
                continue
            all_rows.append({
                'zmw':          zmw,
                'posterior':    row['Top_Posterior'],
                'n_obs':        row.get('N_Observations', np.nan),
                'n_inf':        row.get('Informative_Barcodes', np.nan),
                'n_uninf':      row.get('Uninformative_Barcodes', np.nan),
                'n_extr':       row.get('Extraneous_Barcodes', np.nan),
                'classification': row['Classification'],
                'correct':      (assigned == lookup_map[zmw]),
            })

    df = pd.DataFrame(all_rows)

    # Coerce gate columns to integer where possible (NaN → -1 means "unknown")
    for col in ('n_obs', 'n_inf'):
        df[col] = df[col].fillna(-1).astype(int)

    total = len(df)
    n_correct = df['correct'].sum()
    print(f"\nLoaded {total:,} assignments with ground truth")
    print(f"  Correct:   {n_correct:,} ({n_correct/total*100:.2f}%)")
    print(f"  Incorrect: {total - n_correct:,} ({(total - n_correct)/total*100:.2f}%)")
    print(f"  n_obs range:  {df['n_obs'].min()} – {df['n_obs'].max()}")
    print(f"  n_inf range:  {df['n_inf'].min()} – {df['n_inf'].max()}")

    return df


# =======================
# GRID SEARCH
# =======================

def grid_search(df, posterior_thresholds, min_obs_values, min_inf_values):
    """
    Sweep all combinations of (posterior, min_obs, min_inf).

    Returns a DataFrame with one row per combination containing:
        posterior, min_obs, min_inf, n_kept, pct_kept, n_correct,
        n_errors, accuracy
    """
    n_total = len(df)
    records = []

    total_combos = len(posterior_thresholds) * len(min_obs_values) * len(min_inf_values)
    print(f"\nRunning {total_combos:,} threshold combinations "
          f"({len(posterior_thresholds)} posterior × "
          f"{len(min_obs_values)} min_obs × "
          f"{len(min_inf_values)} min_inf)...")

    for post, min_obs, min_inf in product(posterior_thresholds, min_obs_values, min_inf_values):
        mask = (df['posterior'] >= post)
        if min_obs > 0:
            mask &= (df['n_obs'] >= min_obs)
        if min_inf > 0:
            mask &= (df['n_inf'] >= min_inf)

        kept = df[mask]
        n_kept = len(kept)
        if n_kept == 0:
            continue

        n_correct = kept['correct'].sum()
        n_errors  = n_kept - n_correct
        accuracy  = n_correct / n_kept

        records.append({
            'posterior': round(post, 4),
            'min_obs':   min_obs,
            'min_inf':   min_inf,
            'n_kept':    n_kept,
            'pct_kept':  n_kept / n_total * 100,
            'n_correct': n_correct,
            'n_errors':  n_errors,
            'accuracy':  accuracy,
        })

    result = pd.DataFrame(records)
    print(f"  {len(result):,} non-empty combinations")
    return result


def pareto_frontier(grid_df):
    """
    Return the Pareto-optimal rows: those not dominated on both
    accuracy (higher = better) and pct_kept (higher = better).
    """
    dominated = np.zeros(len(grid_df), dtype=bool)
    acc  = grid_df['accuracy'].values
    kept = grid_df['pct_kept'].values

    for i in range(len(grid_df)):
        # i is dominated if any j has acc[j] >= acc[i] AND kept[j] >= kept[i]
        # with at least one strict
        for j in range(len(grid_df)):
            if i == j:
                continue
            if acc[j] >= acc[i] and kept[j] >= kept[i]:
                if acc[j] > acc[i] or kept[j] > kept[i]:
                    dominated[i] = True
                    break

    return grid_df[~dominated].sort_values('pct_kept', ascending=False).reset_index(drop=True)


def recommend(grid_df, target_accuracy, label="HIGH_CONF"):
    """
    From the grid results, pick the operating point that maximises yield
    subject to achieving target_accuracy. Returns the best row as a dict.
    """
    candidates = grid_df[grid_df['accuracy'] >= target_accuracy]

    if len(candidates) == 0:
        best_acc = grid_df['accuracy'].max()
        print(f"\nWARNING [{label}]: Cannot achieve {target_accuracy*100:.2f}% accuracy "
              f"at any threshold combination. Best achievable: {best_acc*100:.4f}%")
        candidates = grid_df[grid_df['accuracy'] == best_acc]

    # Among candidates, take the one with highest yield; break ties by lowest posterior
    best = (candidates
            .sort_values(['pct_kept', 'posterior'], ascending=[False, True])
            .iloc[0])

    print(f"\n{label} Recommendation (target ≥ {target_accuracy*100:.2f}% accuracy):")
    print(f"  Posterior  ≥ {best['posterior']:.4f}")
    print(f"  Min obs    ≥ {int(best['min_obs'])}")
    print(f"  Min inf    ≥ {int(best['min_inf'])}")
    print(f"  Accuracy:    {best['accuracy']*100:.4f}%  ({int(best['n_errors'])} errors)")
    print(f"  Yield:       {best['pct_kept']:.2f}%  ({int(best['n_kept']):,} assignments)")

    return best.to_dict()


# =======================
# NUMERICAL TABLE OUTPUT
# =======================

def save_numerical_tables(grid_df, pareto_df, recommendations, output_prefix):
    """
    Save numerical summary tables as TSV files and print key ones to stdout.

    Outputs:
      <prefix>_grid_all.tsv       – full grid search results
      <prefix>_pareto.tsv         – Pareto-frontier rows only
      <prefix>_summary_by_post.tsv – accuracy/yield aggregated over obs/inf gates
                                     for a quick 1D view
    """
    # ---- Full grid ----
    grid_out = grid_df.copy()
    grid_out['accuracy_pct'] = (grid_out['accuracy'] * 100).round(4)
    grid_out['pct_kept']     = grid_out['pct_kept'].round(3)
    grid_path = f'{output_prefix}_grid_all.tsv'
    grid_out.to_csv(grid_path, sep='\t', index=False)
    print(f"  Saved full grid ({len(grid_out):,} rows): {grid_path}")

    # ---- Pareto frontier ----
    pareto_out = pareto_df.copy()
    pareto_out['accuracy_pct'] = (pareto_out['accuracy'] * 100).round(4)
    pareto_out['pct_kept']     = pareto_out['pct_kept'].round(3)
    pareto_path = f'{output_prefix}_pareto.tsv'
    pareto_out.to_csv(pareto_path, sep='\t', index=False)
    print(f"  Saved Pareto frontier ({len(pareto_out)} rows): {pareto_path}")

    # ---- Summary by posterior (best accuracy/yield at each posterior across all gates) ----
    summary = (
        grid_df
        .groupby('posterior')
        .agg(
            max_accuracy=('accuracy', 'max'),
            max_yield=('pct_kept', 'max'),
            min_errors=('n_errors', 'min'),
            n_combos=('accuracy', 'count'),
        )
        .reset_index()
    )
    summary['max_accuracy_pct'] = (summary['max_accuracy'] * 100).round(4)
    summary['max_yield']        = summary['max_yield'].round(3)
    summary_path = f'{output_prefix}_summary_by_post.tsv'
    summary.to_csv(summary_path, sep='\t', index=False)
    print(f"  Saved posterior summary ({len(summary)} rows): {summary_path}")

    # ---- Print the Pareto table to stdout (truncated) ----
    print("\n--- Pareto Frontier (top 30 by yield) ---")
    display_cols = ['posterior', 'min_obs', 'min_inf',
                    'accuracy_pct', 'pct_kept', 'n_errors', 'n_kept']
    print(pareto_out[display_cols].head(30).to_string(index=False))

    # ---- Print recommended operating points ----
    print("\n--- Recommended Operating Points ---")
    rec_rows = []
    for label, rec in recommendations.items():
        rec_rows.append({
            'label':        label,
            'posterior':    rec.get('posterior', np.nan),
            'min_obs':      int(rec.get('min_obs', 0)),
            'min_inf':      int(rec.get('min_inf', 0)),
            'accuracy_pct': round(rec.get('accuracy', np.nan) * 100, 4),
            'pct_kept':     round(rec.get('pct_kept', np.nan), 3),
            'n_errors':     int(rec.get('n_errors', 0)),
            'n_kept':       int(rec.get('n_kept', 0)),
        })
    rec_df = pd.DataFrame(rec_rows)
    print(rec_df.to_string(index=False))

    # Save recommendations table too
    rec_path = f'{output_prefix}_recommendations.tsv'
    rec_df.to_csv(rec_path, sep='\t', index=False)
    print(f"  Saved recommendations table: {rec_path}")

    return grid_path, pareto_path, summary_path, rec_path


# =======================
# PLOTTING
# =======================

def plot_1d_sweeps(df, grid_df, recommendations, output_prefix):
    """
    Classic 1D plots from v1: accuracy, yield, and error count vs posterior,
    computed at the recommended min_obs / min_inf gate from the 2D search.
    Also adds a posterior distribution panel.
    """
    rec = recommendations['HIGH_CONF']
    min_obs_fixed = int(rec['min_obs'])
    min_inf_fixed = int(rec['min_inf'])

    # Re-sweep posterior at the recommended gate values
    post_vals = np.arange(0.50, 1.001, 0.005)
    rows = []
    n_total = len(df)
    for p in post_vals:
        mask = (df['posterior'] >= p)
        if min_obs_fixed > 0:
            mask &= (df['n_obs'] >= min_obs_fixed)
        if min_inf_fixed > 0:
            mask &= (df['n_inf'] >= min_inf_fixed)
        kept = df[mask]
        if len(kept) == 0:
            continue
        rows.append({
            'posterior': p,
            'accuracy':  kept['correct'].mean(),
            'pct_kept':  len(kept) / n_total * 100,
            'n_errors':  int((~kept['correct']).sum()),
        })
    sweep = pd.DataFrame(rows)

    hc_post = rec['posterior']
    lc_rec  = recommendations.get('LOW_CONF', {})
    lc_post = lc_rec.get('posterior', None)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"1D Posterior Sweep  (min_obs ≥ {min_obs_fixed}, min_inf ≥ {min_inf_fixed})",
        fontsize=13, fontweight='bold'
    )

    # ---- accuracy ----
    ax = axes[0, 0]
    ax.plot(sweep['posterior'], sweep['accuracy'] * 100, 'steelblue', linewidth=2)
    ax.axhline(99.9, color='red', linestyle='--', alpha=0.5, label='99.9%')
    ax.axvline(hc_post, color='green', linestyle='--', alpha=0.8,
               label=f"HC: {hc_post:.3f}")
    if lc_post:
        ax.axvline(lc_post, color='orange', linestyle='--', alpha=0.8,
                   label=f"LC: {lc_post:.3f}")
    ax.set_xlabel('Posterior Threshold'); ax.set_ylabel('Accuracy (%)')
    ax.set_title('Accuracy vs Posterior')
    ax.set_ylim([95, 100.1]); ax.legend(); ax.grid(True, alpha=0.3)

    # ---- yield ----
    ax = axes[0, 1]
    ax.plot(sweep['posterior'], sweep['pct_kept'], 'steelblue', linewidth=2)
    ax.axvline(hc_post, color='green', linestyle='--', alpha=0.8)
    if lc_post:
        ax.axvline(lc_post, color='orange', linestyle='--', alpha=0.8)
    ax.set_xlabel('Posterior Threshold'); ax.set_ylabel('Yield (%)')
    ax.set_title('Yield vs Posterior')
    ax.grid(True, alpha=0.3)

    # ---- error count ----
    ax = axes[1, 0]
    ax.semilogy(sweep['posterior'], sweep['n_errors'] + 0.1, 'tomato', linewidth=2)
    ax.axvline(hc_post, color='green', linestyle='--', alpha=0.8)
    if lc_post:
        ax.axvline(lc_post, color='orange', linestyle='--', alpha=0.8)
    ax.set_xlabel('Posterior Threshold'); ax.set_ylabel('Errors (log scale)')
    ax.set_title('Error Count vs Posterior')
    ax.grid(True, alpha=0.3)

    # ---- posterior distributions ----
    ax = axes[1, 1]
    bins = np.linspace(0, 1, 60)
    correct   = df[df['correct']]['posterior'].dropna()
    incorrect = df[~df['correct']]['posterior'].dropna()
    ax.hist(correct,   bins=bins, alpha=0.55, color='steelblue', label='Correct',   density=True)
    ax.hist(incorrect, bins=bins, alpha=0.55, color='tomato',    label='Incorrect', density=True)
    ax.axvline(hc_post, color='green',  linestyle='--', alpha=0.8, label=f"HC: {hc_post:.3f}")
    if lc_post:
        ax.axvline(lc_post, color='orange', linestyle='--', alpha=0.8, label=f"LC: {lc_post:.3f}")
    ax.set_xlabel('Posterior'); ax.set_ylabel('Density')
    ax.set_title('Posterior Distribution by Correctness')
    ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = f'{output_prefix}_1d_sweep.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_heatmaps(grid_df, recommendations, output_prefix, row_dim='min_obs'):
    """
    For each unique value of the facet dimension (the *other* dim), draw a
    heatmap of accuracy and yield across (posterior [x], row_dim [y]) space.

    Parameters
    ----------
    row_dim : str
        The gate dimension to place on the **y-axis** of each heatmap panel.
        'min_obs'  → y-axis = min_obs,  facets iterate over min_inf values
        'min_inf'  → y-axis = min_inf,  facets iterate over min_obs values

    The output file is named:
        <prefix>_heatmap_rows<row_dim>.png
    so the two files are always distinguishable regardless of call order.
    """
    # The dimension that defines the facet panels (the one NOT on the y-axis)
    facet_dim  = 'min_inf' if row_dim == 'min_obs' else 'min_obs'
    facet_vals = sorted(grid_df[facet_dim].unique())
    n_facets   = len(facet_vals)

    fig, axes = plt.subplots(
        n_facets, 2,
        figsize=(14, 4 * n_facets),
        squeeze=False
    )
    fig.suptitle(
        f"Accuracy & Yield heatmaps  "
        f"(y-axis = {row_dim},  x-axis = posterior,  facets = {facet_dim})",
        fontsize=12, fontweight='bold'
    )

    rec = recommendations['HIGH_CONF']

    for fi, fval in enumerate(facet_vals):
        sub = grid_df[grid_df[facet_dim] == fval]

        # Pivot to (row_dim × posterior) matrices
        acc_piv  = sub.pivot_table(index=row_dim, columns='posterior',
                                   values='accuracy', aggfunc='mean')
        kept_piv = sub.pivot_table(index=row_dim, columns='posterior',
                                   values='pct_kept', aggfunc='mean')

        for ci, (piv, title, cmap, vmin, vmax) in enumerate([
            (acc_piv,  f'Accuracy  ({facet_dim}={fval})',  'RdYlGn', 0.95, 1.00),
            (kept_piv, f'Yield (%) ({facet_dim}={fval})',  'Blues',  0,    100),
        ]):
            ax = axes[fi][ci]
            im = ax.imshow(
                piv.values,
                aspect='auto', origin='lower',
                cmap=cmap, vmin=vmin, vmax=vmax,
                interpolation='nearest'
            )
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

            # Axis labels — show every Nth tick to avoid crowding
            col_labels = [f"{v:.2f}" for v in piv.columns]
            row_labels  = [str(int(v)) for v in piv.index]
            n_col = len(col_labels)
            step  = max(1, n_col // 10)
            ax.set_xticks(range(0, n_col, step))
            ax.set_xticklabels(col_labels[::step], rotation=45, ha='right', fontsize=7)
            ax.set_yticks(range(len(row_labels)))
            ax.set_yticklabels(row_labels, fontsize=8)
            ax.set_xlabel('Posterior Threshold', fontsize=9)
            ax.set_ylabel(row_dim, fontsize=9)          # ← always the y-axis variable
            ax.set_title(title, fontsize=10)

            # Mark recommended point if it falls in this facet
            if int(rec[facet_dim]) == int(fval):
                try:
                    col_idx = list(piv.columns).index(
                        min(piv.columns, key=lambda x: abs(x - rec['posterior']))
                    )
                    row_idx = list(piv.index).index(
                        min(piv.index, key=lambda x: abs(x - rec[row_dim]))
                    )
                    ax.plot(col_idx, row_idx, 'w*', markersize=14, markeredgecolor='black',
                            markeredgewidth=0.8, label='Recommended')
                    ax.legend(fontsize=8, loc='lower left')
                except (ValueError, KeyError):
                    pass

    plt.tight_layout()
    # File name encodes which variable is on the y-axis → no ambiguity
    path = f'{output_prefix}_heatmap_rows{row_dim}.png'
    plt.savefig(path, dpi=180, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_pareto(pareto_df, all_grid_df, recommendations, output_prefix):
    """
    Scatter of all grid points (accuracy vs yield), highlighting the
    Pareto frontier and the recommended operating point.
    Points are coloured by posterior threshold.
    """
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle("Pareto Frontier: Accuracy vs Yield", fontsize=13, fontweight='bold')

    for ax, colour_by, label in [
        (axes[0], 'posterior', 'Posterior threshold'),
        (axes[1], 'min_obs',   'Min observations'),
    ]:
        sc = ax.scatter(
            all_grid_df['pct_kept'],
            all_grid_df['accuracy'] * 100,
            c=all_grid_df[colour_by],
            cmap='viridis', s=6, alpha=0.35, linewidths=0
        )
        plt.colorbar(sc, ax=ax, label=label)

        # Pareto frontier as a line
        pf = pareto_df.sort_values('pct_kept')
        ax.plot(pf['pct_kept'], pf['accuracy'] * 100,
                'r-', linewidth=1.5, alpha=0.8, label='Pareto frontier', zorder=3)

        # Recommended point
        rec = recommendations['HIGH_CONF']
        ax.scatter(
            rec['pct_kept'], rec['accuracy'] * 100,
            color='red', s=120, marker='*', zorder=5,
            label=f"Recommended\n(p≥{rec['posterior']:.3f}, "
                  f"obs≥{int(rec['min_obs'])}, inf≥{int(rec['min_inf'])})"
        )

        ax.set_xlabel('Yield (%)', fontsize=11)
        ax.set_ylabel('Accuracy (%)', fontsize=11)
        ax.set_ylim([max(95, all_grid_df['accuracy'].min() * 100 - 0.5), 100.1])
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = f'{output_prefix}_pareto.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_obs_marginals(df, output_prefix):
    """
    Show how accuracy and ZMW counts vary with n_obs and n_inf independently,
    to give intuition for where the hard gates should sit.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Marginal Analysis: Observations and Informative Barcodes",
                 fontsize=13, fontweight='bold')

    for row_i, (col, xlabel) in enumerate([
        ('n_obs', 'Total Segments (N_Observations)'),
        ('n_inf', 'Informative Barcodes'),
    ]):
        vals = sorted(df[col].unique())
        vals = [v for v in vals if v >= 0]   # drop -1 (unknown)

        # Accuracy by value
        grp = df[df[col] >= 0].groupby(col)['correct'].agg(['mean', 'sum', 'count']).reset_index()
        grp.columns = [col, 'accuracy', 'n_correct', 'n_total']

        # ---- accuracy curve ----
        ax = axes[row_i][0]
        ax.bar(grp[col], grp['accuracy'] * 100, color='steelblue', alpha=0.7, width=0.8)
        ax.set_xlabel(xlabel); ax.set_ylabel('Accuracy (%)')
        ax.set_title(f'Accuracy by {col}')
        ax.set_ylim([max(0, grp['accuracy'].min() * 100 - 5), 100.5])
        ax.grid(True, alpha=0.3, axis='y')
        # annotate count
        for _, r in grp.iterrows():
            ax.text(r[col], r['accuracy'] * 100 + 0.2, f"n={int(r['n_total'])}", 
                    ha='center', va='bottom', fontsize=6, rotation=90)

        # ---- cumulative yield if gated at >= x ----
        n_total = len(df)
        cum_pct = []
        for v in grp[col]:
            kept = (df[col] >= v).sum()
            cum_pct.append(kept / n_total * 100)

        ax2 = axes[row_i][1]
        ax2b = ax2.twinx()
        ax2.bar(grp[col], grp['accuracy'] * 100, color='steelblue', alpha=0.5,
                width=0.8, label='Accuracy')
        ax2b.plot(grp[col], cum_pct, 'tomato', linewidth=2, marker='o',
                  markersize=4, label='Cumulative yield if gated ≥ x')
        ax2.set_xlabel(xlabel)
        ax2.set_ylabel('Accuracy (%)', color='steelblue')
        ax2.set_ylim([max(0, grp['accuracy'].min() * 100 - 5), 100.5])
        ax2b.set_ylabel('Yield if gated ≥ x (%)', color='tomato')
        ax2b.set_ylim([0, 105])
        ax2.set_title(f'Accuracy & Yield Trade-off at {col} Gate')
        ax2.grid(True, alpha=0.3, axis='y')

        # combined legend
        lines1, labels1 = ax2.get_legend_handles_labels()
        lines2, labels2 = ax2b.get_legend_handles_labels()
        ax2.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc='lower left')

    plt.tight_layout()
    path = f'{output_prefix}_obs_marginals.png'
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# =======================
# COMPARISON (v1 parity)
# =======================

def compare_current_vs_recommended(df, current, recommended):
    """Print a side-by-side comparison of current vs recommended thresholds."""
    print("\n=== CURRENT vs RECOMMENDED ===")

    def _stats(mask, label):
        kept = df[mask]
        if len(kept) == 0:
            print(f"  {label}: 0 assignments")
            return
        print(f"  {label}: {len(kept):,} assignments, "
              f"{kept['correct'].mean()*100:.4f}% accuracy, "
              f"{int((~kept['correct']).sum())} errors")

    rec = recommended['HIGH_CONF']
    lc_rec = recommended.get('LOW_CONF', {})

    current_hc_mask = (
        (df['posterior'] >= current['HIGH_CONF']) &
        (df['n_obs']     >= current.get('min_obs', 0)) &
        (df['n_inf']     >= current.get('min_inf', 0))
    )
    rec_hc_mask = (
        (df['posterior'] >= rec['posterior']) &
        (df['n_obs']     >= int(rec['min_obs'])) &
        (df['n_inf']     >= int(rec['min_inf']))
    )

    print(f"\nCurrent HIGH_CONF "
          f"(post≥{current['HIGH_CONF']}, "
          f"min_obs≥{current.get('min_obs',0)}, "
          f"min_inf≥{current.get('min_inf',0)}):")
    _stats(current_hc_mask, "HIGH_CONF")

    print(f"\nRecommended HIGH_CONF "
          f"(post≥{rec['posterior']:.4f}, "
          f"min_obs≥{int(rec['min_obs'])}, "
          f"min_inf≥{int(rec['min_inf'])}):")
    _stats(rec_hc_mask, "HIGH_CONF")

    delta = int(rec['n_kept']) - current_hc_mask.sum()
    pct_delta = delta / len(df) * 100
    print(f"\n  Δ assignments: {delta:+,} ({pct_delta:+.2f}% of total)")

    if lc_rec:
        current_lc_mask = (
            (df['posterior'] >= current['LOW_CONF']) &
            (df['posterior'] <  current['HIGH_CONF'])
        )
        rec_lc_mask = (
            (df['posterior'] >= lc_rec.get('posterior', 0)) &
            (df['posterior'] <  rec['posterior'])
        )
        print(f"\nCurrent LOW_CONF (post≥{current['LOW_CONF']}):")
        _stats(current_lc_mask, "LOW_CONF")
        print(f"\nRecommended LOW_CONF (post≥{lc_rec.get('posterior',0):.4f}):")
        _stats(rec_lc_mask, "LOW_CONF")


# =======================
# MAIN
# =======================

def main():
    parser = argparse.ArgumentParser(
        description="Multi-dimensional threshold optimisation for Kinnex array assignment",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--assignments', nargs='+', required=True,
                        help='Assignment result .txt files (glob patterns OK)')
    parser.add_argument('--lookups', nargs='+', required=True,
                        help='Ground truth lookup files')
    parser.add_argument('--output', default='threshold_recommendations_v2.json',
                        help='JSON output path')

    # Sweep ranges
    parser.add_argument('--posterior-min',  type=float, default=0.50)
    parser.add_argument('--posterior-max',  type=float, default=1.00)
    parser.add_argument('--posterior-step', type=float, default=0.01,
                        help='Step size for posterior sweep (smaller = finer, slower)')
    parser.add_argument('--max-obs',  type=int, default=6,
                        help='Maximum min_obs gate to test (tests 0 … max-obs)')
    parser.add_argument('--max-inf',  type=int, default=4,
                        help='Maximum min_inf gate to test (tests 0 … max-inf)')
    parser.add_argument('--no-informative-sweep', action='store_true',
                        help='Skip the min_inf dimension (only sweep posterior × min_obs)')

    # Targets
    parser.add_argument('--target-accuracy', type=float, default=0.9990,
                        help='Target accuracy for HIGH_CONF recommendation (0–1)')
    parser.add_argument('--low-conf-target-accuracy', type=float, default=0.95,
                        help='Target accuracy for LOW_CONF recommendation')

    # Current thresholds (for comparison)
    parser.add_argument('--current-high-conf',  type=float, default=0.840)
    parser.add_argument('--current-low-conf',   type=float, default=0.500)
    parser.add_argument('--current-min-obs',    type=int,   default=3)
    parser.add_argument('--current-min-inf',    type=int,   default=0)

    args = parser.parse_args()

    # Resolve file globs
    assignment_files = []
    for pat in args.assignments:
        found = list(Path('.').glob(pat))
        if not found and Path(pat).exists():
            found = [Path(pat)]
        assignment_files.extend(found)

    lookup_files = []
    for pat in args.lookups:
        found = list(Path('.').glob(pat))
        if not found and Path(pat).exists():
            found = [Path(pat)]
        lookup_files.extend(found)

    if not assignment_files:
        print("ERROR: No assignment files found"); sys.exit(1)
    if not lookup_files:
        print("ERROR: No lookup files found"); sys.exit(1)

    print(f"Assignment files: {len(assignment_files)}")
    print(f"Lookup files:     {len(lookup_files)}")

    # ---- Load ----
    df = load_data(assignment_files, lookup_files)

    # ---- Build sweep axes ----
    posterior_thresholds = np.arange(
        args.posterior_min,
        args.posterior_max + args.posterior_step / 2,
        args.posterior_step
    )
    min_obs_values = list(range(0, args.max_obs + 1))
    min_inf_values = [0] if args.no_informative_sweep else list(range(0, args.max_inf + 1))

    # ---- Grid search ----
    grid_df = grid_search(df, posterior_thresholds, min_obs_values, min_inf_values)

    if grid_df.empty:
        print("ERROR: Grid search returned no results"); sys.exit(1)

    # ---- Pareto frontier ----
    pareto_df = pareto_frontier(grid_df)
    print(f"\nPareto frontier: {len(pareto_df)} non-dominated combinations")
    print(pareto_df[['posterior', 'min_obs', 'min_inf',
                      'accuracy', 'pct_kept', 'n_errors']].head(20).to_string(index=False))

    # ---- Recommendations ----
    recommendations = {}
    recommendations['HIGH_CONF'] = recommend(grid_df, args.target_accuracy, "HIGH_CONF")

    # LOW_CONF: sweep posterior only (min_obs=0, min_inf=0) below the HC posterior
    lc_grid = grid_df[
        (grid_df['posterior'] < recommendations['HIGH_CONF']['posterior']) &
        (grid_df['min_obs'] == 0) &
        (grid_df['min_inf'] == 0)
    ]
    if len(lc_grid) > 0:
        recommendations['LOW_CONF'] = recommend(lc_grid, args.low_conf_target_accuracy, "LOW_CONF")
    else:
        recommendations['LOW_CONF'] = {'posterior': 0.5, 'min_obs': 0, 'min_inf': 0}

    # ---- Compare ----
    current = {
        'HIGH_CONF': args.current_high_conf,
        'LOW_CONF':  args.current_low_conf,
        'min_obs':   args.current_min_obs,
        'min_inf':   args.current_min_inf,
    }
    compare_current_vs_recommended(df, current, recommendations)

    # ---- Numerical tables ----
    output_prefix = str(args.output).replace('.json', '')
    print("\nGenerating numerical tables...")
    save_numerical_tables(grid_df, pareto_df, recommendations, output_prefix)

    # ---- Plots ----
    print("\nGenerating plots...")
    plot_obs_marginals(df, output_prefix)
    plot_1d_sweeps(df, grid_df, recommendations, output_prefix)
    plot_pareto(pareto_df, grid_df, recommendations, output_prefix)

    # Heatmap: y-axis = min_obs, facets = min_inf values → file: _heatmap_rowsmin_obs.png
    plot_heatmaps(grid_df, recommendations, output_prefix, row_dim='min_obs')

    # Heatmap: y-axis = min_inf, facets = min_obs values → file: _heatmap_rowsmin_inf.png
    if not args.no_informative_sweep and args.max_inf > 0:
        plot_heatmaps(grid_df, recommendations, output_prefix, row_dim='min_inf')

    # ---- Save JSON ----
    output_data = {
        'recommended_thresholds': recommendations,
        'current_thresholds': current,
        'target_accuracy':         args.target_accuracy,
        'lc_target_accuracy':      args.low_conf_target_accuracy,
        'total_assignments':       len(df),
        'overall_accuracy':        float(df['correct'].mean()),
        'pareto_frontier_size':    len(pareto_df),
        'sweep': {
            'posterior_range':   [float(posterior_thresholds[0]), float(posterior_thresholds[-1])],
            'posterior_step':    args.posterior_step,
            'min_obs_range':     [0, args.max_obs],
            'min_inf_range':     [0, 0 if args.no_informative_sweep else args.max_inf],
            'total_combinations': len(grid_df),
        }
    }

    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)

    print(f"\nResults saved to {args.output}")
    print("\nTo use these thresholds, update assign_kinnex.py:")
    rec = recommendations['HIGH_CONF']
    print(f"  POSTERIOR_HIGH_CONF = {rec['posterior']:.4f}")
    print(f"  MIN_OBS_HIGH_CONF   = {int(rec['min_obs'])}")
    if not args.no_informative_sweep:
        print(f"  MIN_INF_HIGH_CONF   = {int(rec['min_inf'])}  "
              f"(add this constant if adopting the inf gate)")
    lc = recommendations['LOW_CONF']
    print(f"  POSTERIOR_LOW_CONF  = {lc.get('posterior', 0.5):.4f}")


if __name__ == '__main__':
    main()

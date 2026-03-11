#!/usr/bin/env python3
"""
optimize_barcode_weights.py

Uses machine learning (logistic regression or random forest) to find optimal
scoring weights for the Kinnex array assignment classifier. Trains on pools
with known ground truth (from Snakefile_subsample) to minimise misassignment.

Run this before optimize_thresholds.py. Update the weight constants in
assign_kinnex.py with the recommended values from the JSON output.

Usage:
    python optimize_barcode_weights.py \\
        --arrays results/arrays/Pool1_KN_bcM0003.txt \\
        --assignments path/to/production/results/assigned/*.txt \\
        --lookups results/lookup/*_reads.txt \\
        --output optimized_weights.json

Requirements: scikit-learn, pandas, numpy

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
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')


def load_training_data(assignment_files, lookup_files, arrays):
    """
    Load assignment results and ground truth.
    
    Returns:
        training_data: list of dicts with features and labels
    """
    training_data = []
    
    # Build lookup: file -> ground truth mapping
    lookup_map = {}
    for lookup_file in lookup_files:
        df = pd.read_csv(lookup_file, sep='\t', dtype=str)
        # Strip /ccs suffix if present
        df['Read'] = df['Read'].str.replace('/ccs', '', regex=False)
        df = df.rename(columns={'Read': 'ZMW'})
        lookup_map.update(dict(zip(df['ZMW'], df['KinnexLib'])))
    
    print(f"Loaded {len(lookup_map):,} ground truth ZMWs")
    
    # Process assignment files — warn if files were produced with different parameters
    seen_params = {}
    for assign_file in assignment_files:
        pool_name = Path(assign_file).stem   # e.g. TRAIN_3way_KN_bcM0010
        print(f"Processing {assign_file}...")
        header = utils.parse_assignment_header(assign_file)
        param_sig = {k: v for k, v in header.items() if k not in ('run_date', 'git')}
        if seen_params and param_sig != seen_params:
            print(f"WARNING: {Path(assign_file).name} has different scoring parameters "
                  f"than previous files — mixing may produce inconsistent results.", flush=True)
        seen_params = param_sig
        df = utils.load_assignments_df(assign_file)
        
        # Convert numeric columns
        df['Top_Posterior'] = pd.to_numeric(df['Top_Posterior'], errors='coerce')
        df['Specific_Barcodes'] = pd.to_numeric(df['Specific_Barcodes'], errors='coerce')
        df['Shared_Barcodes'] = pd.to_numeric(df['Shared_Barcodes'], errors='coerce')
        df['Discordant_Barcodes'] = pd.to_numeric(df['Discordant_Barcodes'], errors='coerce')
        
        for _, row in df.iterrows():
            zmw = row['ZMW']
            
            # Skip if no ground truth
            if zmw not in lookup_map:
                continue
            
            true_array = lookup_map[zmw]
            assigned_array = row['Assigned_Array']
            classification = row['Classification']
            
            # Skip UNASSIGNED for training (we want to learn from assignments)
            if classification == 'UNASSIGNED' or assigned_array is None:
                continue
            
            # Parse barcodes
            barcodes = row['All_Barcodes'].split(',') if pd.notna(row['All_Barcodes']) else []
            
            # Extract features for this ZMW
            features = extract_features(
                zmw, barcodes, assigned_array, true_array, 
                arrays, classification, row
            )
            
            # Label: 1 if correct, 0 if incorrect
            label = 1 if assigned_array == true_array else 0
            
            training_data.append({
                'features': features,
                'label': label,
                'zmw': zmw,
                'pool': pool_name,
                'classification': classification,
                'assigned': assigned_array,
                'true': true_array,
                'n_specific':   float(row['Specific_Barcodes']),
                'n_shared':     float(row['Shared_Barcodes']),
                'n_discordant': float(row['Discordant_Barcodes']),
                'posterior': float(row['Top_Posterior']),
            })
    
    print(f"Collected {len(training_data):,} training examples")
    return training_data, seen_params


def extract_features(zmw, barcodes, assigned_array, true_array, arrays, classification, row):
    """
    Extract features for ML model.
    
    Features:
    - Raw counts: specific, shared, discordant barcodes
    - Ratios: specific/total, shared/total, discordant/total
    - Segment count
    - Barcode diversity (unique barcodes / total)
    - Evidence strength metrics
    """
    n_segments = len(barcodes)
    n_specific   = row['Specific_Barcodes']
    n_shared     = row['Shared_Barcodes']
    n_discordant = row['Discordant_Barcodes']
    n_total_bc = n_specific + n_shared + n_discordant
    
    # Barcode diversity
    bc_counter = Counter(barcodes)
    n_unique_bc = len(bc_counter)
    max_bc_count = max(bc_counter.values()) if bc_counter else 0
    
    # Calculate what barcodes SHOULD be for true array
    true_array_barcodes = arrays[true_array]['barcodes'] if true_array in arrays else set()
    
    # Count matches to TRUE array
    true_matches = sum(1 for bc in barcodes if bc in true_array_barcodes)
    true_unique_bc = len(set(barcodes) & true_array_barcodes)
    
    # Build feature vector
    features = {
        # Raw counts
        'n_segments':    n_segments,
        'n_specific':    n_specific,
        'n_shared':      n_shared,
        'n_discordant':  n_discordant,

        # Ratios (avoid division by zero)
        'specific_frac':    n_specific   / max(n_total_bc, 1),
        'shared_frac':      n_shared     / max(n_total_bc, 1),
        'discordant_frac':  n_discordant / max(n_total_bc, 1),
        
        # Barcode diversity
        'unique_bc_count': n_unique_bc,
        'bc_diversity': n_unique_bc / max(n_segments, 1),
        'max_bc_frequency': max_bc_count / max(n_segments, 1),
        
        # Evidence strength
        'true_matches': true_matches,
        'true_match_frac': true_matches / max(n_segments, 1),
        'true_unique_bc': true_unique_bc,
        
        # Interaction terms
        'specific_times_segments':    n_specific   * n_segments,
        'discordant_times_segments':  n_discordant * n_segments,
        'specific_minus_discordant':  n_specific   - n_discordant,
        
        # Posterior from current algorithm
        'current_posterior': row['Top_Posterior'],
    }
    
    return features


def train_model(training_data, model_type='logistic'):
    """
    Train ML model to predict correct vs incorrect assignments.
    
    Args:
        training_data: list of dicts with 'features' and 'label'
        model_type: 'logistic' or 'random_forest'
    
    Returns:
        trained model, feature names, performance metrics
    """
    # Convert to arrays
    feature_names = sorted(training_data[0]['features'].keys())
    X = np.array([[d['features'][f] for f in feature_names] for d in training_data])
    y = np.array([d['label'] for d in training_data])
    
    print(f"\nTraining on {len(X):,} examples")
    print(f"  Correct assignments: {sum(y):,} ({sum(y)/len(y)*100:.2f}%)")
    print(f"  Incorrect assignments: {len(y)-sum(y):,} ({(len(y)-sum(y))/len(y)*100:.2f}%)")
    
    # Train model
    if model_type == 'logistic':
        model = LogisticRegression(
            max_iter=1000,
            class_weight='balanced',  # Handle class imbalance
            C=1.0,  # Regularization strength
            random_state=42
        )
    else:  # random_forest
        model = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            class_weight='balanced',
            random_state=42,
            n_jobs=-1
        )
    
    model.fit(X, y)
    
    # Evaluate with cross-validation
    cv_scores = cross_val_score(model, X, y, cv=5, scoring='accuracy')
    
    print(f"\nModel Performance:")
    print(f"  Training accuracy: {model.score(X, y)*100:.2f}%")
    print(f"  Cross-val accuracy: {cv_scores.mean()*100:.2f}% (+/- {cv_scores.std()*100:.2f}%)")
    
    # Feature importance
    if model_type == 'logistic':
        importance = np.abs(model.coef_[0])
    else:
        importance = model.feature_importances_
    
    print(f"\nTop 10 Most Important Features:")
    feature_importance = sorted(zip(feature_names, importance), key=lambda x: x[1], reverse=True)
    for feat, imp in feature_importance[:10]:
        print(f"  {feat:25s}: {imp:.4f}")
    
    return model, feature_names, {
        'train_accuracy': model.score(X, y),
        'cv_mean': cv_scores.mean(),
        'cv_std': cv_scores.std()
    }


def analyze_errors(model, feature_names, training_data, arrays):
    """
    Analyze where the model is making errors.
    """
    X = np.array([[d['features'][f] for f in feature_names] for d in training_data])
    y = np.array([d['label'] for d in training_data])
    
    predictions = model.predict(X)
    probabilities = model.predict_proba(X)[:, 1]  # Probability of correct assignment
    
    print(f"\n=== ERROR ANALYSIS ===")
    
    # False Positives: Model predicts correct, but it's actually wrong
    fp_indices = np.where((predictions == 1) & (y == 0))[0]
    print(f"\nFalse Positives (predicted correct, actually wrong): {len(fp_indices)}")
    if len(fp_indices) > 0:
        print("Examples:")
        for idx in fp_indices[:5]:
            d = training_data[idx]
            print(f"  ZMW: {d['zmw'][:50]}")
            print(f"    True: {d['true']}, Assigned: {d['assigned']}, Class: {d['classification']}")
            print(f"    Model prob: {probabilities[idx]:.3f}")
            print(f"    Features: specific={d['features']['n_specific']}, shared={d['features']['n_shared']}, discordant={d['features']['n_discordant']}")
    
    # False Negatives: Model predicts wrong, but it's actually correct
    fn_indices = np.where((predictions == 0) & (y == 1))[0]
    print(f"\nFalse Negatives (predicted wrong, actually correct): {len(fn_indices)}")
    if len(fn_indices) > 0:
        print("Examples:")
        for idx in fn_indices[:5]:
            d = training_data[idx]
            print(f"  ZMW: {d['zmw'][:50]}")
            print(f"    True: {d['true']}, Assigned: {d['assigned']}, Class: {d['classification']}")
            print(f"    Model prob: {probabilities[idx]:.3f}")
    
    # Analyze by classification
    print(f"\n=== Performance by Classification ===")
    for cls in ['HIGH_CONF', 'LOW_CONF']:
        cls_indices = [i for i, d in enumerate(training_data) if d['classification'] == cls]
        if cls_indices:
            cls_acc = (predictions[cls_indices] == y[cls_indices]).mean()
            cls_total = len(cls_indices)
            cls_correct = sum(y[cls_indices])
            print(f"{cls:12s}: {cls_acc*100:.2f}% accuracy ({cls_correct}/{cls_total} actually correct)")


def extract_weight_recommendations(model, feature_names, model_type='logistic', seen_params=None):
    """
    Extract recommended weight adjustments from trained model.

    seen_params : dict, optional
        Parsed header parameters from the assignment file(s), as returned by
        utils.parse_assignment_header(). When provided, the current
        SPECIFIC_WEIGHT / MAX_SHARED_WEIGHT / DISCORDANT_PENALTY values are
        read from the header so the comparison reflects the actual weights
        used during assignment rather than hardcoded defaults.
    """
    print(f"\n=== WEIGHT RECOMMENDATIONS ===")

    # Baseline weights: prefer values from the assignment file header so the
    # comparison is against what was actually used, not a hardcoded guess.
    if seen_params is None:
        seen_params = {}
    current_specific   = float(seen_params.get('SPECIFIC_WEIGHT',    1.0))
    current_shared     = float(seen_params.get('MAX_SHARED_WEIGHT',  0.2))
    current_discordant = float(seen_params.get('DISCORDANT_PENALTY', -0.10))

    if model_type == 'logistic':
        coef = model.coef_[0]

        # Find coefficients for key features
        feature_dict = dict(zip(feature_names, coef))

        print("\nLogistic Regression Coefficients (higher = more important for correct assignment):")
        print(f"  Specific BC fraction:    {feature_dict.get('specific_frac', 0):+.4f}")
        print(f"  Shared BC fraction:      {feature_dict.get('shared_frac', 0):+.4f}")
        print(f"  Discordant BC fraction:  {feature_dict.get('discordant_frac', 0):+.4f}")
        print(f"  Specific - Discordant:   {feature_dict.get('specific_minus_discordant', 0):+.4f}")
        print(f"  True match fraction:       {feature_dict.get('true_match_frac', 0):+.4f}")

        # Translate to weight recommendations
        print("\nSuggested Weight Adjustments:")

        # Scale based on relative coefficients
        specific_coef   = feature_dict.get('specific_frac', 1.0)
        shared_coef     = feature_dict.get('shared_frac', 0.5)
        discordant_coef = feature_dict.get('discordant_frac', -1.0)

        # Normalize to keep specific at current value as reference
        scale = abs(specific_coef) if abs(specific_coef) > 0 else 1.0

        recommended = {
            'SPECIFIC_WEIGHT':    current_specific,  # Keep as reference
            'MAX_SHARED_WEIGHT':  max(0.01, abs(shared_coef / scale) * current_shared),
            'DISCORDANT_PENALTY': min(-0.01, (discordant_coef / scale) * abs(current_discordant)),
            'CONFIDENCE': 'medium'
        }

        print(f"  SPECIFIC_WEIGHT:    {recommended['SPECIFIC_WEIGHT']:.2f} (unchanged, reference)")
        print(f"  MAX_SHARED_WEIGHT:  {recommended['MAX_SHARED_WEIGHT']:.4f} (current: {current_shared})")
        print(f"  DISCORDANT_PENALTY: {recommended['DISCORDANT_PENALTY']:.4f} (current: {current_discordant})")

    else:  # Random Forest
        # For RF, we can't directly extract weight recommendations
        # But we can use feature importance
        importance = model.feature_importances_
        feature_importance = dict(zip(feature_names, importance))

        print("\nRandom Forest Feature Importance:")
        print(f"  Specific BC fraction:    {feature_importance.get('specific_frac', 0):.4f}")
        print(f"  Shared BC fraction:      {feature_importance.get('shared_frac', 0):.4f}")
        print(f"  Discordant BC fraction:  {feature_importance.get('discordant_frac', 0):.4f}")

        recommended = {
            'SPECIFIC_WEIGHT':    current_specific,
            'MAX_SHARED_WEIGHT':  current_shared,
            'DISCORDANT_PENALTY': current_discordant,
            'CONFIDENCE': 'low',
            'NOTE': 'Random Forest does not provide direct weight translations. Use logistic regression for interpretable weights.'
        }

    return recommended


def plot_model_diagnostics(model, feature_names, training_data, model_type, output_prefix):
    """
    2×2 model diagnostic figure (original plot, now a separate function):
      [0,0] Feature importance / logistic coefficients
      [0,1] Posterior distribution coloured by correctness
      [1,0] Error breakdown by classification tier
      [1,1] Model calibration curve
    """
    import matplotlib.gridspec as gridspec

    X = np.array([[d['features'][f] for f in feature_names] for d in training_data])
    y = np.array([d['label'] for d in training_data])
    probs = model.predict_proba(X)[:, 1]

    fig = plt.figure(figsize=(15, 11))
    fig.suptitle(
        f"Model Diagnostics — {model_type.replace('_', ' ').title()} "
        f"(n={len(training_data):,})",
        fontsize=13, fontweight='bold'
    )
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.35)

    # ── [0,0] Feature importance ──────────────────────────────────────────────
    ax = fig.add_subplot(gs[0, 0])
    if model_type == 'logistic':
        importance = model.coef_[0]
        title = 'Logistic Regression Coefficients\n(positive → correct assignment)'
        colors = ['steelblue' if v >= 0 else 'tomato' for v in importance]
    else:
        importance = model.feature_importances_
        title = 'Random Forest Feature Importance'
        colors = 'steelblue'

    order = np.argsort(np.abs(importance))[::-1][:12]
    labels = [feature_names[i].replace('_', '\n') for i in order]
    vals   = [importance[i] for i in order]
    bar_colors = [colors[i] for i in order] if isinstance(colors, list) else colors

    ax.barh(range(len(vals)), vals[::-1],
            color=bar_colors[::-1] if isinstance(bar_colors, list) else bar_colors)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels[::-1], fontsize=7)
    ax.axvline(0, color='black', linewidth=0.7)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('Coefficient / Importance', fontsize=8)
    ax.grid(True, axis='x', alpha=0.3)

    # ── [0,1] Posterior distribution by correctness ───────────────────────────
    ax = fig.add_subplot(gs[0, 1])
    posteriors = np.array([d['posterior'] for d in training_data])
    bins = np.linspace(0, 1, 50)
    ax.hist(posteriors[y == 1], bins=bins, alpha=0.55, color='steelblue',
            label='Correct', density=True)
    ax.hist(posteriors[y == 0], bins=bins, alpha=0.55, color='tomato',
            label='Incorrect', density=True)
    ax.axvline(0.840, color='green',  linestyle='--', linewidth=1.2,
               label='HC threshold (0.840)')
    ax.axvline(0.500, color='orange', linestyle='--', linewidth=1.2,
               label='LC threshold (0.500)')
    ax.set_xlabel('Current Posterior Score', fontsize=8)
    ax.set_ylabel('Density', fontsize=8)
    ax.set_title('Posterior Distribution by Correctness', fontsize=9)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # ── [1,0] Error breakdown by classification tier ──────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    tiers = ['HIGH_CONF', 'LOW_CONF']
    correct_counts, incorrect_counts = [], []
    for tier in tiers:
        mask = np.array([d['classification'] == tier for d in training_data])
        correct_counts.append(int((y[mask] == 1).sum()))
        incorrect_counts.append(int((y[mask] == 0).sum()))

    x_pos = np.arange(len(tiers))
    w = 0.35
    ax.bar(x_pos - w/2, correct_counts,   w, label='Correct',   color='steelblue', alpha=0.8)
    ax.bar(x_pos + w/2, incorrect_counts, w, label='Incorrect', color='tomato',    alpha=0.8)
    for i, (c, e) in enumerate(zip(correct_counts, incorrect_counts)):
        total = c + e
        acc = c / total * 100 if total > 0 else 0
        ax.text(i, max(c, e) + total * 0.01, f'{acc:.1f}%',
                ha='center', va='bottom', fontsize=8, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(tiers, fontsize=9)
    ax.set_ylabel('ZMW count', fontsize=8)
    ax.set_title('Correct vs Incorrect by Tier', fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(True, axis='y', alpha=0.3)

    # ── [1,1] Model calibration ───────────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    bin_edges = np.linspace(0, 1, 11)
    bin_centres, bin_acc, bin_n = [], [], []
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        mask = (probs >= lo) & (probs < hi)
        if mask.sum() > 0:
            bin_centres.append((lo + hi) / 2)
            bin_acc.append(y[mask].mean())
            bin_n.append(mask.sum())
    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5, label='Perfect calibration')
    sc = ax.scatter(bin_centres, bin_acc, c=bin_n, cmap='Blues',
                    s=80, zorder=5, edgecolors='steelblue', linewidths=0.8)
    plt.colorbar(sc, ax=ax, label='n ZMWs', fraction=0.046, pad=0.04)
    ax.set_xlabel('Model Predicted Probability (correct)', fontsize=8)
    ax.set_ylabel('Actual Fraction Correct', fontsize=8)
    ax.set_title('Model Calibration', fontsize=9)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    out_path = f'{output_prefix}_model_diagnostics.png'
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_weight_sensitivity(training_data, recommended_weights, output_prefix):
    """
    Plot 1: Weight sensitivity sweep.

    Simulate what happens to accuracy and yield as SPECIFIC_WEIGHT,
    MAX_SHARED_WEIGHT, and DISCORDANT_PENALTY each vary around their
    recommended values, holding the other two fixed.  This shows the
    shape of the loss landscape and how robust the recommendation is.

    We approximate the Bayesian score for each ZMW directly from the
    stored specific/shared/discordant counts rather than re-running the full scorer,
    which is fast and sufficient for relative comparison.
    """
    rec = recommended_weights
    inf_ref   = rec.get('SPECIFIC_WEIGHT', 1.0)
    uninf_ref = rec.get('MAX_SHARED_WEIGHT', 0.2)
    extr_ref  = rec.get('DISCORDANT_PENALTY', -0.10)

    # Build a minimal per-ZMW table of counts and ground truth
    records = []
    for d in training_data:
        records.append({
            'n_specific':   d['n_specific'],
            'n_shared':     d['n_shared'],
            'n_discordant': d['n_discordant'],
            'correct':      d['label'],
            'cls':          d['classification'],
        })
    df = pd.DataFrame(records)
    n_total = len(df)

    def score_and_eval(inf_w, uninf_w, extr_w):
        """Return (accuracy_hc, pct_kept_hc) at current HC threshold using given weights."""
        raw = df['n_specific'] * inf_w + df['n_shared'] * uninf_w + df['n_discordant'] * extr_w
        # Approx posterior: softmax with 1 competitor at score=0 (worst case)
        # posterior ≈ exp(raw) / (exp(raw) + exp(0))  → sigmoid
        posterior = 1 / (1 + np.exp(-raw))
        hc_mask = (posterior >= 0.840) & (df['cls'] == 'HIGH_CONF')
        kept = df[hc_mask]
        if len(kept) == 0:
            return np.nan, 0.0
        return kept['correct'].mean() * 100, len(kept) / n_total * 100

    # Sweep ranges: ±50% around recommended, 30 steps each
    sweeps = {
        'SPECIFIC_WEIGHT':    np.linspace(max(0.1, inf_ref * 0.5),   inf_ref * 1.5,   30),
        'MAX_SHARED_WEIGHT':  np.linspace(max(0.01, uninf_ref * 0.5), uninf_ref * 1.5, 30),
        'DISCORDANT_PENALTY': np.linspace(extr_ref * 1.5,  min(-0.01, extr_ref * 0.5), 30),
    }

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle('Weight Sensitivity Sweep\n'
                 '(each weight varied ±50% while others held at recommended value)',
                 fontsize=12, fontweight='bold')

    param_labels = {
        'SPECIFIC_WEIGHT':    ('SPECIFIC_WEIGHT',    inf_ref,   'steelblue'),
        'MAX_SHARED_WEIGHT':  ('MAX_SHARED_WEIGHT',  uninf_ref, 'darkorange'),
        'DISCORDANT_PENALTY': ('DISCORDANT_PENALTY', extr_ref,  'firebrick'),
    }

    for col, (param, (label, ref_val, color)) in enumerate(param_labels.items()):
        vals = sweeps[param]
        accs, yields = [], []
        for v in vals:
            iw = v        if param == 'SPECIFIC_WEIGHT'    else inf_ref
            uw = v        if param == 'MAX_SHARED_WEIGHT'  else uninf_ref
            ew = v        if param == 'DISCORDANT_PENALTY' else extr_ref
            acc, yld = score_and_eval(iw, uw, ew)
            accs.append(acc)
            yields.append(yld)

        accs   = np.array(accs)
        yields = np.array(yields)

        # Accuracy row
        ax = axes[0, col]
        ax.plot(vals, accs, color=color, linewidth=2)
        ax.axvline(ref_val, color='black', linestyle='--', linewidth=1,
                   label=f'Recommended ({ref_val:.3f})')
        ax.set_xlabel(label, fontsize=9)
        ax.set_ylabel('HC Accuracy (%)', fontsize=9)
        ax.set_title(f'Accuracy vs {label}', fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

        # Yield row
        ax = axes[1, col]
        ax.plot(vals, yields, color=color, linewidth=2)
        ax.axvline(ref_val, color='black', linestyle='--', linewidth=1)
        ax.set_xlabel(label, fontsize=9)
        ax.set_ylabel('HC Yield (%)', fontsize=9)
        ax.set_title(f'Yield vs {label}', fontsize=10)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = f'{output_prefix}_weight_sensitivity.png'
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_error_anatomy(training_data, output_prefix):
    """
    Plot 2: Error anatomy breakdown per library/pool.

    For each source pool, show a stacked bar of HIGH_CONF errors broken
    down by likely cause:
      - Low specific count     (n_specific < 2)
      - High discordant count  (n_discordant > n_specific)
      - Ambiguous (shared dominates: n_shared > n_specific and n_discordant <= n_specific)
      - Other / unexplained

    Correct HIGH_CONF calls shown as a separate bar for reference scale.
    """
    # Categorise errors
    records = []
    for d in training_data:
        if d['classification'] != 'HIGH_CONF':
            continue
        n_specific, n_shared, n_discordant = d['n_specific'], d['n_shared'], d['n_discordant']
        correct = bool(d['label'])
        if correct:
            cause = 'correct'
        elif n_specific < 2:
            cause = 'low_specific'
        elif n_discordant > n_specific:
            cause = 'discordant_dominated'
        elif n_shared > n_specific:
            cause = 'shared_dominated'
        else:
            cause = 'other'
        records.append({'pool': d['pool'], 'cause': cause})

    if not records:
        print("  Skipping error anatomy plot: no HIGH_CONF data")
        return

    df = pd.DataFrame(records)
    pools = sorted(df['pool'].unique())
    causes = ['correct', 'low_specific', 'discordant_dominated', 'shared_dominated', 'other']
    cause_colors = {
        'correct':              'steelblue',
        'low_specific':         '#F4A261',
        'discordant_dominated': '#E63946',
        'shared_dominated':     '#9B5DE5',
        'other':                '#AAAAAA',
    }
    cause_labels = {
        'correct':              'Correct',
        'low_specific':         'Error: low specific count (n_specific < 2)',
        'discordant_dominated': 'Error: discordant dominant (n_discordant > n_specific)',
        'shared_dominated':     'Error: shared dominant (n_shared > n_specific)',
        'other':                'Error: other',
    }

    counts = df.groupby(['pool', 'cause']).size().unstack(fill_value=0)
    for c in causes:
        if c not in counts.columns:
            counts[c] = 0
    counts = counts[causes].reindex(pools, fill_value=0)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('HIGH_CONF Error Anatomy by Pool', fontsize=13, fontweight='bold')

    # Absolute counts
    ax = axes[0]
    bottom = np.zeros(len(pools))
    x_pos = np.arange(len(pools))
    for cause in causes:
        vals = counts[cause].values.astype(float)
        ax.bar(x_pos, vals, bottom=bottom, color=cause_colors[cause],
               label=cause_labels[cause], alpha=0.85)
        bottom += vals
    ax.set_xticks(x_pos)
    ax.set_xticklabels([p.replace('_KN_', '\n') for p in pools], fontsize=8)
    ax.set_ylabel('ZMW count')
    ax.set_title('Absolute counts')
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(True, axis='y', alpha=0.3)

    # Fraction of errors only (exclude correct bar)
    ax = axes[1]
    error_causes = [c for c in causes if c != 'correct']
    err_counts = counts[error_causes].copy()
    err_totals = err_counts.sum(axis=1).replace(0, np.nan)
    err_frac = err_counts.div(err_totals, axis=0).fillna(0)

    bottom = np.zeros(len(pools))
    for cause in error_causes:
        vals = err_frac[cause].values
        bars = ax.bar(x_pos, vals, bottom=bottom, color=cause_colors[cause],
                      label=cause_labels[cause], alpha=0.85)
        # Annotate non-trivial fractions
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v > 0.08:
                ax.text(xi, b + v / 2, f'{v*100:.0f}%',
                        ha='center', va='center', fontsize=7, color='white', fontweight='bold')
        bottom += vals
    # Annotate total error count above each bar
    for xi, total in enumerate(err_totals.fillna(0).astype(int)):
        ax.text(xi, 1.01, f'n={total}', ha='center', va='bottom', fontsize=7)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([p.replace('_KN_', '\n') for p in pools], fontsize=8)
    ax.set_ylabel('Fraction of errors')
    ax.set_ylim(0, 1.12)
    ax.set_title('Error type breakdown (fraction)')
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    out_path = f'{output_prefix}_error_anatomy.png'
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_specific_shared_scatter(training_data, output_prefix):
    """
    Plot 3: Specific vs shared barcode count scatter.

    Each ZMW is a point at (n_shared, n_specific), coloured by correctness.
    Incorrect assignments clustering in the high-shared / low-specific region
    indicate MAX_SHARED_WEIGHT is too permissive.

    Separate panels for HIGH_CONF and LOW_CONF.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Specific vs Shared Barcode Counts by Correctness',
                 fontsize=12, fontweight='bold')

    for ax, tier in zip(axes, ['HIGH_CONF', 'LOW_CONF']):
        sub = [d for d in training_data if d['classification'] == tier]
        if not sub:
            ax.set_visible(False)
            continue

        n_shared   = np.array([d['n_shared']   for d in sub])
        n_specific = np.array([d['n_specific'] for d in sub])
        correct    = np.array([d['label']      for d in sub], dtype=bool)

        # Jitter slightly so overlapping points are visible
        jitter = 0.15
        rng = np.random.default_rng(42)
        jx = rng.uniform(-jitter, jitter, size=len(n_shared))
        jy = rng.uniform(-jitter, jitter, size=len(n_specific))

        ax.scatter(n_shared[correct]    + jx[correct],
                   n_specific[correct]  + jy[correct],
                   c='steelblue', alpha=0.25, s=8, label='Correct', rasterized=True)
        ax.scatter(n_shared[~correct]   + jx[~correct],
                   n_specific[~correct] + jy[~correct],
                   c='tomato', alpha=0.6, s=12, label='Incorrect', rasterized=True)

        # Reference lines: y = x (specific == shared) and y = 2x
        max_val = max(n_shared.max(), n_specific.max()) + 1
        xs = np.array([0, max_val])
        ax.plot(xs, xs,     'k--', linewidth=0.8, alpha=0.5, label='specific = shared')
        ax.plot(xs, xs * 2, 'k:',  linewidth=0.8, alpha=0.4, label='specific = 2× shared')

        n_err = int((~correct).sum())
        n_tot = len(sub)
        ax.set_xlabel('Shared Barcode Count', fontsize=9)
        ax.set_ylabel('Specific Barcode Count', fontsize=9)
        ax.set_title(f'{tier}  (n={n_tot:,}, {n_err} errors)', fontsize=10)
        ax.legend(fontsize=8, markerscale=2)
        ax.grid(True, alpha=0.2)

    plt.tight_layout()
    out_path = f'{output_prefix}_specific_shared_scatter.png'
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_per_pool_accuracy(training_data, output_prefix):
    """
    Plot 4: Per-pool accuracy summary.

    Grouped bar chart: for each training pool, show accuracy for
    HIGH_CONF and LOW_CONF separately, with n= annotated.
    Lets you see whether weight changes help representative pools
    or only stress-test pools.
    """
    records = []
    for d in training_data:
        records.append({
            'pool': d['pool'],
            'cls':  d['classification'],
            'correct': d['label'],
        })
    df = pd.DataFrame(records)

    pools = sorted(df['pool'].unique())
    tiers = ['HIGH_CONF', 'LOW_CONF']
    tier_colors = {'HIGH_CONF': 'steelblue', 'LOW_CONF': 'darkorange'}

    x_pos  = np.arange(len(pools))
    width  = 0.35
    offsets = {'HIGH_CONF': -width / 2, 'LOW_CONF': width / 2}

    fig, ax = plt.subplots(figsize=(max(8, len(pools) * 1.8), 6))
    fig.suptitle('Per-Pool Accuracy by Confidence Tier', fontsize=12, fontweight='bold')

    for tier in tiers:
        accs, ns = [], []
        for pool in pools:
            sub = df[(df['pool'] == pool) & (df['cls'] == tier)]
            accs.append(sub['correct'].mean() * 100 if len(sub) > 0 else np.nan)
            ns.append(len(sub))

        bars = ax.bar(x_pos + offsets[tier], accs, width,
                      label=tier, color=tier_colors[tier], alpha=0.8)
        for xi, (acc, n) in enumerate(zip(accs, ns)):
            if not np.isnan(acc):
                ax.text(xi + offsets[tier], acc + 0.2, f'{acc:.1f}%\nn={n:,}',
                        ha='center', va='bottom', fontsize=7)

    ax.axhline(99.9, color='red',  linestyle='--', linewidth=1, alpha=0.6, label='99.9% target')
    ax.axhline(95.0, color='gray', linestyle=':', linewidth=1, alpha=0.6, label='95% (LC target)')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([p.replace('_KN_', '\n') for p in pools], fontsize=8)
    ax.set_ylabel('Accuracy (%)')
    ax.set_ylim([max(0, df.groupby(['pool','cls'])['correct'].mean().min() * 100 - 5), 101])
    ax.legend(fontsize=9)
    ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    out_path = f'{output_prefix}_per_pool_accuracy.png'
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_score_separation(training_data, output_prefix):
    """
    Plot 5: Score separation — posterior gap between top and second-best library.

    The posterior of the best library is p1; by the softmax identity, the
    remaining probability mass (1 - p1) is shared among all other libraries.
    With N libraries the expected second-best is (1 - p1) / (N - 1), giving:

        gap = p1 - (1 - p1) / (N - 1)

    With only 2 libraries this simplifies to gap = 2*p1 - 1.

    A well-tuned weight set should produce a bimodal gap distribution:
    correct assignments cluster at large gap, errors at small gap.
    """
    # Infer N_libraries per ZMW from the pool: count unique 'true' values per pool
    pool_n_libs: dict[str, int] = {}
    from collections import Counter
    pool_libs_counter: dict[str, set] = defaultdict(set)
    for d in training_data:
        pool_libs_counter[d['pool']].add(d['true'])
    for pool, libs in pool_libs_counter.items():
        pool_n_libs[pool] = max(2, len(libs))  # at least 2 to avoid div/0

    gaps, correct_arr, tier_arr = [], [], []
    for d in training_data:
        p1 = d['posterior']
        n  = pool_n_libs.get(d['pool'], 2)
        gap = p1 - (1 - p1) / max(n - 1, 1)
        gaps.append(gap)
        correct_arr.append(d['label'])
        tier_arr.append(d['classification'])

    gaps      = np.array(gaps)
    correct   = np.array(correct_arr, dtype=bool)
    tiers     = np.array(tier_arr)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle('Score Separation: Posterior Gap (top vs second-best library)',
                 fontsize=12, fontweight='bold')

    bins = np.linspace(gaps.min(), 1.0, 60)

    for ax, tier in zip(axes, ['HIGH_CONF', 'LOW_CONF']):
        mask = tiers == tier
        g_corr = gaps[mask &  correct]
        g_err  = gaps[mask & ~correct]

        ax.hist(g_corr, bins=bins, alpha=0.55, color='steelblue',
                label=f'Correct (n={len(g_corr):,})',   density=True)
        ax.hist(g_err,  bins=bins, alpha=0.65, color='tomato',
                label=f'Incorrect (n={len(g_err):,})', density=True)

        # Mark median of each group
        if len(g_corr) > 0:
            ax.axvline(np.median(g_corr), color='steelblue', linestyle='--',
                       linewidth=1.5, label=f'Median correct ({np.median(g_corr):.2f})')
        if len(g_err) > 0:
            ax.axvline(np.median(g_err), color='tomato', linestyle='--',
                       linewidth=1.5, label=f'Median incorrect ({np.median(g_err):.2f})')

        ax.set_xlabel('Posterior Gap  p₁ − (1−p₁)/(N−1)', fontsize=9)
        ax.set_ylabel('Density', fontsize=9)
        ax.set_title(tier, fontsize=11)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out_path = f'{output_prefix}_score_separation.png'
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Optimize barcode weights using machine learning",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--arrays', required=True, help='Arrays definition file (arrays.txt)')
    parser.add_argument('--assignments', nargs='+', required=True, 
                        help='Assignment result files (can use wildcards)')
    parser.add_argument('--lookups', nargs='+', required=True,
                        help='Ground truth lookup files (can use wildcards)')
    parser.add_argument('--output', default='optimized_weights.json',
                        help='Output file for recommended weights')
    parser.add_argument('--plot-prefix', default=None,
                        help='Prefix for output plot (e.g. results/weights). '
                             'Produces <prefix>_weight_analysis.png. '
                             'Defaults to the output path with .json stripped.')
    parser.add_argument('--model', choices=['logistic', 'random_forest'], default='logistic',
                        help='ML model type (logistic is more interpretable)')
    
    args = parser.parse_args()
    
    # Expand wildcards
    assignment_files = []
    for pattern in args.assignments:
        assignment_files.extend(Path('.').glob(pattern))
    
    lookup_files = []
    for pattern in args.lookups:
        lookup_files.extend(Path('.').glob(pattern))
    
    if not assignment_files:
        print(f"ERROR: No assignment files found matching {args.assignments}")
        sys.exit(1)
    
    if not lookup_files:
        print(f"ERROR: No lookup files found matching {args.lookups}")
        sys.exit(1)
    
    print(f"Found {len(assignment_files)} assignment files")
    print(f"Found {len(lookup_files)} lookup files")
    
    # Load data
    arrays = utils.read_arrays(args.arrays)
    print(f"Loaded {len(arrays)} array definitions")
    
    training_data, seen_params = load_training_data(assignment_files, lookup_files, arrays)
    
    if len(training_data) < 100:
        print(f"WARNING: Only {len(training_data)} training examples. Need at least 100 for reliable training.")
        sys.exit(1)
    
    # Train model
    model, feature_names, metrics = train_model(training_data, model_type=args.model)
    
    # Analyze errors
    analyze_errors(model, feature_names, training_data, arrays)
    
    # Extract recommendations
    recommended_weights = extract_weight_recommendations(model, feature_names, args.model, seen_params)

    # Plot diagnostic figures
    plot_prefix = args.plot_prefix or str(Path(args.output).with_suffix(''))
    print("\nGenerating plots...")
    plot_model_diagnostics(model, feature_names, training_data, args.model, plot_prefix)
    plot_weight_sensitivity(training_data, recommended_weights, plot_prefix)
    plot_error_anatomy(training_data, plot_prefix)
    plot_specific_shared_scatter(training_data, plot_prefix)
    plot_per_pool_accuracy(training_data, plot_prefix)
    plot_score_separation(training_data, plot_prefix)

    # Save results
    output_data = {
        'recommended_weights': recommended_weights,
        'model_type': args.model,
        'performance': metrics,
        'training_examples': len(training_data),
        'feature_importance': dict(zip(feature_names, 
            np.abs(model.coef_[0]) if args.model == 'logistic' else model.feature_importances_))
    }
    
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nResults saved to {args.output}")
    print(f"\nTo use these weights, update your assign_kinnex script:")
    print(f"  SPECIFIC_WEIGHT    = {recommended_weights['SPECIFIC_WEIGHT']}")
    print(f"  MAX_SHARED_WEIGHT  = {recommended_weights['MAX_SHARED_WEIGHT']:.3f}")
    print(f"  # And change line with 'score -= 1.0' to:")
    print(f"  score += {recommended_weights['DISCORDANT_PENALTY']:.3f}  # (was -1.0)")


if __name__ == '__main__':
    main()

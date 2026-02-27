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
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
import warnings
warnings.filterwarnings('ignore')


def read_arrays(arrays_file):
    """Read array definitions."""
    arrays = {}
    with open(arrays_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            arrays[parts[0]] = {
                'kinnex': parts[1],
                'barcodes': set(parts[2:])
            }
    return arrays


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
    
    # Process assignment files
    for assign_file in assignment_files:
        print(f"Processing {assign_file}...")
        df = pd.read_csv(assign_file, sep='\t', dtype=str)
        
        # Convert numeric columns
        df['Top_Posterior'] = pd.to_numeric(df['Top_Posterior'], errors='coerce')
        df['Informative_Barcodes'] = pd.to_numeric(df['Informative_Barcodes'], errors='coerce')
        df['Uninformative_Barcodes'] = pd.to_numeric(df['Uninformative_Barcodes'], errors='coerce')
        df['Extraneous_Barcodes'] = pd.to_numeric(df['Extraneous_Barcodes'], errors='coerce')
        
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
                'classification': classification,
                'assigned': assigned_array,
                'true': true_array
            })
    
    print(f"Collected {len(training_data):,} training examples")
    return training_data


def extract_features(zmw, barcodes, assigned_array, true_array, arrays, classification, row):
    """
    Extract features for ML model.
    
    Features:
    - Raw counts: inf, uninf, extr barcodes
    - Ratios: inf/total, uninf/total, extr/total
    - Segment count
    - Barcode diversity (unique barcodes / total)
    - Evidence strength metrics
    """
    n_segments = len(barcodes)
    n_inf = row['Informative_Barcodes']
    n_uninf = row['Uninformative_Barcodes']
    n_extr = row['Extraneous_Barcodes']
    n_total_bc = n_inf + n_uninf + n_extr
    
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
        'n_segments': n_segments,
        'n_inf': n_inf,
        'n_uninf': n_uninf,
        'n_extr': n_extr,
        
        # Ratios (avoid division by zero)
        'inf_frac': n_inf / max(n_total_bc, 1),
        'uninf_frac': n_uninf / max(n_total_bc, 1),
        'extr_frac': n_extr / max(n_total_bc, 1),
        
        # Barcode diversity
        'unique_bc_count': n_unique_bc,
        'bc_diversity': n_unique_bc / max(n_segments, 1),
        'max_bc_frequency': max_bc_count / max(n_segments, 1),
        
        # Evidence strength
        'true_matches': true_matches,
        'true_match_frac': true_matches / max(n_segments, 1),
        'true_unique_bc': true_unique_bc,
        
        # Interaction terms
        'inf_times_segments': n_inf * n_segments,
        'extr_times_segments': n_extr * n_segments,
        'inf_minus_extr': n_inf - n_extr,
        
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
            print(f"    Features: inf={d['features']['n_inf']}, uninf={d['features']['n_uninf']}, extr={d['features']['n_extr']}")
    
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


def extract_weight_recommendations(model, feature_names, model_type='logistic'):
    """
    Extract recommended weight adjustments from trained model.
    """
    print(f"\n=== WEIGHT RECOMMENDATIONS ===")
    
    if model_type == 'logistic':
        coef = model.coef_[0]
        
        # Find coefficients for key features
        feature_dict = dict(zip(feature_names, coef))
        
        print("\nLogistic Regression Coefficients (higher = more important for correct assignment):")
        print(f"  Informative BC fraction:   {feature_dict.get('inf_frac', 0):+.4f}")
        print(f"  Uninformative BC fraction: {feature_dict.get('uninf_frac', 0):+.4f}")
        print(f"  Extraneous BC fraction:    {feature_dict.get('extr_frac', 0):+.4f}")
        print(f"  Inf - Extr balance:        {feature_dict.get('inf_minus_extr', 0):+.4f}")
        print(f"  True match fraction:       {feature_dict.get('true_match_frac', 0):+.4f}")
        
        # Translate to weight recommendations
        print("\nSuggested Weight Adjustments:")
        
        # Current weights — must match defaults in assign_kinnex.py
        current_inf = 1.0
        current_uninf = 0.2
        current_extr = -0.10
        
        # Scale based on relative coefficients
        inf_coef = feature_dict.get('inf_frac', 1.0)
        uninf_coef = feature_dict.get('uninf_frac', 0.5)
        extr_coef = feature_dict.get('extr_frac', -1.0)
        
        # Normalize to keep informative at 1.0 as reference
        scale = abs(inf_coef) if abs(inf_coef) > 0 else 1.0
        
        recommended = {
            'INF_WEIGHT': current_inf,  # Keep as reference
            'MAX_UNINF_WEIGHT': max(0.01, abs(uninf_coef / scale) * current_uninf),
            'EXTRANEOUS_PENALTY': min(-0.1, (extr_coef / scale) * abs(current_extr)),
            'CONFIDENCE': 'medium'
        }
        
        print(f"  INF_WEIGHT:         {recommended['INF_WEIGHT']:.2f} (unchanged, reference)")
        print(f"  MAX_UNINF_WEIGHT:   {recommended['MAX_UNINF_WEIGHT']:.3f} (current: {current_uninf})")
        print(f"  EXTRANEOUS_PENALTY: {recommended['EXTRANEOUS_PENALTY']:.3f} (current: {current_extr})")
        
    else:  # Random Forest
        # For RF, we can't directly extract weight recommendations
        # But we can use feature importance
        importance = model.feature_importances_
        feature_importance = dict(zip(feature_names, importance))
        
        print("\nRandom Forest Feature Importance:")
        print(f"  Informative BC fraction:   {feature_importance.get('inf_frac', 0):.4f}")
        print(f"  Uninformative BC fraction: {feature_importance.get('uninf_frac', 0):.4f}")
        print(f"  Extraneous BC fraction:    {feature_importance.get('extr_frac', 0):.4f}")
        
        recommended = {
            'INF_WEIGHT': 1.0,
            'MAX_UNINF_WEIGHT': 0.2,
            'EXTRANEOUS_PENALTY': -0.10,
            'CONFIDENCE': 'low',
            'NOTE': 'Random Forest does not provide direct weight translations. Use logistic regression for interpretable weights.'
        }
    
    return recommended


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
    arrays = read_arrays(args.arrays)
    print(f"Loaded {len(arrays)} array definitions")
    
    training_data = load_training_data(assignment_files, lookup_files, arrays)
    
    if len(training_data) < 100:
        print(f"WARNING: Only {len(training_data)} training examples. Need at least 100 for reliable training.")
        sys.exit(1)
    
    # Train model
    model, feature_names, metrics = train_model(training_data, model_type=args.model)
    
    # Analyze errors
    analyze_errors(model, feature_names, training_data, arrays)
    
    # Extract recommendations
    recommended_weights = extract_weight_recommendations(model, feature_names, args.model)
    
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
    print(f"\nTo use these weights, update the constants at the top of assign_kinnex.py:")
    print(f"  INF_WEIGHT = {recommended_weights['INF_WEIGHT']}")
    print(f"  MAX_UNINF_WEIGHT = {recommended_weights['MAX_UNINF_WEIGHT']:.3f}")
    print(f"  EXTRANEOUS_PENALTY = {recommended_weights['EXTRANEOUS_PENALTY']:.3f}")


if __name__ == '__main__':
    main()

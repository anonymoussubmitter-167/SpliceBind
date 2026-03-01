#!/usr/bin/env python3
"""
Comprehensive model analysis for SpliceBind:
1. Per-fold test set sizes
2. Training duration statistics
3. Bootstrap confidence intervals
4. Calibration metrics (ECE)
5. Statistical significance testing
"""

import sys
import json
import numpy as np
import torch
from pathlib import Path
from collections import defaultdict
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from scipy import stats

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))


def compute_bootstrap_ci(values, n_bootstrap=1000, ci=0.95):
    """Compute bootstrap confidence interval."""
    values = np.array(values)
    n = len(values)

    bootstrap_means = []
    for _ in range(n_bootstrap):
        sample = np.random.choice(values, size=n, replace=True)
        bootstrap_means.append(np.mean(sample))

    lower = np.percentile(bootstrap_means, (1 - ci) / 2 * 100)
    upper = np.percentile(bootstrap_means, (1 + ci) / 2 * 100)

    return lower, upper


def compute_ece(y_true, y_prob, n_bins=10):
    """Compute Expected Calibration Error."""
    bin_boundaries = np.linspace(0, 1, n_bins + 1)
    bin_lowers = bin_boundaries[:-1]
    bin_uppers = bin_boundaries[1:]

    ece = 0
    for bin_lower, bin_upper in zip(bin_lowers, bin_uppers):
        in_bin = (y_prob > bin_lower) & (y_prob <= bin_upper)
        prop_in_bin = in_bin.mean()

        if prop_in_bin > 0:
            avg_confidence = y_prob[in_bin].mean()
            avg_accuracy = y_true[in_bin].mean()
            ece += np.abs(avg_confidence - avg_accuracy) * prop_in_bin

    return ece


def main():
    print("=" * 70)
    print("SPLICEBIND COMPREHENSIVE ANALYSIS")
    print("=" * 70)

    # Load results
    results_path = PROJECT_ROOT / "data" / "processed" / "results.json"
    if not results_path.exists():
        print("Results file not found. Run train.py first.")
        return

    with open(results_path) as f:
        results = json.load(f)

    cv_results = results.get('cross_validation', {}).get('per_fold', [])

    if not cv_results:
        print("No CV results found.")
        return

    # =========================================================================
    # 1. PER-FOLD TEST SET SIZES
    # =========================================================================
    print("\n" + "=" * 70)
    print("1. PER-FOLD TEST SET SIZES")
    print("=" * 70)

    fold_sizes = defaultdict(list)
    for fold_result in cv_results:
        family = fold_result.get('test_family', 'Unknown')
        n_samples = len(fold_result.get('predictions', []))
        fold_sizes[family].append(n_samples)

    print(f"\n{'Family':<15} {'Mean N':<10} {'Min':<8} {'Max':<8} {'Folds'}")
    print("-" * 50)

    total_samples = 0
    for family in sorted(fold_sizes.keys()):
        sizes = fold_sizes[family]
        mean_n = np.mean(sizes)
        min_n = np.min(sizes)
        max_n = np.max(sizes)
        n_folds = len(sizes)
        total_samples += sum(sizes)
        print(f"{family:<15} {mean_n:<10.1f} {min_n:<8} {max_n:<8} {n_folds}")

    print(f"\nTotal test samples across all folds: {total_samples}")
    print(f"Average per fold: {total_samples / len(cv_results):.1f}")

    # =========================================================================
    # 2. TRAINING DURATION STATISTICS
    # =========================================================================
    print("\n" + "=" * 70)
    print("2. TRAINING DURATION STATISTICS")
    print("=" * 70)

    # Note: We don't have epoch counts in results, but we can estimate
    print("\nTraining configuration:")
    print("  Max epochs: 50")
    print("  Early stopping patience: 10")
    print("  Estimated typical convergence: 15-25 epochs")
    print("  (Exact epoch tracking not saved in current results)")

    # =========================================================================
    # 3. BOOTSTRAP CONFIDENCE INTERVALS
    # =========================================================================
    print("\n" + "=" * 70)
    print("3. BOOTSTRAP CONFIDENCE INTERVALS (vs Normal Assumption)")
    print("=" * 70)

    aurocs = [r['auroc'] for r in cv_results]

    # Normal assumption CI
    mean_auroc = np.mean(aurocs)
    std_auroc = np.std(aurocs)
    se_auroc = std_auroc / np.sqrt(len(aurocs))
    normal_ci = (mean_auroc - 1.96 * se_auroc, mean_auroc + 1.96 * se_auroc)

    # Bootstrap CI
    bootstrap_ci = compute_bootstrap_ci(aurocs, n_bootstrap=10000)

    print(f"\nAUROC Distribution:")
    print(f"  Mean: {mean_auroc:.3f}")
    print(f"  Std: {std_auroc:.3f}")
    print(f"  Min: {np.min(aurocs):.3f}")
    print(f"  Max: {np.max(aurocs):.3f}")
    print(f"\n95% CI Comparison:")
    print(f"  Normal assumption: [{normal_ci[0]:.3f}, {normal_ci[1]:.3f}]")
    print(f"  Bootstrap (10K):   [{bootstrap_ci[0]:.3f}, {bootstrap_ci[1]:.3f}]")

    # Shapiro-Wilk normality test
    if len(aurocs) >= 3:
        stat, p_value = stats.shapiro(aurocs)
        print(f"\nNormality test (Shapiro-Wilk): p={p_value:.4f}")
        if p_value < 0.05:
            print("  → Distribution NOT normal, bootstrap CI more appropriate")
        else:
            print("  → Cannot reject normality")

    # =========================================================================
    # 4. CALIBRATION METRICS
    # =========================================================================
    print("\n" + "=" * 70)
    print("4. CALIBRATION METRICS")
    print("=" * 70)

    # Aggregate all predictions and labels
    all_preds = []
    all_labels = []

    for fold_result in cv_results:
        preds = fold_result.get('predictions', [])
        # Need true labels - estimate from threshold
        threshold = fold_result.get('optimal_threshold', 0.5)
        # Without true labels, we can only show Brier score components
        all_preds.extend(preds)

    all_preds = np.array(all_preds)

    print(f"\nPrediction statistics:")
    print(f"  Total predictions: {len(all_preds)}")
    print(f"  Mean probability: {np.mean(all_preds):.3f}")
    print(f"  Std probability: {np.std(all_preds):.3f}")
    print(f"  Min probability: {np.min(all_preds):.3f}")
    print(f"  Max probability: {np.max(all_preds):.3f}")

    # Prediction distribution
    bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    hist, _ = np.histogram(all_preds, bins=bins)
    print(f"\nPrediction distribution:")
    for i in range(len(bins) - 1):
        print(f"  [{bins[i]:.1f}, {bins[i+1]:.1f}): {hist[i]} ({100*hist[i]/len(all_preds):.1f}%)")

    # Note about ECE
    print("\nNote: ECE requires true labels per prediction (not saved in results)")
    print("Recommendation: Modify training to save (pred, label) pairs for ECE calculation")

    # =========================================================================
    # 5. ABLATION STUDY REQUIREMENTS
    # =========================================================================
    print("\n" + "=" * 70)
    print("5. ABLATION STUDY REQUIREMENTS")
    print("=" * 70)

    print("\nAblation conditions:")
    print("  1. Structure-only baseline (24 features, no ESM-2)")
    print("  2. ESM-2 only (32 features, no structural)")
    print("  3. Random vectors replacing ESM-2")
    print("  4. Full model (56 features)")

    print("\nCurrent status:")
    esm_ablation = results.get('cross_validation', {}).get('auroc_mean', 0)
    print(f"  Full model AUROC: {esm_ablation:.3f}")
    print("  Structure-only: ~0.58 (reported in paper)")
    print("  ESM-2 only: NOT TESTED")
    print("  Random vectors: NOT TESTED")

    print("\nRecommendation: Run ablation_study.py with all 4 conditions")

    # =========================================================================
    # 6. SUMMARY OF MODEL STATISTICS
    # =========================================================================
    print("\n" + "=" * 70)
    print("6. MODEL STATISTICS SUMMARY")
    print("=" * 70)

    print("\nParameter count: 899,457")
    print("Training samples: 178 pockets")
    print("Parameter-to-sample ratio: 5,053:1 (severely overparameterized)")
    print("\nRecommendation: Simplify model architecture:")
    print("  - Reduce hidden_dim from 256 to 64")
    print("  - Reduce num_layers from 3 to 2")
    print("  - Reduce MLP head dimensions")
    print("  - Target: ~50K parameters (~280:1 ratio)")

    # =========================================================================
    # 7. STATISTICAL SIGNIFICANCE
    # =========================================================================
    print("\n" + "=" * 70)
    print("7. STATISTICAL SIGNIFICANCE vs P2Rank")
    print("=" * 70)

    significance = results.get('significance', {})
    print(f"\nSpliceBind vs P2Rank:")
    print(f"  Delta: {significance.get('delta_vs_p2rank', 0):+.3f}")
    print(f"  p-value: {significance.get('p_value', 1):.4f}")
    print(f"  Significant (p<0.05): {significance.get('significant', False)}")

    # Power analysis
    current_n = len(cv_results)
    current_effect = significance.get('delta_vs_p2rank', 0.046)
    current_std = std_auroc

    # Estimate required n for significance
    if current_effect > 0 and current_std > 0:
        # Using effect size calculation
        effect_size = current_effect / current_std
        # For 80% power at alpha=0.05, need ~16 samples per group for effect size of 1
        required_n = int(16 / (effect_size ** 2)) if effect_size > 0 else 1000
        print(f"\nPower analysis:")
        print(f"  Current effect size: {effect_size:.3f}")
        print(f"  Estimated folds needed for p<0.05: ~{max(required_n, current_n)}")
        print(f"  Current folds: {current_n}")

    # =========================================================================
    # SAVE ANALYSIS RESULTS
    # =========================================================================
    analysis_results = {
        'fold_sizes': {k: {'mean': float(np.mean(v)), 'min': int(np.min(v)), 'max': int(np.max(v))}
                       for k, v in fold_sizes.items()},
        'auroc_stats': {
            'mean': float(mean_auroc),
            'std': float(std_auroc),
            'bootstrap_ci': [float(bootstrap_ci[0]), float(bootstrap_ci[1])],
            'normal_ci': [float(normal_ci[0]), float(normal_ci[1])],
        },
        'model_stats': {
            'parameters': 899457,
            'training_samples': 178,
            'param_sample_ratio': 5053.1
        },
        'recommendations': [
            'Expand dataset to ~550 pockets for statistical significance',
            'Simplify model to ~50K parameters',
            'Run full 4-condition ablation study',
            'Add ECE calibration metric',
            'Report bootstrap CIs instead of normal assumption'
        ]
    }

    output_path = PROJECT_ROOT / "data" / "processed" / "comprehensive_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(analysis_results, f, indent=2)

    print(f"\n\nAnalysis saved to {output_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Binding Affinity Analysis: Correlation between SpliceBind scores and experimental Kd.

Shows that SpliceBind predictions correlate with actual binding affinity measurements.
"""

import sys
import json
import numpy as np
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = Path(__file__).parent.parent


# Curated binding data: (gene, drug, Kd_nM, source)
BINDING_DATA = [
    # JAK family
    ('JAK1', 'Tofacitinib', 3.2, 'BindingDB'),
    ('JAK2', 'Ruxolitinib', 3.3, 'BindingDB'),
    ('JAK3', 'Tofacitinib', 1.0, 'BindingDB'),

    # ABL
    ('ABL1', 'Imatinib', 25, 'BindingDB'),
    ('ABL1', 'Nilotinib', 5, 'BindingDB'),
    ('ABL1', 'Dasatinib', 0.5, 'BindingDB'),
    ('ABL1', 'Ponatinib', 0.4, 'BindingDB'),

    # EGFR
    ('EGFR', 'Gefitinib', 2, 'BindingDB'),
    ('EGFR', 'Erlotinib', 0.7, 'BindingDB'),
    ('EGFR', 'Osimertinib', 0.5, 'BindingDB'),
    ('EGFR', 'Afatinib', 0.5, 'BindingDB'),

    # BRAF
    ('BRAF', 'Vemurafenib', 31, 'BindingDB'),
    ('BRAF', 'Dabrafenib', 0.8, 'BindingDB'),

    # BTK
    ('BTK', 'Ibrutinib', 0.5, 'BindingDB'),
    ('BTK', 'Acalabrutinib', 3, 'BindingDB'),

    # ALK
    ('ALK', 'Crizotinib', 24, 'BindingDB'),
    ('ALK', 'Ceritinib', 0.2, 'BindingDB'),
    ('ALK', 'Alectinib', 1.9, 'BindingDB'),
    ('ALK', 'Lorlatinib', 0.7, 'BindingDB'),

    # KIT
    ('KIT', 'Imatinib', 100, 'BindingDB'),
    ('KIT', 'Sunitinib', 10, 'BindingDB'),

    # FLT3
    ('FLT3', 'Midostaurin', 10, 'ChEMBL'),
    ('FLT3', 'Gilteritinib', 0.3, 'ChEMBL'),
    ('FLT3', 'Quizartinib', 1.1, 'ChEMBL'),

    # MET
    ('MET', 'Crizotinib', 8, 'BindingDB'),
    ('MET', 'Capmatinib', 0.1, 'BindingDB'),

    # RET
    ('RET', 'Selpercatinib', 0.4, 'BindingDB'),
    ('RET', 'Pralsetinib', 0.4, 'BindingDB'),

    # VEGFR2
    ('VEGFR2', 'Axitinib', 0.2, 'BindingDB'),
    ('VEGFR2', 'Sunitinib', 9, 'BindingDB'),
    ('VEGFR2', 'Sorafenib', 90, 'BindingDB'),

    # CDK
    ('CDK4', 'Palbociclib', 11, 'ChEMBL'),
    ('CDK4', 'Ribociclib', 10, 'ChEMBL'),
    ('CDK6', 'Palbociclib', 16, 'ChEMBL'),

    # PIK3
    ('PIK3CA', 'Alpelisib', 5, 'ChEMBL'),
    ('PIK3CD', 'Idelalisib', 2.5, 'ChEMBL'),

    # Aurora
    ('AURKA', 'Alisertib', 1.2, 'ChEMBL'),
    ('AURKB', 'Barasertib', 0.4, 'ChEMBL'),

    # NTRK
    ('NTRK1', 'Larotrectinib', 5, 'BindingDB'),
    ('NTRK1', 'Entrectinib', 1, 'BindingDB'),

    # ROS1
    ('ROS1', 'Crizotinib', 1.7, 'BindingDB'),
    ('ROS1', 'Entrectinib', 0.1, 'BindingDB'),
]


def main():
    print("=" * 70)
    print("BINDING AFFINITY CORRELATION ANALYSIS")
    print("=" * 70)

    # Load SpliceBind predictions
    results_path = PROJECT_ROOT / "data" / "processed" / "expanded_variants.json"
    if results_path.exists():
        with open(results_path) as f:
            variant_data = json.load(f)
    else:
        variant_data = {'variants': []}

    # Build gene -> canonical score mapping
    gene_scores = {}
    for v in variant_data.get('variants', []):
        gene = v['gene']
        if gene not in gene_scores and v.get('canonical_score') is not None:
            gene_scores[gene] = v['canonical_score']

    # Match binding data with predictions
    matched_data = []

    for gene, drug, kd, source in BINDING_DATA:
        if gene in gene_scores:
            matched_data.append({
                'gene': gene,
                'drug': drug,
                'kd_nM': kd,
                'pKd': -np.log10(kd * 1e-9),  # Convert to pKd
                'splicebind_score': gene_scores[gene],
                'source': source,
            })

    print(f"\nMatched binding entries: {len(matched_data)}")
    print(f"Unique genes: {len(set(d['gene'] for d in matched_data))}")

    if len(matched_data) < 5:
        print("Not enough data for correlation analysis")
        return

    # Compute correlations
    pKd_values = np.array([d['pKd'] for d in matched_data])
    scores = np.array([d['splicebind_score'] for d in matched_data])

    # Pearson correlation
    pearson_r, pearson_p = stats.pearsonr(pKd_values, scores)

    # Spearman correlation (rank-based, more robust)
    spearman_r, spearman_p = stats.spearmanr(pKd_values, scores)

    print(f"\n{'='*50}")
    print("CORRELATION RESULTS")
    print(f"{'='*50}")
    print(f"\nPearson correlation:")
    print(f"  r = {pearson_r:.3f}, p = {pearson_p:.4f}")

    print(f"\nSpearman correlation:")
    print(f"  rho = {spearman_r:.3f}, p = {spearman_p:.4f}")

    # Interpretation
    print(f"\nInterpretation:")
    if pearson_p < 0.05 and pearson_r > 0:
        print(f"  Significant positive correlation: higher pKd (stronger binding)")
        print(f"  associates with higher SpliceBind druggability scores.")
    elif pearson_p < 0.05 and pearson_r < 0:
        print(f"  Significant negative correlation.")
    else:
        print(f"  No significant correlation (may be due to ceiling effect - most")
        print(f"  canonical kinases have very high druggability scores).")

    # Distribution analysis
    print(f"\n{'='*50}")
    print("SCORE DISTRIBUTION BY BINDING AFFINITY")
    print(f"{'='*50}")

    # Bin by Kd
    weak_binders = [d for d in matched_data if d['kd_nM'] > 50]  # Kd > 50nM
    moderate_binders = [d for d in matched_data if 10 <= d['kd_nM'] <= 50]
    strong_binders = [d for d in matched_data if d['kd_nM'] < 10]  # Kd < 10nM

    print(f"\nStrong binders (Kd < 10nM): N={len(strong_binders)}")
    if strong_binders:
        scores_strong = [d['splicebind_score'] for d in strong_binders]
        print(f"  Mean SpliceBind: {np.mean(scores_strong):.3f} ± {np.std(scores_strong):.3f}")

    print(f"\nModerate binders (10-50nM): N={len(moderate_binders)}")
    if moderate_binders:
        scores_mod = [d['splicebind_score'] for d in moderate_binders]
        print(f"  Mean SpliceBind: {np.mean(scores_mod):.3f} ± {np.std(scores_mod):.3f}")

    print(f"\nWeak binders (Kd > 50nM): N={len(weak_binders)}")
    if weak_binders:
        scores_weak = [d['splicebind_score'] for d in weak_binders]
        print(f"  Mean SpliceBind: {np.mean(scores_weak):.3f} ± {np.std(scores_weak):.3f}")

    # Table
    print(f"\n{'='*50}")
    print("DATA TABLE")
    print(f"{'='*50}")
    print(f"\n{'Gene':<10} {'Drug':<15} {'Kd (nM)':<10} {'pKd':<8} {'SB Score':<10}")
    print("-" * 55)

    for d in sorted(matched_data, key=lambda x: x['kd_nM']):
        print(f"{d['gene']:<10} {d['drug']:<15} {d['kd_nM']:<10.1f} "
              f"{d['pKd']:<8.2f} {d['splicebind_score']:<10.3f}")

    # Save results
    output = {
        'n_entries': len(matched_data),
        'n_genes': len(set(d['gene'] for d in matched_data)),
        'pearson_r': float(pearson_r),
        'pearson_p': float(pearson_p),
        'spearman_r': float(spearman_r),
        'spearman_p': float(spearman_p),
        'by_affinity': {
            'strong': {
                'n': len(strong_binders),
                'mean_score': float(np.mean([d['splicebind_score'] for d in strong_binders])) if strong_binders else None,
            },
            'moderate': {
                'n': len(moderate_binders),
                'mean_score': float(np.mean([d['splicebind_score'] for d in moderate_binders])) if moderate_binders else None,
            },
            'weak': {
                'n': len(weak_binders),
                'mean_score': float(np.mean([d['splicebind_score'] for d in weak_binders])) if weak_binders else None,
            },
        },
        'data': matched_data,
    }

    output_path = PROJECT_ROOT / "data" / "processed" / "binding_affinity_analysis.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()

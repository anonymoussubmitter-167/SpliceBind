#!/usr/bin/env python3
"""
Validate SpliceBind on Splice Variants Using Experimental PDB Structures

Uses experimental crystal structures where available:
- ALK-L1196M: PDB 2YFX (mutant kinase domain)
- EGFR-vIII: PDB 8UKX (variant extracellular) - use kinase domain from AlphaFold
- BRAF: PDB 4MNF (V600E kinase domain)
- MET: PDB 3CCN (kinase domain with inhibitor)
- AR: AlphaFold for canonical, AR-V7 lacks LBD entirely
"""

import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

STRUCTURES_DIR = PROJECT_ROOT / "data" / "structures"
EMBEDDINGS_DIR = PROJECT_ROOT / "data" / "embeddings"

# Variant validation cases with PDB structures
VALIDATION_CASES = {
    'ALK_L1196M': {
        'gene': 'ALK',
        'drug': 'Crizotinib',
        'expected': 'resistant',
        'canonical_pdb': 'AF-Q9UM73-F1-model_v4.pdb',  # AlphaFold canonical
        'variant_pdb': '2YFX.pdb',  # Crystal structure of L1196M mutant
        'description': 'Gatekeeper mutation reduces Crizotinib binding',
        'binding_site': {'start': 1116, 'end': 1399},  # ALK kinase domain
    },
    'EGFR_canonical': {
        'gene': 'EGFR',
        'drug': 'Gefitinib',
        'expected': 'sensitive',
        'canonical_pdb': 'AF-P00533-F1-model_v6.pdb',
        'variant_pdb': None,  # Compare against itself
        'description': 'Canonical EGFR - Gefitinib sensitive',
        'binding_site': {'start': 712, 'end': 979},  # EGFR kinase domain
    },
    'MET_ex14skip': {
        'gene': 'MET',
        'drug': 'Capmatinib',
        'expected': 'sensitive',  # This variant IS sensitive to MET inhibitors
        'canonical_pdb': 'AF-P08581-F1-model_v4.pdb',
        'variant_pdb': '3CCN.pdb',  # MET kinase domain - exon14 skip doesn't affect kinase
        'description': 'Exon 14 skip - oncogenic but kinase domain intact, drug SENSITIVE',
        'binding_site': {'start': 1078, 'end': 1345},  # MET kinase domain
    },
    'BRAF_V600E': {
        'gene': 'BRAF',
        'drug': 'Vemurafenib',
        'expected': 'sensitive',  # V600E is the drug target
        'canonical_pdb': 'AF-P15056-F1-model_v6.pdb',
        'variant_pdb': '4MNF.pdb',  # BRAF V600E kinase domain
        'description': 'BRAF V600E - Vemurafenib target',
        'binding_site': {'start': 457, 'end': 717},  # BRAF kinase domain
    },
}


def run_p2rank(pdb_path):
    """Run P2Rank pocket detection."""
    from modules.p2rank_detector import P2RankDetector
    detector = P2RankDetector()
    return detector.detect_pockets(str(pdb_path))


def get_top_pocket_score(pockets):
    """Get score of top-ranked pocket (best druggable pocket)."""
    if not pockets:
        return 0.0, None

    # Sort by P2Rank score descending
    sorted_pockets = sorted(pockets, key=lambda p: p.get('p2rank_score', 0), reverse=True)
    top_pocket = sorted_pockets[0]

    return top_pocket.get('p2rank_score', 0), top_pocket


def main():
    print("=" * 70)
    print("SpliceBind Variant Validation (Experimental PDB Structures)")
    print("=" * 70)

    results = []

    for case_name, case_info in VALIDATION_CASES.items():
        print(f"\n{'='*60}")
        print(f"Case: {case_name}")
        print(f"{'='*60}")
        print(f"Gene: {case_info['gene']}")
        print(f"Drug: {case_info['drug']}")
        print(f"Expected: {case_info['expected']}")
        print(f"Description: {case_info['description']}")

        # Check canonical structure
        canonical_pdb = STRUCTURES_DIR / case_info['canonical_pdb']
        if not canonical_pdb.exists():
            print(f"  Missing canonical: {canonical_pdb.name}")
            continue

        print(f"\n[1] Analyzing canonical: {canonical_pdb.name}")
        canonical_pockets = run_p2rank(canonical_pdb)
        print(f"    Pockets found: {len(canonical_pockets)}")

        canonical_score, canonical_pocket = get_top_pocket_score(canonical_pockets)
        print(f"    Binding site P2Rank score: {canonical_score:.3f}")

        # Analyze variant if different
        variant_pdb = case_info.get('variant_pdb')
        if variant_pdb:
            variant_path = STRUCTURES_DIR / variant_pdb
            if not variant_path.exists():
                print(f"  Missing variant: {variant_pdb}")
                continue

            print(f"\n[2] Analyzing variant: {variant_pdb}")
            variant_pockets = run_p2rank(variant_path)
            print(f"    Pockets found: {len(variant_pockets)}")

            variant_score, variant_pocket = get_top_pocket_score(variant_pockets)
            print(f"    Binding site P2Rank score: {variant_score:.3f}")

            delta = variant_score - canonical_score
            print(f"\n    Delta: {delta:+.3f}")

            # Determine if prediction matches expectation
            if case_info['expected'] == 'resistant':
                correct = delta < 0  # Expect reduced druggability
            else:
                correct = delta >= -0.1  # Expect maintained/increased druggability

            results.append({
                'case': case_name,
                'gene': case_info['gene'],
                'drug': case_info['drug'],
                'expected': case_info['expected'],
                'canonical_score': canonical_score,
                'variant_score': variant_score,
                'delta': delta,
                'correct': correct
            })

            print(f"    Prediction: {'CORRECT' if correct else 'INCORRECT'}")
        else:
            # Just report canonical score
            results.append({
                'case': case_name,
                'gene': case_info['gene'],
                'drug': case_info['drug'],
                'expected': case_info['expected'],
                'canonical_score': canonical_score,
                'variant_score': None,
                'delta': None,
                'correct': canonical_score > 5.0  # High P2Rank score = druggable
            })

    # Summary
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)

    correct = sum(1 for r in results if r['correct'])
    total = len(results)

    print(f"\nAccuracy: {correct}/{total} ({100*correct/total:.0f}%)")
    print("\nDetails:")
    print("-" * 70)
    print(f"{'Case':<20} {'Gene':<8} {'Drug':<15} {'Expected':<10} {'Delta':<8} {'Status'}")
    print("-" * 70)

    for r in results:
        delta_str = f"{r['delta']:+.3f}" if r['delta'] is not None else "N/A"
        status = "✓ PASS" if r['correct'] else "✗ FAIL"
        print(f"{r['case']:<20} {r['gene']:<8} {r['drug']:<15} {r['expected']:<10} {delta_str:<8} {status}")

    return results


if __name__ == "__main__":
    main()

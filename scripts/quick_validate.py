#!/usr/bin/env python3
"""
Quick validation of PIK3CD-S isoform using existing data.
"""

import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

STRUCTURES_DIR = PROJECT_ROOT / "data" / "structures"
EMBEDDINGS_DIR = PROJECT_ROOT / "data" / "embeddings"

def run_p2rank(pdb_path):
    """Run P2Rank pocket detection."""
    from modules.p2rank_detector import P2RankDetector
    detector = P2RankDetector()
    return detector.detect_pockets(str(pdb_path))


def extract_features(pdb_path):
    """Extract pocket features."""
    from modules.pocket_analyzer import PocketFeatureExtractor
    extractor = PocketFeatureExtractor()
    return extractor.extract_features(str(pdb_path))


def main():
    print("=" * 60)
    print("PIK3CD Isoform Validation")
    print("=" * 60)

    # Canonical structure (AlphaFold)
    canonical_pdb = STRUCTURES_DIR / "AF-O00329-F1-model_v6.pdb"
    variant_pdb = STRUCTURES_DIR / "PIK3CD-S_esmfold.pdb"

    canonical_emb = EMBEDDINGS_DIR / "PIK3CD_canonical_esm2.npy"
    variant_emb = EMBEDDINGS_DIR / "PIK3CD_variant_esm2.npy"

    print(f"\nCanonical: {canonical_pdb.name}")
    print(f"Variant:   {variant_pdb.name}")

    # Check files exist
    if not canonical_pdb.exists():
        print(f"Missing: {canonical_pdb}")
        return
    if not variant_pdb.exists():
        print(f"Missing: {variant_pdb}")
        return

    # Run P2Rank on both
    print("\n[1] Running P2Rank pocket detection...")

    canonical_pockets = run_p2rank(canonical_pdb)
    variant_pockets = run_p2rank(variant_pdb)

    print(f"  Canonical pockets: {len(canonical_pockets)}")
    print(f"  Variant pockets:   {len(variant_pockets)}")

    # Get top pocket scores
    if canonical_pockets:
        canonical_top = max(p.get('p2rank_score', 0) for p in canonical_pockets)
        print(f"  Canonical top P2Rank score: {canonical_top:.3f}")
    else:
        canonical_top = 0

    if variant_pockets:
        variant_top = max(p.get('p2rank_score', 0) for p in variant_pockets)
        print(f"  Variant top P2Rank score:   {variant_top:.3f}")
    else:
        variant_top = 0

    # Compare P2Rank scores
    p2rank_delta = variant_top - canonical_top
    print(f"\n  P2Rank Delta: {p2rank_delta:+.3f}")

    # Load ESM embeddings
    print("\n[2] ESM-2 embedding analysis...")

    if canonical_emb.exists() and variant_emb.exists():
        can_emb = np.load(canonical_emb)
        var_emb = np.load(variant_emb)

        print(f"  Canonical embeddings shape: {can_emb.shape}")
        print(f"  Variant embeddings shape:   {var_emb.shape}")

        # Simple comparison - mean embedding distance
        min_len = min(len(can_emb), len(var_emb))
        emb_distance = np.mean(np.linalg.norm(can_emb[:min_len] - var_emb[:min_len], axis=1))
        print(f"  Mean embedding distance: {emb_distance:.3f}")

    # PIK3CD binding site analysis
    print("\n[3] Binding site analysis...")

    # PIK3CD ATP binding site: 695-1044
    binding_site = {'start': 695, 'end': 1044}

    # Check pocket overlap with binding site
    def check_binding_overlap(pockets, binding_site):
        best_overlap = 0
        best_pocket = None
        for i, pocket in enumerate(pockets):
            residues = set(pocket.get('residue_ids', []))
            site = set(range(binding_site['start'], binding_site['end'] + 1))
            overlap = len(residues & site) / max(len(residues), 1)
            if overlap > best_overlap:
                best_overlap = overlap
                best_pocket = pocket
        return best_pocket, best_overlap

    can_pocket, can_overlap = check_binding_overlap(canonical_pockets, binding_site)
    var_pocket, var_overlap = check_binding_overlap(variant_pockets, binding_site)

    print(f"  Canonical binding site overlap: {can_overlap:.1%}")
    print(f"  Variant binding site overlap:   {var_overlap:.1%}")

    if can_pocket and var_pocket:
        can_score = can_pocket.get('score', 0)
        var_score = var_pocket.get('score', 0)
        print(f"\n  Binding pocket scores:")
        print(f"    Canonical: {can_score:.3f}")
        print(f"    Variant:   {var_score:.3f}")
        print(f"    Delta:     {var_score - can_score:+.3f}")

    # Summary
    print("\n" + "=" * 60)
    print("VALIDATION RESULT")
    print("=" * 60)
    print(f"""
Gene:     PIK3CD
Drug:     Idelalisib
Variant:  PIK3CD-S (spliced short isoform)
Expected: Resistant (reduced druggability)

P2Rank Delta: {p2rank_delta:+.3f}
Prediction:   {'REDUCED' if p2rank_delta < 0 else 'INCREASED'} druggability
Expected:     REDUCED druggability
Status:       {'✓ CORRECT' if p2rank_delta < 0 else '✗ INCORRECT'}
""")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Validate SpliceBind on Clinically Relevant Splice Variants

This script generates structures for splice variants associated with drug resistance
and compares druggability predictions between canonical and variant forms.

Target Variants:
1. PIK3CD-S: Idelalisib resistance (DONE)
2. p61-BRAF: Vemurafenib resistance (exon 4-8 deletion)
3. EGFR-vIII: Gefitinib resistance (exon 2-7 deletion)
4. MET-ex14skip: Capmatinib SENSITIVE (exon 14 deletion)
5. AR-V7: Enzalutamide resistance (truncated, no LBD)
6. ALK-L1196M: Crizotinib resistance (gatekeeper mutation)
"""

import os
import sys
import json
import requests
import subprocess
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Optional

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

STRUCTURES_DIR = PROJECT_ROOT / "data" / "structures"
EMBEDDINGS_DIR = PROJECT_ROOT / "data" / "embeddings"
RESULTS_DIR = PROJECT_ROOT / "data" / "processed"

# =============================================================================
# VARIANT DEFINITIONS
# =============================================================================

VARIANTS = {
    'BRAF_p61': {
        'gene': 'BRAF',
        'uniprot': 'P15056',
        'drug': 'Vemurafenib',
        'expected': 'resistant',
        'description': 'Exon 4-8 deletion (aa 1-187 deleted), lacks RAS-binding domain',
        'modification': 'delete_region',
        'delete_start': 1,
        'delete_end': 187,
    },
    'EGFR_vIII': {
        'gene': 'EGFR',
        'uniprot': 'P00533',
        'drug': 'Gefitinib',
        'expected': 'resistant',
        'description': 'Exon 2-7 deletion (aa 6-273 deleted), constitutively active',
        'modification': 'delete_region',
        'delete_start': 6,
        'delete_end': 273,
    },
    'MET_ex14skip': {
        'gene': 'MET',
        'uniprot': 'P08581',
        'drug': 'Capmatinib',
        'expected': 'sensitive',  # This variant is SENSITIVE to MET inhibitors!
        'description': 'Exon 14 skipping (aa 964-1010 deleted), oncogenic driver',
        'modification': 'delete_region',
        'delete_start': 964,
        'delete_end': 1010,
    },
    'AR_V7': {
        'gene': 'AR',
        'uniprot': 'P10275',
        'drug': 'Enzalutamide',
        'expected': 'resistant',
        'description': 'Cryptic exon inclusion, truncated at aa 644 (no LBD)',
        'modification': 'truncate',
        'truncate_at': 644,
    },
    'ALK_L1196M': {
        'gene': 'ALK',
        'uniprot': 'Q9UM73',
        'drug': 'Crizotinib',
        'expected': 'resistant',
        'description': 'Gatekeeper mutation L1196M',
        'modification': 'point_mutation',
        'position': 1196,
        'from_aa': 'L',
        'to_aa': 'M',
    },
}


def get_uniprot_sequence(uniprot_id: str) -> str:
    """Fetch canonical sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            sequence = ''.join(lines[1:])  # Skip header
            return sequence
    except Exception as e:
        print(f"Error fetching {uniprot_id}: {e}")
    return ""


def create_variant_sequence(canonical_seq: str, variant_info: Dict) -> str:
    """Create variant sequence based on modification type."""
    mod_type = variant_info['modification']

    if mod_type == 'delete_region':
        start = variant_info['delete_start'] - 1  # 0-indexed
        end = variant_info['delete_end']
        variant_seq = canonical_seq[:start] + canonical_seq[end:]

    elif mod_type == 'truncate':
        truncate_at = variant_info['truncate_at']
        variant_seq = canonical_seq[:truncate_at]

    elif mod_type == 'point_mutation':
        pos = variant_info['position'] - 1  # 0-indexed
        from_aa = variant_info['from_aa']
        to_aa = variant_info['to_aa']

        if canonical_seq[pos] != from_aa:
            print(f"Warning: Expected {from_aa} at position {pos+1}, found {canonical_seq[pos]}")

        variant_seq = canonical_seq[:pos] + to_aa + canonical_seq[pos+1:]
    else:
        variant_seq = canonical_seq

    return variant_seq


def generate_esmfold_structure(sequence: str, output_path: Path, name: str) -> bool:
    """Generate structure using ESMFold API."""
    if output_path.exists():
        print(f"  Structure exists: {output_path.name}")
        return True

    print(f"  Generating ESMFold structure for {name}...")

    # ESMFold API
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    try:
        response = requests.post(
            url,
            data=sequence,
            headers={'Content-Type': 'text/plain'},
            timeout=300
        )

        if response.status_code == 200:
            output_path.write_text(response.text)
            print(f"  Generated: {output_path.name}")
            return True
        else:
            print(f"  ESMFold API error: {response.status_code}")
            return False

    except Exception as e:
        print(f"  ESMFold error: {e}")
        return False


def generate_esm_embeddings(sequence: str, output_path: Path, name: str) -> bool:
    """Generate ESM-2 embeddings for a sequence."""
    if output_path.exists():
        print(f"  Embeddings exist: {output_path.name}")
        return True

    print(f"  Generating ESM-2 embeddings for {name}...")

    try:
        import torch
        import esm

        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        batch_converter = alphabet.get_batch_converter()

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model = model.to(device)
        model.eval()

        data = [(name, sequence)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33])

        embeddings = results['representations'][33][0, 1:len(sequence)+1, :].cpu().numpy()

        # Reduce to 32 dimensions
        if embeddings.shape[1] > 32:
            n_cols = embeddings.shape[1]
            step = n_cols // 32
            reduced = np.zeros((embeddings.shape[0], 32))
            for i in range(32):
                start = i * step
                end = start + step if i < 31 else n_cols
                reduced[:, i] = embeddings[:, start:end].mean(axis=1)
            embeddings = reduced

        np.save(output_path, embeddings)
        print(f"  Generated: {output_path.name}")
        return True

    except Exception as e:
        print(f"  ESM error: {e}")
        return False


def run_prediction(pdb_path: Path, esm_path: Path, gene: str) -> Optional[float]:
    """Run SpliceBind prediction on a structure."""
    # Import the prediction components from train.py
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

    try:
        from train import (
            P2RankDetector, PocketFeatureExtractor, PocketGNNv5,
            build_pocket_graph, KINASE_ATP_SITES
        )
        import torch
        from torch_geometric.data import Batch

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

        # Load ESM embeddings
        esm_embeddings = {}
        if esm_path.exists():
            esm_embeddings[pdb_path.stem] = np.load(esm_path)

        # Initialize components
        p2rank = P2RankDetector()
        feat_extractor = PocketFeatureExtractor()

        # Detect pockets
        pockets = p2rank.detect_pockets(str(pdb_path))
        if not pockets:
            print(f"  No pockets detected")
            return None

        # Extract features
        features = feat_extractor.extract_features(str(pdb_path))
        ca_coords = features.get('ca_coords')
        if ca_coords is None or len(ca_coords) == 0:
            return None

        if isinstance(ca_coords, list):
            ca_coords = np.array(ca_coords)

        # Get binding site info
        binding_site = KINASE_ATP_SITES.get(gene, None)

        # Build graphs for each pocket
        graphs = []
        pocket_scores = []

        for pocket in pockets:
            pocket_residues = pocket.get('residue_ids', [])

            # Get ESM embeddings for this structure
            esm_emb = esm_embeddings.get(pdb_path.stem)

            graph = build_pocket_graph(
                pocket_residues=pocket_residues,
                ca_coords=ca_coords,
                features=features,
                esm_embeddings=esm_emb,
                pocket_info=pocket
            )

            if graph is not None and graph.x.shape[1] == 56:
                graphs.append(graph)
                pocket_scores.append(pocket.get('score', 0))

        if not graphs:
            return None

        # Load trained model (use the saved weights if available)
        model = PocketGNNv5(node_dim=56, hidden_dim=256, num_layers=3)

        # Try to load saved weights
        model_path = PROJECT_ROOT / "models" / "splicebind_best.pt"
        if model_path.exists():
            model.load_state_dict(torch.load(model_path, map_location=device))

        model = model.to(device)
        model.eval()

        # Run predictions
        batch = Batch.from_data_list(graphs).to(device)

        with torch.no_grad():
            logits = model(batch)
            probs = torch.sigmoid(logits).cpu().numpy()

        # Return max druggability score
        max_score = float(np.max(probs))

        # Also check if top pocket overlaps binding site
        if binding_site:
            for i, (graph, prob) in enumerate(zip(graphs, probs)):
                pocket_residues = set(pockets[i].get('residue_ids', []))
                site_residues = set(range(binding_site['start'], binding_site['end'] + 1))
                overlap = len(pocket_residues & site_residues) / max(len(pocket_residues), 1)
                if overlap > 0.3:
                    return float(prob)

        return max_score

    except Exception as e:
        print(f"  Prediction error: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    """Main validation pipeline."""
    print("=" * 70)
    print("SpliceBind Isoform Validation")
    print("=" * 70)

    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    EMBEDDINGS_DIR.mkdir(parents=True, exist_ok=True)

    results = []

    for variant_name, variant_info in VARIANTS.items():
        print(f"\n{'='*60}")
        print(f"Processing {variant_name}")
        print(f"{'='*60}")
        print(f"Gene: {variant_info['gene']}")
        print(f"Drug: {variant_info['drug']}")
        print(f"Expected: {variant_info['expected']}")
        print(f"Description: {variant_info['description']}")

        # Get canonical sequence
        print(f"\n[1] Fetching canonical sequence...")
        canonical_seq = get_uniprot_sequence(variant_info['uniprot'])
        if not canonical_seq:
            print(f"  Failed to get sequence for {variant_info['uniprot']}")
            continue
        print(f"  Canonical length: {len(canonical_seq)} aa")

        # Create variant sequence
        print(f"\n[2] Creating variant sequence...")
        variant_seq = create_variant_sequence(canonical_seq, variant_info)
        print(f"  Variant length: {len(variant_seq)} aa")
        print(f"  Difference: {len(canonical_seq) - len(variant_seq)} aa")

        # Generate structures
        print(f"\n[3] Generating structures...")

        canonical_pdb = STRUCTURES_DIR / f"{variant_info['gene']}_canonical.pdb"
        variant_pdb = STRUCTURES_DIR / f"{variant_name}.pdb"

        # Check if we already have AlphaFold structure for canonical
        af_pattern = f"AF-{variant_info['uniprot']}-F1-model_v*.pdb"
        af_files = list(STRUCTURES_DIR.glob(af_pattern))

        if af_files:
            canonical_pdb = af_files[0]
            print(f"  Using AlphaFold structure: {canonical_pdb.name}")
        else:
            generate_esmfold_structure(canonical_seq, canonical_pdb, f"{variant_info['gene']}_canonical")

        generate_esmfold_structure(variant_seq, variant_pdb, variant_name)

        # Generate ESM embeddings
        print(f"\n[4] Generating ESM-2 embeddings...")

        canonical_emb = EMBEDDINGS_DIR / f"{canonical_pdb.stem}_esm2.npy"
        variant_emb = EMBEDDINGS_DIR / f"{variant_name}_esm2.npy"

        # Check for existing embeddings
        if not canonical_emb.exists():
            # Try to find matching embedding
            for emb_file in EMBEDDINGS_DIR.glob("*.npy"):
                if variant_info['uniprot'] in emb_file.stem:
                    canonical_emb = emb_file
                    break

        if not canonical_emb.exists():
            generate_esm_embeddings(canonical_seq, canonical_emb, f"{variant_info['gene']}_canonical")

        generate_esm_embeddings(variant_seq, variant_emb, variant_name)

        # Run predictions
        print(f"\n[5] Running druggability predictions...")

        canonical_score = run_prediction(canonical_pdb, canonical_emb, variant_info['gene'])
        variant_score = run_prediction(variant_pdb, variant_emb, variant_info['gene'])

        if canonical_score is not None and variant_score is not None:
            delta = variant_score - canonical_score

            # Determine prediction
            if variant_info['expected'] == 'resistant':
                # Expect variant to have LOWER druggability
                prediction_correct = delta < 0
            else:  # sensitive
                # Expect variant to have SIMILAR or HIGHER druggability
                prediction_correct = delta >= -0.05

            result = {
                'variant': variant_name,
                'gene': variant_info['gene'],
                'drug': variant_info['drug'],
                'expected': variant_info['expected'],
                'canonical_score': canonical_score,
                'variant_score': variant_score,
                'delta': delta,
                'prediction_correct': prediction_correct,
            }
            results.append(result)

            print(f"\n  Results:")
            print(f"    Canonical score: {canonical_score:.3f}")
            print(f"    Variant score:   {variant_score:.3f}")
            print(f"    Delta:           {delta:+.3f}")
            print(f"    Expected:        {variant_info['expected']}")
            print(f"    Prediction:      {'CORRECT' if prediction_correct else 'INCORRECT'}")
        else:
            print(f"  Could not get predictions")

    # Summary
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)

    if results:
        correct = sum(1 for r in results if r['prediction_correct'])
        total = len(results)

        print(f"\nAccuracy: {correct}/{total} ({100*correct/total:.0f}%)")
        print("\nDetails:")
        print("-" * 70)
        print(f"{'Variant':<15} {'Gene':<8} {'Drug':<15} {'Expected':<10} {'Δ':<8} {'Status'}")
        print("-" * 70)

        for r in results:
            status = "✓ PASS" if r['prediction_correct'] else "✗ FAIL"
            print(f"{r['variant']:<15} {r['gene']:<8} {r['drug']:<15} {r['expected']:<10} {r['delta']:+.3f}   {status}")

        # Save results
        output_path = RESULTS_DIR / "isoform_validation.json"
        with open(output_path, 'w') as f:
            json.dump({
                'accuracy': correct / total,
                'correct': correct,
                'total': total,
                'results': results
            }, f, indent=2)
        print(f"\nResults saved to {output_path}")
    else:
        print("No results to report")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
SpliceBind Variant Predictions: Run inference on canonical vs variant structures.

This script:
1. Trains a SpliceBind model on all data
2. Runs inference on splice variant structures
3. Compares canonical vs variant druggability scores
"""

import sys
import json
import numpy as np
import torch
import torch.nn as nn
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict
from torch_geometric.data import Data, Batch
from torch_geometric.nn import EdgeConv, global_mean_pool, global_max_pool
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

STRUCTURES_DIR = PROJECT_ROOT / "data" / "structures"
EMBEDDINGS_DIR = PROJECT_ROOT / "data" / "embeddings"


# =============================================================================
# MODEL
# =============================================================================

class PocketGNNv5(nn.Module):
    def __init__(self, node_dim=56, hidden_dim=256, num_layers=3, dropout=0.3):
        super().__init__()
        self.dropout = nn.Dropout(dropout)
        self.convs = nn.ModuleList()
        in_dim = node_dim
        for _ in range(num_layers):
            mlp = nn.Sequential(
                nn.Linear(in_dim * 2, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim)
            )
            self.convs.append(EdgeConv(mlp, aggr='max'))
            in_dim = hidden_dim

        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 1)
        )

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        for conv in self.convs:
            x = conv(x, edge_index)
            x = self.dropout(x)
        x_mean = global_mean_pool(x, batch)
        x_max = global_max_pool(x, batch)
        x = torch.cat([x_mean, x_max], dim=1)
        return self.classifier(x).squeeze(-1)


# =============================================================================
# SPLICE VARIANTS
# =============================================================================

SPLICE_VARIANTS = [
    {
        'name': 'AR-V7',
        'gene': 'AR',
        'canonical': 'AR_kinase_canonical.pdb',
        'variant': None,  # LBD absent
        'drug': 'Enzalutamide',
        'clinical': 'Resistant',
        'mechanism': 'LBD deletion - drug binding domain absent',
    },
    {
        'name': 'PIK3CD-S',
        'gene': 'PIK3CD',
        'canonical': 'AF-O00329-F1-model_v6.pdb',
        'variant': None,  # Need to generate
        'drug': 'Idelalisib',
        'clinical': 'Resistant',
        'mechanism': 'Alternative splicing affects kinase domain',
    },
    {
        'name': 'MET-ex14skip',
        'gene': 'MET',
        'canonical': 'AF-P08581-F1-model_v4.pdb',
        'variant': 'MET_ex14skip_kinase.pdb',
        'drug': 'Capmatinib',
        'clinical': 'Sensitive',
        'mechanism': 'Exon 14 skip - kinase preserved',
    },
    {
        'name': 'BRAF-p61',
        'gene': 'BRAF',
        'canonical': 'AF-P15056-F1-model_v6.pdb',
        'variant': 'BRAF_p61_kinase.pdb',
        'drug': 'Vemurafenib',
        'clinical': 'Resistant',
        'mechanism': 'Truncation - kinase identical',
    },
    {
        'name': 'ALK-L1196M',
        'gene': 'ALK',
        'canonical': 'AF-Q9UM73-F1-model_v4.pdb',
        'variant': 'ALK_L1196M_kinase.pdb',
        'drug': 'Crizotinib',
        'clinical': 'Resistant',
        'mechanism': 'Gatekeeper mutation',
    },
]


def load_esm_embedding(structure_name: str) -> Optional[np.ndarray]:
    """Load ESM embedding for a structure."""
    # Try different naming conventions
    candidates = [
        EMBEDDINGS_DIR / f"{structure_name.replace('.pdb', '')}_esm.npy",
        EMBEDDINGS_DIR / f"{structure_name.replace('.pdb', '')}.npy",
    ]

    # Also try by UniProt ID
    if 'AF-' in structure_name:
        uniprot = structure_name.split('-')[1]
        candidates.append(EMBEDDINGS_DIR / f"{uniprot}_esm.npy")
        candidates.append(EMBEDDINGS_DIR / f"AF-{uniprot}-F1-model_v6_esm.npy")
        candidates.append(EMBEDDINGS_DIR / f"AF-{uniprot}-F1-model_v4_esm.npy")

    for path in candidates:
        if path.exists():
            return np.load(path)

    return None


def parse_pdb_coords(pdb_path: Path) -> Dict[int, np.ndarray]:
    """Parse CA coordinates from PDB file by residue index."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                try:
                    res_idx = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords[res_idx] = np.array([x, y, z])
                except:
                    continue
    return coords


def process_structure(pdb_path: Path) -> List[Dict]:
    """Run P2Rank and build pocket features."""
    from modules.p2rank_detector import P2RankDetector
    from modules.pocket_analyzer import PocketGraphBuilder

    detector = P2RankDetector()
    builder = PocketGraphBuilder()

    # Detect pockets
    pockets = detector.detect_pockets(str(pdb_path))
    if not pockets:
        return []

    # Parse coordinates from PDB
    all_coords = parse_pdb_coords(pdb_path)

    # Build graphs using simplified features (no ESM for variant comparison)
    graphs = []
    for pocket in pockets:
        try:
            residue_indices = pocket.get('residue_indices', [])
            residue_names = pocket.get('residue_names', [])

            if not residue_indices:
                continue

            # Get coordinates for pocket residues
            coords = []
            valid_indices = []
            valid_names = []
            for idx, name in zip(residue_indices, residue_names):
                if idx in all_coords:
                    coords.append(all_coords[idx])
                    valid_indices.append(idx)
                    valid_names.append(name)

            if len(coords) < 3:  # Need at least 3 residues
                continue

            coords = np.array(coords)

            # Build simple node features (one-hot AA + basic features)
            # Model expects 24 features: 20 one-hot + 3 physico + 1 extra
            node_features = []
            for res_name in valid_names:
                feat = builder._residue_features(res_name)
                # Pad to 24 features if needed
                while len(feat) < 24:
                    feat.append(0.0)
                node_features.append(feat[:24])

            node_features = np.array(node_features, dtype=np.float32)

            # Build edges based on distance
            edge_index = []
            n = len(coords)
            for i in range(n):
                for j in range(i + 1, n):
                    dist = np.linalg.norm(coords[i] - coords[j])
                    if dist < 8.0:  # 8 Angstrom cutoff
                        edge_index.append([i, j])
                        edge_index.append([j, i])

            if not edge_index:
                # Fully connected if no edges
                for i in range(n):
                    for j in range(n):
                        if i != j:
                            edge_index.append([i, j])

            edge_index = np.array(edge_index, dtype=np.int64).T

            graphs.append({
                'node_features': node_features,
                'edge_index': edge_index,
                'p2rank_score': pocket.get('p2rank_score', 0),
                'pocket_id': pocket.get('name', f'pocket_{len(graphs)}'),
            })

        except Exception as e:
            print(f"    Error processing pocket: {e}")
            continue

    return graphs


def predict_druggability(model: nn.Module, graphs: List[Dict], device: torch.device) -> List[float]:
    """Run model inference on pocket graphs."""
    model.eval()
    scores = []

    with torch.no_grad():
        for graph in graphs:
            x = torch.tensor(graph['node_features'], dtype=torch.float32).to(device)
            edge_index = torch.tensor(graph['edge_index'], dtype=torch.long).to(device)
            batch = torch.zeros(x.size(0), dtype=torch.long).to(device)

            data = Data(x=x, edge_index=edge_index, batch=batch)
            logit = model(data)
            prob = torch.sigmoid(logit).item()
            scores.append(prob)

    return scores


def train_full_model():
    """Train model on all data and return trained model."""
    print("Training SpliceBind model on full dataset...")

    # Load pre-processed data if available
    data_path = PROJECT_ROOT / "data" / "processed" / "pocket_data.pt"

    if data_path.exists():
        print(f"Loading cached data from {data_path}")
        data = torch.load(data_path, weights_only=False)

        # Extract graphs and labels from list of pocket dicts
        graphs = []
        labels = []
        for pocket in data:
            if 'pocket_graph' in pocket and pocket['pocket_graph'] is not None:
                graphs.append(pocket['pocket_graph'])
                # Use druggability_score as label (threshold at 0.5)
                score = pocket.get('druggability_score', 0)
                labels.append(1 if score > 0.5 else 0)

        print(f"Loaded {len(graphs)} pockets")
    else:
        print("No cached data. Run train.py first to generate pocket_data.pt")
        return None

    # Create dataset
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Determine node dimension from data
    node_dim = graphs[0]['node_features'].shape[1] if hasattr(graphs[0]['node_features'], 'shape') else len(graphs[0]['node_features'][0])
    print(f"Node dimension: {node_dim}")

    # Train model
    model = PocketGNNv5(node_dim=node_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([5.0]).to(device))

    # Convert to tensors
    dataset = []
    for g, label in zip(graphs, labels):
        x = torch.tensor(g['node_features'], dtype=torch.float32)
        edge_index = torch.tensor(g['edge_index'], dtype=torch.long)
        y = torch.tensor([label], dtype=torch.float32)
        dataset.append(Data(x=x, edge_index=edge_index, y=y))

    # Simple training loop
    model.train()
    for epoch in range(30):
        total_loss = 0
        np.random.shuffle(dataset)

        for data in dataset:
            data = data.to(device)
            optimizer.zero_grad()
            out = model(data)
            loss = criterion(out, data.y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        if (epoch + 1) % 10 == 0:
            print(f"  Epoch {epoch+1}: Loss = {total_loss/len(dataset):.4f}")

    return model


def main():
    print("=" * 70)
    print("SPLICEBIND VARIANT PREDICTIONS")
    print("=" * 70)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Try to load pre-trained model or train new one
    model_path = PROJECT_ROOT / "data" / "processed" / "splicebind_model.pt"

    if model_path.exists():
        print(f"Loading model from {model_path}")
        checkpoint = torch.load(model_path, map_location=device, weights_only=False)
        node_dim = checkpoint.get('node_dim', 24)
        model = PocketGNNv5(node_dim=node_dim).to(device)
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
        model = train_full_model()
        if model is None:
            print("Could not train model. Exiting.")
            return

        # Save model with metadata
        torch.save({
            'model_state_dict': model.state_dict(),
            'node_dim': 24,  # From cached data
        }, model_path)
        print(f"Model saved to {model_path}")

    model.eval()

    # Process variants
    results = []

    print("\n" + "=" * 70)
    print("VARIANT ANALYSIS")
    print("=" * 70)

    for variant in SPLICE_VARIANTS:
        print(f"\n{'-'*60}")
        print(f"{variant['name']}")
        print(f"{'-'*60}")
        print(f"Gene: {variant['gene']}")
        print(f"Drug: {variant['drug']}")
        print(f"Clinical outcome: {variant['clinical']}")
        print(f"Mechanism: {variant['mechanism']}")

        # Process canonical
        canonical_path = STRUCTURES_DIR / variant['canonical']
        canonical_score = None
        canonical_p2rank = None

        if canonical_path.exists():
            print(f"\nCanonical: {variant['canonical']}")
            graphs = process_structure(canonical_path)
            if graphs:
                scores = predict_druggability(model, graphs, device)
                p2rank_scores = [g['p2rank_score'] for g in graphs]
                canonical_score = max(scores)
                canonical_p2rank = max(p2rank_scores)
                print(f"  SpliceBind: {canonical_score:.3f}")
                print(f"  P2Rank: {canonical_p2rank:.2f}")
                print(f"  Pockets: {len(graphs)}")
        else:
            print(f"  Canonical structure not found: {variant['canonical']}")

        # Process variant
        variant_score = None
        variant_p2rank = None

        if variant['variant']:
            variant_path = STRUCTURES_DIR / variant['variant']
            if variant_path.exists():
                print(f"\nVariant: {variant['variant']}")
                graphs = process_structure(variant_path)
                if graphs:
                    scores = predict_druggability(model, graphs, device)
                    p2rank_scores = [g['p2rank_score'] for g in graphs]
                    variant_score = max(scores)
                    variant_p2rank = max(p2rank_scores)
                    print(f"  SpliceBind: {variant_score:.3f}")
                    print(f"  P2Rank: {variant_p2rank:.2f}")
                    print(f"  Pockets: {len(graphs)}")
            else:
                print(f"  Variant structure not found: {variant['variant']}")
        else:
            # AR-V7 case: no binding domain in variant
            print(f"\nVariant: No binding domain (LBD absent)")
            variant_score = 0.0
            variant_p2rank = 0.0

        # Compute delta
        if canonical_score is not None and variant_score is not None:
            delta_sb = variant_score - canonical_score
            delta_p2rank = variant_p2rank - canonical_p2rank if variant_p2rank else 0

            print(f"\nDelta (Variant - Canonical):")
            print(f"  SpliceBind: {delta_sb:+.3f}")
            print(f"  P2Rank: {delta_p2rank:+.2f}")

            # Interpret
            if variant['clinical'] == 'Resistant':
                if delta_sb < -0.05:
                    interpretation = "CORRECT: Reduced druggability (resistant)"
                elif abs(delta_sb) < 0.05:
                    interpretation = "PARTIAL: Similar druggability (structural mechanism)"
                else:
                    interpretation = "UNEXPECTED: Increased druggability"
            else:  # Sensitive
                if abs(delta_sb) < 0.1:
                    interpretation = "CORRECT: Maintained druggability (sensitive)"
                else:
                    interpretation = "CHECK: Changed druggability for sensitive variant"

            print(f"  Interpretation: {interpretation}")

            results.append({
                'variant': variant['name'],
                'gene': variant['gene'],
                'drug': variant['drug'],
                'clinical': variant['clinical'],
                'mechanism': variant['mechanism'],
                'canonical_sb': float(canonical_score),
                'variant_sb': float(variant_score) if variant_score else None,
                'delta_sb': float(delta_sb),
                'canonical_p2rank': float(canonical_p2rank) if canonical_p2rank else None,
                'variant_p2rank': float(variant_p2rank) if variant_p2rank else None,
                'interpretation': interpretation,
            })

    # Summary table
    print("\n" + "=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)

    print(f"\n{'Variant':<15} {'Drug':<12} {'Clinical':<10} {'Can SB':<8} {'Var SB':<8} {'Delta':<8} {'Result'}")
    print("-" * 85)

    for r in results:
        var_sb = f"{r['variant_sb']:.3f}" if r['variant_sb'] is not None else "N/A"
        print(f"{r['variant']:<15} {r['drug']:<12} {r['clinical']:<10} "
              f"{r['canonical_sb']:.3f}   {var_sb:<8} {r['delta_sb']:+.3f}   "
              f"{'PASS' if 'CORRECT' in r['interpretation'] else 'CHECK'}")

    # Save results
    output_path = PROJECT_ROOT / "data" / "processed" / "variant_predictions.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()

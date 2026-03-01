#!/usr/bin/env python3
"""
Full SpliceBind validation on splice variants.
Trains model and runs predictions on canonical vs variant structures.
"""

import sys
import numpy as np
import torch
import torch.nn as nn
from pathlib import Path
from typing import Dict, List, Optional

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

STRUCTURES_DIR = PROJECT_ROOT / "data" / "structures"
EMBEDDINGS_DIR = PROJECT_ROOT / "data" / "embeddings"


# Binding sites
KINASE_ATP_SITES = {
    'PIK3CD': {'start': 695, 'end': 1044},
    'BRAF': {'start': 457, 'end': 717},
    'EGFR': {'start': 712, 'end': 979},
    'ALK': {'start': 1116, 'end': 1399},
    'MET': {'start': 1078, 'end': 1345},
}


def load_esm_embeddings(embeddings_dir: Path) -> Dict[str, np.ndarray]:
    """Load all ESM-2 embeddings."""
    esm_embeddings = {}
    if embeddings_dir.exists():
        for emb_file in embeddings_dir.glob("*_esm2.npy"):
            protein_id = emb_file.stem.replace("_esm2", "")
            try:
                esm_embeddings[protein_id] = np.load(emb_file)
            except Exception:
                pass
    return esm_embeddings


def main():
    print("=" * 70)
    print("SpliceBind Full Validation")
    print("=" * 70)

    # Import components
    from modules.p2rank_detector import P2RankDetector
    from modules.structure_generator import StructureFeatureExtractor
    from modules.pocket_analyzer import PocketGraphBuilder
    from train import PocketGNNv5
    from torch_geometric.data import Batch, Data

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Load ESM embeddings
    print("\n[1] Loading ESM-2 embeddings...")
    esm_embeddings = load_esm_embeddings(EMBEDDINGS_DIR)
    print(f"    Loaded {len(esm_embeddings)} embeddings")

    # Also load variant-specific embeddings
    for emb_file in EMBEDDINGS_DIR.glob("*.npy"):
        if emb_file.stem not in esm_embeddings:
            try:
                esm_embeddings[emb_file.stem] = np.load(emb_file)
            except Exception:
                pass

    # Initialize components
    p2rank = P2RankDetector()
    feat_extractor = StructureFeatureExtractor()
    graph_builder = PocketGraphBuilder()

    # Validation cases
    CASES = {
        'PIK3CD_canonical': {
            'pdb': 'AF-O00329-F1-model_v6.pdb',
            'gene': 'PIK3CD',
            'expected': 'druggable',
        },
        'PIK3CD_variant': {
            'pdb': 'PIK3CD-S_esmfold.pdb',
            'gene': 'PIK3CD',
            'expected': 'reduced',
        },
        'BRAF_canonical': {
            'pdb': 'AF-P15056-F1-model_v6.pdb',
            'gene': 'BRAF',
            'expected': 'druggable',
        },
        'EGFR_canonical': {
            'pdb': 'AF-P00533-F1-model_v6.pdb',
            'gene': 'EGFR',
            'expected': 'druggable',
        },
        'ALK_canonical': {
            'pdb': 'AF-Q9UM73-F1-model_v4.pdb',
            'gene': 'ALK',
            'expected': 'druggable',
        },
    }

    # Train model first
    print("\n[2] Training model...")
    from train import main as train_main

    # Run training to get model
    import io
    import contextlib

    # Suppress output during training
    with contextlib.redirect_stdout(io.StringIO()):
        with contextlib.redirect_stderr(io.StringIO()):
            pass  # Skip full training for now

    # Load pre-computed pocket data instead
    pocket_data_path = PROJECT_ROOT / "data" / "processed" / "pocket_data.pt"
    if pocket_data_path.exists():
        print("    Loading pre-computed pocket data...")
        pocket_data = torch.load(pocket_data_path, weights_only=False)
        # pocket_data is a list of dicts with 'pocket_graph' key
        graphs = [item.get('pocket_graph') for item in pocket_data if item.get('pocket_graph')]
        labels = [item.get('label', 0) for item in pocket_data]
        print(f"    Loaded {len(graphs)} graphs")

        # Quick model training
        if graphs and labels:
            from torch_geometric.data import Data as PyGData
            from torch_geometric.loader import DataLoader

            # Convert to PyG format
            pyg_graphs = []
            for g, label in zip(graphs, labels):
                if isinstance(g, dict):
                    data = PyGData(
                        x=torch.tensor(g['node_features'], dtype=torch.float32),
                        edge_index=torch.tensor(g['edge_index'], dtype=torch.long),
                        y=torch.tensor([label], dtype=torch.float32)
                    )
                    if data.x.shape[1] == 56:
                        pyg_graphs.append(data)
                elif hasattr(g, 'x') and g.x.shape[1] == 56:
                    g.y = torch.tensor([label], dtype=torch.float32)
                    pyg_graphs.append(g)

            print(f"    Valid graphs: {len(pyg_graphs)}")

            if pyg_graphs:
                # Simple train/val split
                n_train = int(0.8 * len(pyg_graphs))
                train_graphs = pyg_graphs[:n_train]
                val_graphs = pyg_graphs[n_train:]

                train_loader = DataLoader(train_graphs, batch_size=32, shuffle=True)

                # Train model
                model = PocketGNNv5(node_dim=56, hidden_dim=256, num_layers=3)
                model = model.to(device)
                optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
                criterion = nn.BCEWithLogitsLoss()

                model.train()
                for epoch in range(30):
                    for batch in train_loader:
                        batch = batch.to(device)
                        optimizer.zero_grad()
                        logits = model(batch)
                        loss = criterion(logits.squeeze(), batch.y)
                        loss.backward()
                        optimizer.step()

                print(f"    Training complete")
    else:
        print("    No pre-computed data, using untrained model")
        model = PocketGNNv5(node_dim=56, hidden_dim=256, num_layers=3)
        model = model.to(device)

    model.eval()

    # Run predictions
    print("\n[3] Running predictions...")

    results = {}

    for case_name, case_info in CASES.items():
        pdb_path = STRUCTURES_DIR / case_info['pdb']
        if not pdb_path.exists():
            print(f"    {case_name}: Missing {pdb_path.name}")
            continue

        print(f"\n    {case_name}:")

        # Get pockets
        pockets = p2rank.detect_pockets(str(pdb_path))
        if not pockets:
            print(f"      No pockets detected")
            continue
        print(f"      Pockets: {len(pockets)}")

        # Get features
        features = feat_extractor.extract_features(str(pdb_path))
        ca_coords = features.get('ca_coords')
        if ca_coords is None:
            print(f"      No CA coords")
            continue
        if isinstance(ca_coords, list):
            ca_coords = np.array(ca_coords)

        # Find ESM embeddings
        pdb_stem = pdb_path.stem
        esm_emb = esm_embeddings.get(pdb_stem)
        if esm_emb is None:
            # Try alternate names
            for key in esm_embeddings:
                if case_info['gene'] in key or pdb_stem.split('-')[1] in key:
                    esm_emb = esm_embeddings[key]
                    break

        if esm_emb is None:
            print(f"      No ESM embeddings found")
            continue

        # Get binding site
        binding_site = KINASE_ATP_SITES.get(case_info['gene'])

        # Build graphs for binding site pockets
        graphs = []
        pocket_infos = []

        for pocket in pockets:
            pocket_residues = pocket.get('residue_indices', [])

            # Check binding site overlap
            if binding_site:
                site_residues = set(range(binding_site['start'], binding_site['end'] + 1))
                overlap = len(set(pocket_residues) & site_residues) / max(len(pocket_residues), 1)
                if overlap < 0.1:
                    continue

            graph_dict = graph_builder.build_pocket_graph(pocket, ca_coords, esm_emb)
            if graph_dict is None:
                continue

            # Convert to PyG Data object
            graph = Data(
                x=torch.tensor(graph_dict['node_features'], dtype=torch.float32),
                edge_index=torch.tensor(graph_dict['edge_index'], dtype=torch.long),
            )

            if graph.x.shape[1] == 56:
                graphs.append(graph)
                pocket_infos.append(pocket)

        if not graphs:
            print(f"      No valid binding site pockets")
            continue

        # Run model
        batch = Batch.from_data_list(graphs).to(device)
        with torch.no_grad():
            logits = model(batch)
            probs = torch.sigmoid(logits).cpu().numpy()

        # Get best binding site pocket score
        best_score = float(np.max(probs))
        best_idx = int(np.argmax(probs))
        p2rank_score = pocket_infos[best_idx].get('p2rank_score', 0)

        results[case_name] = {
            'score': best_score,
            'p2rank_score': p2rank_score,
            'gene': case_info['gene'],
            'expected': case_info['expected'],
        }

        print(f"      Druggability: {best_score:.3f}")
        print(f"      P2Rank: {p2rank_score:.2f}")

    # Summary
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)

    print(f"\n{'Case':<20} {'Gene':<8} {'Score':<8} {'P2Rank':<8} {'Expected'}")
    print("-" * 60)

    for case_name, r in results.items():
        print(f"{case_name:<20} {r['gene']:<8} {r['score']:.3f}    {r['p2rank_score']:.2f}     {r['expected']}")

    # PIK3CD comparison
    if 'PIK3CD_canonical' in results and 'PIK3CD_variant' in results:
        delta = results['PIK3CD_variant']['score'] - results['PIK3CD_canonical']['score']
        print(f"\nPIK3CD Delta: {delta:+.3f}")
        if delta < 0:
            print("CORRECT: Variant shows reduced druggability")
        else:
            print("INCORRECT: Expected reduced druggability")

    return results


if __name__ == "__main__":
    main()

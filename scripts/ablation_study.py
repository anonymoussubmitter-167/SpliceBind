#!/usr/bin/env python3
"""
Ablation Study for SpliceBind
Testing 4 conditions:
1. Structure-only baseline (24 features, no ESM-2)
2. ESM-2 only (32 features, no structural)
3. Random vectors replacing ESM-2 (56 features with random ESM)
4. Full model (56 features)
"""

import sys
import json
import numpy as np
import torch
import torch.nn as nn
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
from torch_geometric.data import Data, Batch
from torch_geometric.nn import global_mean_pool, global_max_pool
from torch_geometric.nn.conv import MessagePassing
from scipy import stats

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))


class EdgeConv(MessagePassing):
    """Edge convolution layer for point cloud processing."""

    def __init__(self, in_channels, out_channels):
        super().__init__(aggr='max')
        self.mlp = nn.Sequential(
            nn.Linear(in_channels * 2, out_channels),
            nn.ReLU(),
            nn.Linear(out_channels, out_channels)
        )

    def forward(self, x, edge_index):
        return self.propagate(edge_index, x=x)

    def message(self, x_i, x_j):
        return self.mlp(torch.cat([x_i, x_j - x_i], dim=-1))


class PocketGNN(nn.Module):
    """Pocket GNN model with configurable input dimension."""

    def __init__(self, node_dim, hidden_dim=256, num_layers=3, dropout=0.3):
        super().__init__()

        self.convs = nn.ModuleList()
        self.convs.append(EdgeConv(node_dim, hidden_dim))
        for _ in range(num_layers - 1):
            self.convs.append(EdgeConv(hidden_dim, hidden_dim))

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
            x = torch.relu(x)

        # Global pooling
        x_mean = global_mean_pool(x, batch)
        x_max = global_max_pool(x, batch)
        x = torch.cat([x_mean, x_max], dim=1)

        return self.classifier(x).squeeze(-1)


class FocalLoss(nn.Module):
    """Focal loss for handling class imbalance."""

    def __init__(self, gamma=2.0, pos_weight=1.0):
        super().__init__()
        self.gamma = gamma
        self.pos_weight = pos_weight

    def forward(self, logits, targets):
        bce = nn.functional.binary_cross_entropy_with_logits(
            logits, targets, pos_weight=torch.tensor(self.pos_weight).to(logits.device),
            reduction='none'
        )
        pt = torch.exp(-bce)
        focal_loss = ((1 - pt) ** self.gamma) * bce
        return focal_loss.mean()


def load_pocket_data():
    """Load pocket data from processed files."""
    pocket_path = PROJECT_ROOT / "data" / "processed" / "pocket_data.pt"
    if not pocket_path.exists():
        raise FileNotFoundError(f"Pocket data not found at {pocket_path}")

    data = torch.load(pocket_path, weights_only=False)

    # Load ESM embeddings
    embeddings_dir = PROJECT_ROOT / "data" / "embeddings"
    esm_embeddings = {}
    if embeddings_dir.exists():
        for emb_file in embeddings_dir.glob("*_esm2.npy"):
            protein_id = emb_file.stem.replace("_esm2", "")
            try:
                esm_embeddings[protein_id] = np.load(emb_file)
            except Exception:
                pass

    # Gene to family mapping
    GENE_TO_FAMILY = {
        'PIM1': 'PIM', 'PIM2': 'PIM', 'PIM3': 'PIM',
        'AKT1': 'AKT', 'AKT2': 'AKT', 'AKT3': 'AKT',
        'PIK3CA': 'PIK3', 'PIK3CB': 'PIK3', 'PIK3CG': 'PIK3', 'PIK3CD': 'PIK3',
        'CLK1': 'CLK', 'CLK2': 'CLK', 'CLK3': 'CLK', 'CLK4': 'CLK',
        'JAK1': 'JAK', 'JAK2': 'JAK', 'JAK3': 'JAK', 'TYK2': 'JAK',
        'ABL1': 'ABL', 'ABL2': 'ABL',
        'EGFR': 'EGFR', 'ERBB2': 'EGFR', 'ERBB3': 'EGFR', 'ERBB4': 'EGFR',
        'BRAF': 'RAF', 'RAF1': 'RAF', 'ARAF': 'RAF',
        'SRC': 'SRC', 'LCK': 'SRC', 'FYN': 'SRC', 'YES1': 'SRC', 'HCK': 'SRC',
        'CDK1': 'CDK', 'CDK2': 'CDK', 'CDK4': 'CDK', 'CDK6': 'CDK',
        'MAPK1': 'MAPK', 'MAPK3': 'MAPK', 'MAPK14': 'MAPK',
        'AURKA': 'AURK', 'AURKB': 'AURK',
        'PLK1': 'PLK', 'PLK4': 'PLK',
        'FGFR1': 'FGFR', 'FGFR2': 'FGFR', 'FGFR3': 'FGFR', 'FGFR4': 'FGFR',
        'KDR': 'VEGFR', 'FLT1': 'VEGFR', 'FLT4': 'VEGFR',
        'PDGFRA': 'PDGFR', 'PDGFRB': 'PDGFR',
        'MET': 'MET', 'ALK': 'ALK',
        'CHEK1': 'CHEK', 'CHEK2': 'CHEK',
        'WEE1': 'WEE', 'GSK3B': 'GSK', 'GSK3A': 'GSK',
        'ROCK1': 'ROCK', 'ROCK2': 'ROCK',
        'BTK': 'TEC', 'ITK': 'TEC', 'TXK': 'TEC',
        'SYK': 'SYK', 'ZAP70': 'SYK',
        'MTOR': 'MTOR',
    }

    # Assign family and ESM embeddings to pockets
    for pocket in data:
        gene = pocket.get('gene', 'UNKNOWN')
        pocket['family'] = GENE_TO_FAMILY.get(gene, 'OTHER')

        # Try to find ESM embedding for this pocket
        uniprot = pocket.get('uniprot', '')
        esm_emb = None

        # Try different key patterns
        for key in esm_embeddings:
            if uniprot in key:
                esm_emb = esm_embeddings[key]
                break

        pocket['esm_embedding'] = esm_emb

    return data


def assign_labels_by_druggability(pockets, threshold_high=0.7, threshold_low=0.3):
    """Assign binary labels based on druggability score."""
    for pocket in pockets:
        score = pocket.get('druggability_score', 0.5)
        if score >= threshold_high:
            pocket['label'] = 1
        elif score <= threshold_low:
            pocket['label'] = 0
        else:
            pocket['label'] = None  # Ambiguous, will be filtered

    return [p for p in pockets if p.get('label') is not None]


def create_ablation_graphs(pockets, condition):
    """Create graphs with modified features based on ablation condition.

    Conditions:
    - 'full': All 56 features (24 structural + 32 ESM-2)
    - 'structure_only': Only 24 structural features
    - 'esm_only': Only 32 ESM-2 features
    - 'random_esm': 24 structural + 32 random features
    """
    graphs = []
    labels = []
    families = []

    for pocket in pockets:
        # Get the pocket graph from the correct key
        pocket_graph = pocket.get('pocket_graph')
        if pocket_graph is None:
            continue

        # Get node features
        node_features = pocket_graph.get('node_features')
        if node_features is None or len(node_features) == 0:
            continue

        if isinstance(node_features, np.ndarray):
            node_features = torch.tensor(node_features, dtype=torch.float32)

        n_nodes = node_features.shape[0]
        struct_dim = node_features.shape[1]  # Should be 24

        # Get ESM embedding if available
        esm_emb = pocket.get('esm_embedding')
        has_esm = esm_emb is not None

        # Create features based on condition
        if condition == 'full':
            if has_esm:
                # Reduce ESM to 32 dims using binned pooling
                esm_reduced = np.zeros((n_nodes, 32))
                for i in range(n_nodes):
                    # Use mean pooling across bins
                    for j in range(32):
                        start = j * 40
                        end = start + 40
                        if i < len(esm_emb) and end <= esm_emb.shape[1]:
                            esm_reduced[i, j] = esm_emb[i, start:end].mean()
                esm_tensor = torch.tensor(esm_reduced, dtype=torch.float32)
                new_x = torch.cat([node_features, esm_tensor], dim=1)
            else:
                # No ESM available, skip this pocket for full model
                continue

        elif condition == 'structure_only':
            new_x = node_features

        elif condition == 'esm_only':
            if has_esm:
                esm_reduced = np.zeros((n_nodes, 32))
                for i in range(n_nodes):
                    for j in range(32):
                        start = j * 40
                        end = start + 40
                        if i < len(esm_emb) and end <= esm_emb.shape[1]:
                            esm_reduced[i, j] = esm_emb[i, start:end].mean()
                new_x = torch.tensor(esm_reduced, dtype=torch.float32)
            else:
                continue

        elif condition == 'random_esm':
            # Keep structural features, add random ESM
            random_esm = torch.randn(n_nodes, 32)
            new_x = torch.cat([node_features, random_esm], dim=1)

        else:
            raise ValueError(f"Unknown condition: {condition}")

        # Get edge info
        edge_index = pocket_graph.get('edge_index')
        edge_attr = pocket_graph.get('edge_attr')

        if isinstance(edge_index, np.ndarray):
            edge_index = torch.tensor(edge_index, dtype=torch.long)
        if edge_attr is not None and isinstance(edge_attr, np.ndarray):
            edge_attr = torch.tensor(edge_attr, dtype=torch.float32)

        new_graph = Data(
            x=new_x.float(),
            edge_index=edge_index,
            edge_attr=edge_attr,
            num_nodes=n_nodes
        )

        graphs.append(new_graph)
        labels.append(pocket.get('label', 0))
        families.append(pocket.get('family', 'OTHER'))

    return graphs, np.array(labels), np.array(families)


def train_and_evaluate(graphs, labels, families, node_dim, device, n_seeds=3):
    """Train and evaluate model with cross-validation."""
    all_results = []

    for seed in range(n_seeds):
        np.random.seed(seed)
        torch.manual_seed(seed)

        kfold = GroupKFold(n_splits=5)

        for fold, (train_idx, test_idx) in enumerate(kfold.split(graphs, labels, families)):
            train_labels = labels[train_idx]
            test_labels = labels[test_idx]

            # Skip degenerate folds
            if len(np.unique(train_labels)) < 2 or len(np.unique(test_labels)) < 2:
                continue

            # Split train into train/val
            val_size = max(1, len(train_idx) // 5)
            val_idx = train_idx[:val_size]
            train_idx_final = train_idx[val_size:]

            train_graphs = [graphs[i] for i in train_idx_final]
            val_graphs = [graphs[i] for i in val_idx]
            test_graphs = [graphs[i] for i in test_idx]

            train_labels_final = labels[train_idx_final]
            val_labels_fold = labels[val_idx]

            # Initialize model
            model = PocketGNN(node_dim=node_dim).to(device)

            # Calculate class weights
            pos_rate = train_labels_final.mean()
            pos_weight = min((1 - pos_rate) / max(pos_rate, 0.01), 10.0)

            criterion = FocalLoss(gamma=2.0, pos_weight=pos_weight)
            optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)

            # Training
            best_val_auroc = 0
            best_state = None
            patience = 0

            train_labels_t = torch.tensor(train_labels_final, dtype=torch.float32, device=device)
            val_labels_t = torch.tensor(val_labels_fold, dtype=torch.float32, device=device)

            for epoch in range(50):
                model.train()
                train_batch = Batch.from_data_list(train_graphs).to(device)

                optimizer.zero_grad()
                logits = model(train_batch)
                loss = criterion(logits, train_labels_t)
                loss.backward()
                optimizer.step()

                # Validation
                model.eval()
                with torch.no_grad():
                    val_batch = Batch.from_data_list(val_graphs).to(device)
                    val_logits = model(val_batch)
                    val_probs = torch.sigmoid(val_logits).cpu().numpy()

                if len(np.unique(val_labels_fold)) > 1:
                    val_auroc = roc_auc_score(val_labels_fold, val_probs)
                else:
                    val_auroc = 0.5

                if val_auroc > best_val_auroc:
                    best_val_auroc = val_auroc
                    best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
                    patience = 0
                else:
                    patience += 1

                if patience >= 10:
                    break

            # Load best model and evaluate
            if best_state:
                model.load_state_dict({k: v.to(device) for k, v in best_state.items()})

            model.eval()
            with torch.no_grad():
                test_batch = Batch.from_data_list(test_graphs).to(device)
                test_logits = model(test_batch)
                test_probs = torch.sigmoid(test_logits).cpu().numpy()

            if len(np.unique(test_labels)) > 1:
                auroc = roc_auc_score(test_labels, test_probs)
                auprc = average_precision_score(test_labels, test_probs)
            else:
                auroc = 0.5
                auprc = 0.5

            all_results.append({
                'seed': seed,
                'fold': fold,
                'auroc': auroc,
                'auprc': auprc,
                'n_test': len(test_idx)
            })

    return all_results


def main():
    print("=" * 70)
    print("ABLATION STUDY - SpliceBind Feature Analysis")
    print("=" * 70)

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Device: {device}")

    # Load data
    print("\nLoading pocket data...")
    pockets = load_pocket_data()
    print(f"Loaded {len(pockets)} pockets")

    # Count ESM embeddings
    esm_count = sum(1 for p in pockets if p.get('esm_embedding') is not None)
    print(f"Pockets with ESM embeddings: {esm_count}")

    # Assign labels based on druggability
    pockets = assign_labels_by_druggability(pockets)
    print(f"Pockets with valid labels: {len(pockets)}")
    print(f"Positive rate: {np.mean([p['label'] for p in pockets]):.1%}")

    # Define ablation conditions
    conditions = {
        'full': {'desc': 'Full model (56 features)', 'dim': 56},
        'structure_only': {'desc': 'Structure only (24 features)', 'dim': 24},
        'esm_only': {'desc': 'ESM-2 only (32 features)', 'dim': 32},
        'random_esm': {'desc': 'Random ESM (56 features)', 'dim': 56}
    }

    results = {
        'timestamp': datetime.now().isoformat(),
        'conditions': {}
    }

    for condition, info in conditions.items():
        print(f"\n{'=' * 70}")
        print(f"Testing: {info['desc']}")
        print("=" * 70)

        # Create graphs with modified features
        graphs, labels, families = create_ablation_graphs(pockets, condition)

        if len(graphs) == 0:
            print(f"No valid graphs for condition {condition}")
            continue

        print(f"Created {len(graphs)} graphs")
        print(f"Positive rate: {labels.mean():.1%}")

        # Determine actual feature dimension
        actual_dim = graphs[0].x.shape[1]
        print(f"Feature dimension: {actual_dim}")

        # Train and evaluate
        cv_results = train_and_evaluate(graphs, labels, families, actual_dim, device)

        if len(cv_results) == 0:
            print(f"No valid CV results for condition {condition}")
            continue

        aurocs = [r['auroc'] for r in cv_results]
        auprcs = [r['auprc'] for r in cv_results]

        results['conditions'][condition] = {
            'description': info['desc'],
            'feature_dim': actual_dim,
            'n_folds': len(cv_results),
            'n_graphs': len(graphs),
            'auroc_mean': float(np.mean(aurocs)),
            'auroc_std': float(np.std(aurocs)),
            'auprc_mean': float(np.mean(auprcs)),
            'auprc_std': float(np.std(auprcs)),
            'per_fold': cv_results
        }

        print(f"AUROC: {np.mean(aurocs):.3f} ± {np.std(aurocs):.3f}")
        print(f"AUPRC: {np.mean(auprcs):.3f} ± {np.std(auprcs):.3f}")

    # Statistical comparisons
    print("\n" + "=" * 70)
    print("STATISTICAL COMPARISONS")
    print("=" * 70)

    comparisons = []

    if 'full' in results['conditions'] and 'structure_only' in results['conditions']:
        full_aurocs = [r['auroc'] for r in results['conditions']['full']['per_fold']]
        struct_aurocs = [r['auroc'] for r in results['conditions']['structure_only']['per_fold']]

        # Paired t-test (where folds match)
        n = min(len(full_aurocs), len(struct_aurocs))
        if n > 2:
            t_stat, p_value = stats.ttest_rel(full_aurocs[:n], struct_aurocs[:n])
            delta = np.mean(full_aurocs) - np.mean(struct_aurocs)

            comp = {
                'comparison': 'Full vs Structure-only',
                'delta': float(delta),
                'p_value': float(p_value),
                'significant': p_value < 0.05
            }
            comparisons.append(comp)
            print(f"\nFull vs Structure-only: Δ={delta:+.3f}, p={p_value:.4f}")

    if 'full' in results['conditions'] and 'random_esm' in results['conditions']:
        full_aurocs = [r['auroc'] for r in results['conditions']['full']['per_fold']]
        random_aurocs = [r['auroc'] for r in results['conditions']['random_esm']['per_fold']]

        n = min(len(full_aurocs), len(random_aurocs))
        if n > 2:
            t_stat, p_value = stats.ttest_rel(full_aurocs[:n], random_aurocs[:n])
            delta = np.mean(full_aurocs) - np.mean(random_aurocs)

            comp = {
                'comparison': 'Full vs Random-ESM',
                'delta': float(delta),
                'p_value': float(p_value),
                'significant': p_value < 0.05
            }
            comparisons.append(comp)
            print(f"Full vs Random-ESM: Δ={delta:+.3f}, p={p_value:.4f}")

    if 'esm_only' in results['conditions'] and 'structure_only' in results['conditions']:
        esm_aurocs = [r['auroc'] for r in results['conditions']['esm_only']['per_fold']]
        struct_aurocs = [r['auroc'] for r in results['conditions']['structure_only']['per_fold']]

        n = min(len(esm_aurocs), len(struct_aurocs))
        if n > 2:
            t_stat, p_value = stats.ttest_rel(esm_aurocs[:n], struct_aurocs[:n])
            delta = np.mean(esm_aurocs) - np.mean(struct_aurocs)

            comp = {
                'comparison': 'ESM-only vs Structure-only',
                'delta': float(delta),
                'p_value': float(p_value),
                'significant': p_value < 0.05
            }
            comparisons.append(comp)
            print(f"ESM-only vs Structure-only: Δ={delta:+.3f}, p={p_value:.4f}")

    results['comparisons'] = comparisons

    # Summary table
    print("\n" + "=" * 70)
    print("ABLATION STUDY SUMMARY")
    print("=" * 70)
    print(f"\n{'Condition':<25} {'Features':<10} {'AUROC':<18} {'AUPRC':<18}")
    print("-" * 70)

    for condition in ['full', 'structure_only', 'esm_only', 'random_esm']:
        if condition in results['conditions']:
            c = results['conditions'][condition]
            print(f"{c['description']:<25} {c['feature_dim']:<10} "
                  f"{c['auroc_mean']:.3f} ± {c['auroc_std']:.3f}    "
                  f"{c['auprc_mean']:.3f} ± {c['auprc_std']:.3f}")

    # ESM contribution
    if 'full' in results['conditions'] and 'structure_only' in results['conditions']:
        esm_contrib = results['conditions']['full']['auroc_mean'] - results['conditions']['structure_only']['auroc_mean']
        print(f"\nESM-2 contribution to AUROC: {esm_contrib:+.3f}")

    # Save results
    output_path = PROJECT_ROOT / "data" / "processed" / "ablation_results.json"
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=float)

    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()

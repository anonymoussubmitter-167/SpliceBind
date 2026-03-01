#!/usr/bin/env python3
"""
External Validation: Hold-out kinase family evaluation.

Train on N-5 families, test on 5 held-out families.
More rigorous than GroupKFold (complete family exclusion).
"""

import sys
import json
import numpy as np
import torch
import torch.nn as nn
from pathlib import Path
from typing import Dict, List
from collections import defaultdict
from sklearn.metrics import roc_auc_score
from torch_geometric.data import Data
from torch_geometric.nn import EdgeConv, global_mean_pool, global_max_pool
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))


class PocketGNNv5(nn.Module):
    def __init__(self, node_dim=24, hidden_dim=256, num_layers=3, dropout=0.3):
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


def load_data():
    """Load pocket data with family annotations."""
    data_path = PROJECT_ROOT / "data" / "processed" / "pocket_data.pt"
    data = torch.load(data_path, weights_only=False)

    # Group by family (extract from gene name)
    GENE_TO_FAMILY = {
        'JAK1': 'JAK', 'JAK2': 'JAK', 'JAK3': 'JAK', 'TYK2': 'JAK',
        'EGFR': 'EGFR', 'ERBB2': 'EGFR', 'ERBB4': 'EGFR',
        'ABL1': 'ABL', 'ABL2': 'ABL',
        'BRAF': 'RAF', 'RAF1': 'RAF', 'ARAF': 'RAF',
        'PIK3CA': 'PIK3', 'PIK3CB': 'PIK3', 'PIK3CD': 'PIK3', 'PIK3CG': 'PIK3',
        'ALK': 'ALK', 'ROS1': 'ALK', 'LTK': 'ALK',
        'MET': 'MET', 'MST1R': 'MET',
        'KIT': 'RTK', 'FLT3': 'RTK', 'PDGFRA': 'RTK', 'PDGFRB': 'RTK', 'CSF1R': 'RTK',
        'VEGFR1': 'VEGFR', 'VEGFR2': 'VEGFR', 'VEGFR3': 'VEGFR',
        'SRC': 'SRC', 'LCK': 'SRC', 'LYN': 'SRC', 'FYN': 'SRC',
        'BTK': 'TEC', 'ITK': 'TEC', 'TEC': 'TEC',
        'AURKA': 'AURK', 'AURKB': 'AURK', 'AURKC': 'AURK',
        'CDK1': 'CDK', 'CDK2': 'CDK', 'CDK4': 'CDK', 'CDK6': 'CDK',
        'RET': 'RET',
        'NTRK1': 'NTRK', 'NTRK2': 'NTRK', 'NTRK3': 'NTRK',
        'FGFR1': 'FGFR', 'FGFR2': 'FGFR', 'FGFR3': 'FGFR', 'FGFR4': 'FGFR',
        'MEK1': 'MEK', 'MEK2': 'MEK',
    }

    samples = []
    for pocket in data:
        if pocket.get('pocket_graph') is None:
            continue

        gene = pocket.get('gene', 'Unknown')
        family = GENE_TO_FAMILY.get(gene, gene)  # Use gene name if not in mapping

        graph = pocket['pocket_graph']
        score = pocket.get('druggability_score', 0)
        label = 1 if score > 0.5 else 0

        samples.append({
            'gene': gene,
            'family': family,
            'graph': graph,
            'label': label,
        })

    return samples


def train_model(train_data: List[Dict], device: torch.device, epochs: int = 30) -> nn.Module:
    """Train model on training data."""
    model = PocketGNNv5(node_dim=24).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)
    criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor([5.0]).to(device))

    dataset = []
    for s in train_data:
        g = s['graph']
        x = torch.tensor(g['node_features'], dtype=torch.float32)
        edge_index = torch.tensor(g['edge_index'], dtype=torch.long)
        y = torch.tensor([s['label']], dtype=torch.float32)
        dataset.append(Data(x=x, edge_index=edge_index, y=y))

    model.train()
    for epoch in range(epochs):
        np.random.shuffle(dataset)
        total_loss = 0
        for data in dataset:
            data = data.to(device)
            optimizer.zero_grad()
            out = model(data)
            loss = criterion(out, data.y)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

    return model


def evaluate_model(model: nn.Module, test_data: List[Dict], device: torch.device) -> Dict:
    """Evaluate model on test data."""
    model.eval()

    predictions = []
    labels = []

    with torch.no_grad():
        for s in test_data:
            g = s['graph']
            x = torch.tensor(g['node_features'], dtype=torch.float32).to(device)
            edge_index = torch.tensor(g['edge_index'], dtype=torch.long).to(device)
            batch = torch.zeros(x.size(0), dtype=torch.long).to(device)

            data = Data(x=x, edge_index=edge_index, batch=batch)
            logit = model(data)
            prob = torch.sigmoid(logit).item()

            predictions.append(prob)
            labels.append(s['label'])

    if len(set(labels)) < 2:
        auroc = 0.5
    else:
        auroc = roc_auc_score(labels, predictions)

    return {
        'auroc': auroc,
        'n_samples': len(test_data),
        'n_positive': sum(labels),
    }


def main():
    print("=" * 70)
    print("EXTERNAL VALIDATION: Hold-out Family Evaluation")
    print("=" * 70)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Load data
    samples = load_data()
    print(f"Total samples: {len(samples)}")

    # Group by family
    family_samples = defaultdict(list)
    for s in samples:
        family_samples[s['family']].append(s)

    print(f"\nFamilies: {len(family_samples)}")
    for fam, samps in sorted(family_samples.items(), key=lambda x: -len(x[1])):
        print(f"  {fam}: {len(samps)} samples")

    # Select top 5 families for holdout testing
    top_families = sorted(family_samples.keys(), key=lambda x: -len(family_samples[x]))[:5]
    print(f"\nHold-out families: {top_families}")

    results = []

    for holdout_family in top_families:
        print(f"\n{'='*60}")
        print(f"Hold-out: {holdout_family}")
        print(f"{'='*60}")

        # Split data
        train_data = [s for s in samples if s['family'] != holdout_family]
        test_data = [s for s in samples if s['family'] == holdout_family]

        print(f"  Train: {len(train_data)} samples ({len(set(s['family'] for s in train_data))} families)")
        print(f"  Test: {len(test_data)} samples")

        # Train and evaluate
        model = train_model(train_data, device)
        result = evaluate_model(model, test_data, device)

        print(f"  AUROC: {result['auroc']:.3f}")
        print(f"  Positive rate: {result['n_positive']}/{result['n_samples']}")

        results.append({
            'holdout_family': holdout_family,
            **result,
        })

    # Summary
    print("\n" + "=" * 70)
    print("EXTERNAL VALIDATION SUMMARY")
    print("=" * 70)

    aurocs = [r['auroc'] for r in results]
    print(f"\n{'Family':<12} {'AUROC':<10} {'N':<8} {'Pos':<8}")
    print("-" * 40)

    for r in results:
        print(f"{r['holdout_family']:<12} {r['auroc']:<10.3f} {r['n_samples']:<8} {r['n_positive']:<8}")

    print("-" * 40)
    print(f"{'Mean':<12} {np.mean(aurocs):<10.3f}")
    print(f"{'Std':<12} {np.std(aurocs):<10.3f}")

    # Save results
    output_path = PROJECT_ROOT / "data" / "processed" / "external_validation.json"
    with open(output_path, 'w') as f:
        json.dump({
            'method': 'hold-out family validation',
            'holdout_families': top_families,
            'mean_auroc': float(np.mean(aurocs)),
            'std_auroc': float(np.std(aurocs)),
            'results': results,
        }, f, indent=2)

    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()

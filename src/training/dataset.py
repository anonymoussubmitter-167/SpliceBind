"""Dataset classes for SpliceBind training.

Handles creation of pocket graph datasets from processed data,
including train/val/test splitting by kinase family.
"""

import os
import json
import logging
import numpy as np
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Data, InMemoryDataset
from typing import Dict, List, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)


class PocketDruggabilityDataset(Dataset):
    """Dataset for pocket druggability prediction.

    Each sample is a pocket graph with a druggability label.
    """

    def __init__(self, pocket_graphs: List[Dict], labels: List[float],
                 metadata: Optional[List[Dict]] = None):
        self.pocket_graphs = pocket_graphs
        self.labels = labels
        self.metadata = metadata or [{}] * len(labels)
        assert len(pocket_graphs) == len(labels)

    def __len__(self):
        return len(self.pocket_graphs)

    def __getitem__(self, idx) -> Data:
        graph = self.pocket_graphs[idx]
        label = self.labels[idx]

        data = Data(
            x=graph["node_features"].float(),
            edge_index=graph["edge_index"].long(),
            edge_attr=graph["edge_attr"].float(),
            num_nodes=graph["num_nodes"],
            y=torch.tensor([label], dtype=torch.float32),
        )

        # Store metadata
        if self.metadata[idx]:
            data.gene = self.metadata[idx].get("gene", "")
            data.family = self.metadata[idx].get("family", "")
            data.pocket_idx = self.metadata[idx].get("pocket_idx", 0)

        return data


class IsoformPairDataset(Dataset):
    """Dataset for isoform pair comparison.

    Each sample is a pair of pocket graphs (canonical, isoform)
    with a pocket conservation label.
    """

    def __init__(self, canonical_graphs: List[Dict], isoform_graphs: List[Dict],
                 similarity_labels: List[float],
                 structural_features: List[np.ndarray],
                 metadata: Optional[List[Dict]] = None):
        self.canonical_graphs = canonical_graphs
        self.isoform_graphs = isoform_graphs
        self.similarity_labels = similarity_labels
        self.structural_features = structural_features
        self.metadata = metadata or [{}] * len(similarity_labels)

    def __len__(self):
        return len(self.similarity_labels)

    def __getitem__(self, idx) -> Dict:
        can_graph = self.canonical_graphs[idx]
        iso_graph = self.isoform_graphs[idx]

        can_data = Data(
            x=can_graph["node_features"].float(),
            edge_index=can_graph["edge_index"].long(),
            edge_attr=can_graph["edge_attr"].float(),
            num_nodes=can_graph["num_nodes"],
        )

        iso_data = Data(
            x=iso_graph["node_features"].float(),
            edge_index=iso_graph["edge_index"].long(),
            edge_attr=iso_graph["edge_attr"].float(),
            num_nodes=iso_graph["num_nodes"],
        )

        return {
            "canonical": can_data,
            "isoform": iso_data,
            "similarity": torch.tensor([self.similarity_labels[idx]], dtype=torch.float32),
            "structural_features": torch.tensor(self.structural_features[idx], dtype=torch.float32),
        }


def create_training_data_from_structures(structure_dir: str,
                                          binding_data_path: str,
                                          output_dir: str) -> Dict:
    """Create training datasets from protein structures and binding data.

    This is the main data preparation function that:
    1. Loads protein structures
    2. Detects pockets using GeometricPocketDetector
    3. Builds pocket graphs using PocketGraphBuilder
    4. Labels pockets as druggable/non-druggable based on binding data
    5. Splits by kinase family

    Returns dict with train/val/test datasets.
    """
    import pandas as pd
    from ..modules.pocket_analyzer import GeometricPocketDetector, PocketGraphBuilder

    os.makedirs(output_dir, exist_ok=True)
    cache_path = os.path.join(output_dir, "training_data.pt")

    if os.path.exists(cache_path):
        logger.info(f"Loading cached training data from {cache_path}")
        return torch.load(cache_path, weights_only=False)

    # Load binding data
    binding_df = pd.read_csv(binding_data_path)
    logger.info(f"Loaded binding data: {len(binding_df)} entries")

    # Initialize modules
    detector = GeometricPocketDetector(contact_distance=6.0, min_pocket_size=5)
    graph_builder = PocketGraphBuilder(contact_distance=6.0)

    all_graphs = []
    all_labels = []
    all_metadata = []

    # Process each structure
    structure_files = [f for f in os.listdir(structure_dir) if f.endswith(".pdb")]
    logger.info(f"Processing {len(structure_files)} structure files...")

    for pdb_file in structure_files:
        pdb_path = os.path.join(structure_dir, pdb_file)

        # Extract UniProt ID from filename (AF-XXXXX-F1-model_v4.pdb)
        parts = pdb_file.split("-")
        if len(parts) >= 2:
            uniprot_id = parts[1]
        else:
            continue

        # Find gene name for this UniProt ID
        gene_match = binding_df[binding_df["uniprot"] == uniprot_id]
        if gene_match.empty:
            continue
        gene = gene_match.iloc[0]["gene"]
        family = gene_match.iloc[0]["family"]

        logger.info(f"Processing {gene} ({uniprot_id})...")

        # Detect pockets
        pockets = detector.detect_pockets(pdb_path)
        if not pockets:
            logger.warning(f"No pockets detected for {gene}")
            continue

        # Parse coordinates for graph building
        coords = []
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    coords.append([
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ])
        coords = np.array(coords)

        # Label pockets
        has_binding = len(gene_match) > 0 and gene_match["pKd"].max() > 5.0

        for pi, pocket in enumerate(pockets):
            # Build pocket graph
            pocket_graph = graph_builder.build_pocket_graph(pocket, coords)
            if pocket_graph is None:
                continue

            # Label: top pocket of a drug target = druggable
            # Other pockets = non-druggable (simplified)
            if pi == 0 and has_binding:
                label = 1.0
            elif pi == 0:
                label = 0.5  # might be druggable but unknown
            else:
                label = 0.0  # non-primary pocket

            all_graphs.append(pocket_graph)
            all_labels.append(label)
            all_metadata.append({
                "gene": gene,
                "family": family,
                "uniprot": uniprot_id,
                "pocket_idx": pi,
                "pocket_score": pocket.get("druggability_score", 0),
                "num_residues": pocket.get("num_residues", 0),
            })

    logger.info(f"Total pocket graphs: {len(all_graphs)}")
    logger.info(f"  Druggable (label=1): {sum(1 for l in all_labels if l == 1.0)}")
    logger.info(f"  Unknown (label=0.5): {sum(1 for l in all_labels if l == 0.5)}")
    logger.info(f"  Non-druggable (label=0): {sum(1 for l in all_labels if l == 0.0)}")

    # Split by family
    train_families = {"PIM", "AKT", "PIK3"}
    test_families = {"CLK", "JAK"}

    train_indices = [i for i, m in enumerate(all_metadata) if m["family"] in train_families]
    test_indices = [i for i, m in enumerate(all_metadata) if m["family"] in test_families]

    # Split train into train/val
    np.random.seed(42)
    np.random.shuffle(train_indices)
    val_size = max(1, len(train_indices) // 5)
    val_indices = train_indices[:val_size]
    train_indices = train_indices[val_size:]

    def make_dataset(indices):
        return PocketDruggabilityDataset(
            pocket_graphs=[all_graphs[i] for i in indices],
            labels=[all_labels[i] for i in indices],
            metadata=[all_metadata[i] for i in indices],
        )

    datasets = {
        "train": make_dataset(train_indices),
        "val": make_dataset(val_indices),
        "test": make_dataset(test_indices),
        "all_graphs": all_graphs,
        "all_labels": all_labels,
        "all_metadata": all_metadata,
    }

    # Log split info
    logger.info(f"Data splits:")
    logger.info(f"  Train: {len(train_indices)} samples")
    logger.info(f"  Val: {len(val_indices)} samples")
    logger.info(f"  Test: {len(test_indices)} samples")

    # Save
    torch.save(datasets, cache_path)
    logger.info(f"Saved training data to {cache_path}")

    return datasets


def create_synthetic_training_data(num_samples: int = 500, seed: int = 42) -> Dict:
    """Create synthetic training data for initial model development.

    Generates realistic pocket graphs with known druggability labels
    for testing the training pipeline before real data is available.
    """
    np.random.seed(seed)

    all_graphs = []
    all_labels = []
    all_metadata = []

    families = ["PIM", "AKT", "PIK3", "CLK", "JAK"]
    genes_per_family = {
        "PIM": ["PIM1", "PIM2", "PIM3"],
        "AKT": ["AKT1", "AKT2", "AKT3"],
        "PIK3": ["PIK3CA", "PIK3CB", "PIK3CG", "PIK3CD"],
        "CLK": ["CLK1", "CLK2", "CLK3"],
        "JAK": ["JAK1", "JAK2", "JAK3"],
    }

    node_feature_dim = 24  # 20 (AA one-hot) + 3 (properties) + 1 (pLDDT)

    samples_per_gene = num_samples // sum(len(g) for g in genes_per_family.values())

    for family in families:
        for gene in genes_per_family[family]:
            for i in range(max(1, samples_per_gene)):
                # Generate pocket of random size
                num_nodes = np.random.randint(6, 25)

                # Node features: one-hot AA + properties + pLDDT
                aa_indices = np.random.randint(0, 20, num_nodes)
                one_hot = np.eye(20)[aa_indices]
                properties = np.random.randn(num_nodes, 3) * 0.5
                plddt = np.random.uniform(50, 95, (num_nodes, 1)) / 100.0
                node_features = np.concatenate([one_hot, properties, plddt], axis=1).astype(np.float32)

                # Edge construction: random graph with spatial structure
                num_edges = np.random.randint(num_nodes, num_nodes * 3)
                edge_index = np.random.randint(0, num_nodes, (2, num_edges))
                # Remove self-loops
                mask = edge_index[0] != edge_index[1]
                edge_index = edge_index[:, mask]

                edge_attr = np.random.rand(edge_index.shape[1], 2).astype(np.float32)
                edge_attr[:, 0] = np.clip(edge_attr[:, 0], 0, 1)  # normalized distance
                edge_attr[:, 1] = np.clip(edge_attr[:, 1], 0, 1)  # sequence separation

                # Druggability label: influenced by pocket properties
                # More hydrophobic, moderate size → more druggable
                hydro_score = properties[:, 0].mean()
                size_score = 1.0 if 10 <= num_nodes <= 20 else 0.5
                plddt_score = plddt.mean()

                # Primary pocket for known drug targets is druggable
                is_primary = (i == 0)
                base_drug = 0.7 if is_primary else 0.3
                label = min(1.0, max(0.0,
                    base_drug + 0.1 * hydro_score + 0.1 * size_score + 0.1 * plddt_score
                    + np.random.normal(0, 0.05)
                ))

                # Binarize
                label = 1.0 if label > 0.5 else 0.0

                graph = {
                    "node_features": torch.tensor(node_features),
                    "edge_index": torch.tensor(edge_index, dtype=torch.long),
                    "edge_attr": torch.tensor(edge_attr),
                    "num_nodes": num_nodes,
                }

                all_graphs.append(graph)
                all_labels.append(label)
                all_metadata.append({
                    "gene": gene,
                    "family": family,
                    "pocket_idx": i,
                    "synthetic": True,
                })

    # Split by family
    train_families = {"PIM", "AKT", "PIK3"}
    test_families = {"CLK", "JAK"}

    train_indices = [i for i, m in enumerate(all_metadata) if m["family"] in train_families]
    test_indices = [i for i, m in enumerate(all_metadata) if m["family"] in test_families]

    np.random.shuffle(train_indices)
    val_size = max(1, len(train_indices) // 5)
    val_indices = train_indices[:val_size]
    train_indices = train_indices[val_size:]

    def make_dataset(indices):
        return PocketDruggabilityDataset(
            pocket_graphs=[all_graphs[i] for i in indices],
            labels=[all_labels[i] for i in indices],
            metadata=[all_metadata[i] for i in indices],
        )

    datasets = {
        "train": make_dataset(train_indices),
        "val": make_dataset(val_indices),
        "test": make_dataset(test_indices),
    }

    logger.info(f"Created synthetic training data:")
    logger.info(f"  Train: {len(train_indices)} samples (families: {train_families})")
    logger.info(f"  Val: {len(val_indices)} samples")
    logger.info(f"  Test: {len(test_indices)} samples (families: {test_families})")
    logger.info(f"  Label distribution: {sum(all_labels):.0f} positive / {len(all_labels) - sum(all_labels):.0f} negative")

    return datasets

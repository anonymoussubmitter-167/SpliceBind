"""Module 3: Druggability Predictor (GNN)

Graph Neural Network for predicting pocket druggability.
Uses edge-conditioned message passing on pocket residue graphs
with MC Dropout for uncertainty quantification.
"""

import logging
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Dict, List, Optional, Tuple
from torch_geometric.nn import (
    GCNConv, GATConv, global_mean_pool, global_max_pool, global_add_pool
)
from torch_geometric.data import Data, Batch

logger = logging.getLogger(__name__)


class EdgeConvLayer(nn.Module):
    """Edge-conditioned convolution layer for pocket graphs."""

    def __init__(self, in_channels: int, out_channels: int, edge_dim: int = 2):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(2 * in_channels + edge_dim, out_channels),
            nn.BatchNorm1d(out_channels),
            nn.ReLU(),
            nn.Linear(out_channels, out_channels),
        )
        self.residual = nn.Linear(in_channels, out_channels) if in_channels != out_channels else nn.Identity()

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor,
                edge_attr: torch.Tensor) -> torch.Tensor:
        row, col = edge_index
        # Aggregate neighbor features with edge features
        x_i = x[row]
        x_j = x[col]

        # Handle case where edge_attr might be empty
        if edge_attr.shape[0] == 0:
            return self.residual(x)

        edge_features = torch.cat([x_i, x_j, edge_attr], dim=-1)
        messages = self.mlp(edge_features)

        # Aggregate by mean
        out = torch.zeros_like(x[:, :messages.shape[-1]])
        if out.shape[-1] != messages.shape[-1]:
            out = torch.zeros(x.shape[0], messages.shape[-1], device=x.device)

        out = out.scatter_reduce(0, row.unsqueeze(-1).expand_as(messages), messages, reduce='mean')

        return out + self.residual(x)


class PocketGNN(nn.Module):
    """GNN for pocket druggability prediction.

    Architecture:
    - Stack of EdgeConv layers
    - Global pooling (mean + max concatenation)
    - MLP classifier with MC Dropout
    """

    def __init__(self, node_feature_dim: int = 24, edge_feature_dim: int = 2,
                 hidden_dim: int = 256, num_layers: int = 3,
                 dropout: float = 0.2, pooling: str = "mean_max",
                 num_classes: int = 1):
        super().__init__()

        self.node_feature_dim = node_feature_dim
        self.hidden_dim = hidden_dim
        self.dropout = dropout
        self.pooling = pooling

        # Input projection
        self.input_proj = nn.Sequential(
            nn.Linear(node_feature_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
        )

        # GNN layers
        self.gnn_layers = nn.ModuleList()
        self.gnn_norms = nn.ModuleList()
        for i in range(num_layers):
            self.gnn_layers.append(EdgeConvLayer(hidden_dim, hidden_dim, edge_feature_dim))
            self.gnn_norms.append(nn.BatchNorm1d(hidden_dim))

        # Pooling
        pool_dim = hidden_dim * 2 if pooling == "mean_max" else hidden_dim

        # MLP classifier
        self.classifier = nn.Sequential(
            nn.Linear(pool_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, num_classes),
        )

        # Pocket embedding output (for multi-task learning)
        self.pocket_embedding = nn.Linear(pool_dim, 128)

    def forward(self, data: Data) -> Dict[str, torch.Tensor]:
        """Forward pass.

        Args:
            data: PyG Data object with node_features (x), edge_index, edge_attr, batch

        Returns:
            Dict with:
            - druggability_score: (batch_size, 1) druggability prediction
            - pocket_embedding: (batch_size, 128) pocket representation
            - node_embeddings: (num_nodes, hidden_dim) per-node embeddings
        """
        x = data.x
        edge_index = data.edge_index
        edge_attr = data.edge_attr if data.edge_attr is not None else torch.zeros(edge_index.shape[1], 2, device=x.device)
        batch = data.batch if hasattr(data, 'batch') and data.batch is not None else torch.zeros(x.shape[0], dtype=torch.long, device=x.device)

        # Input projection
        x = self.input_proj(x)

        # GNN layers with residual connections
        for layer, norm in zip(self.gnn_layers, self.gnn_norms):
            x_res = x
            x = layer(x, edge_index, edge_attr)
            x = norm(x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
            x = x + x_res  # Residual

        node_embeddings = x

        # Global pooling
        if self.pooling == "mean_max":
            x_mean = global_mean_pool(x, batch)
            x_max = global_max_pool(x, batch)
            graph_embedding = torch.cat([x_mean, x_max], dim=-1)
        elif self.pooling == "mean":
            graph_embedding = global_mean_pool(x, batch)
        elif self.pooling == "max":
            graph_embedding = global_max_pool(x, batch)
        else:
            graph_embedding = global_mean_pool(x, batch)

        # Druggability prediction (logits - apply sigmoid externally for inference)
        druggability_logits = self.classifier(graph_embedding)
        druggability_score = torch.sigmoid(druggability_logits)

        # Pocket embedding for multi-task learning
        pocket_emb = self.pocket_embedding(graph_embedding)

        return {
            "druggability_logits": druggability_logits,
            "druggability_score": druggability_score,
            "pocket_embedding": pocket_emb,
            "node_embeddings": node_embeddings,
        }


class MCDropoutPredictor:
    """Monte Carlo Dropout for uncertainty estimation."""

    def __init__(self, model: PocketGNN, num_samples: int = 20):
        self.model = model
        self.num_samples = num_samples

    @torch.no_grad()
    def predict_with_uncertainty(self, data: Data) -> Dict[str, float]:
        """Run multiple forward passes with dropout for uncertainty estimation."""
        self.model.train()  # Keep dropout active

        predictions = []
        embeddings = []

        for _ in range(self.num_samples):
            output = self.model(data)
            predictions.append(output["druggability_score"].cpu().numpy())
            embeddings.append(output["pocket_embedding"].cpu().numpy())

        predictions = np.array(predictions).squeeze()
        embeddings = np.array(embeddings)

        # Statistics
        mean_pred = predictions.mean(axis=0)
        std_pred = predictions.std(axis=0)
        mean_embedding = embeddings.mean(axis=0)

        # Ensure scalar outputs for single samples
        if mean_pred.ndim == 0:
            mean_pred = float(mean_pred)
            std_pred = float(std_pred)

        return {
            "mean_prediction": mean_pred,
            "std_prediction": std_pred,
            "confidence": 1.0 - float(np.mean(std_pred)),  # Higher confidence = lower std
            "all_predictions": predictions.tolist(),
            "mean_embedding": mean_embedding.squeeze().tolist(),
        }


class PocketSimilarityNetwork(nn.Module):
    """Network for predicting similarity between two pocket embeddings.

    Used for cross-isoform pocket comparison (Task 2 in multi-task learning).
    """

    def __init__(self, embedding_dim: int = 128, hidden_dim: int = 64):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(embedding_dim * 2, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1),
            nn.Sigmoid(),
        )

    def forward(self, emb1: torch.Tensor, emb2: torch.Tensor) -> torch.Tensor:
        """Predict similarity score between two pocket embeddings."""
        combined = torch.cat([emb1, emb2], dim=-1)
        return self.network(combined)


class DeltaDruggabilityPredictor(nn.Module):
    """Predict change in druggability between canonical and isoform pockets.

    Takes pocket embeddings from both versions and predicts ΔDruggability.
    """

    def __init__(self, embedding_dim: int = 128, hidden_dim: int = 64):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(embedding_dim * 2 + 3, hidden_dim),  # +3 for structural features
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, 1),
            nn.Tanh(),  # ΔDruggability in [-1, 1]
        )

    def forward(self, canonical_emb: torch.Tensor, isoform_emb: torch.Tensor,
                structural_features: torch.Tensor) -> torch.Tensor:
        """Predict ΔDruggability.

        Args:
            canonical_emb: (batch, 128) canonical pocket embedding
            isoform_emb: (batch, 128) isoform pocket embedding
            structural_features: (batch, 3) [length_diff, rmsd, plddt_drop]
        """
        combined = torch.cat([canonical_emb, isoform_emb, structural_features], dim=-1)
        return self.network(combined)


def create_pocket_data(pocket_graph: Dict, label: Optional[float] = None) -> Data:
    """Convert pocket graph dict to PyG Data object."""
    data = Data(
        x=pocket_graph["node_features"],
        edge_index=pocket_graph["edge_index"],
        edge_attr=pocket_graph["edge_attr"],
        num_nodes=pocket_graph["num_nodes"],
    )

    if label is not None:
        data.y = torch.tensor([label], dtype=torch.float32)

    return data

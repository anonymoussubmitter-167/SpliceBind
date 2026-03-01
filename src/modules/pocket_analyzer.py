"""Module 2: Binding Pocket Analyzer

Detects and analyzes binding pockets on protein structures,
computes pocket descriptors, and performs cross-isoform pocket
alignment and conservation scoring.
"""

import os
import logging
import subprocess
import tempfile
import numpy as np
import torch
from typing import Dict, List, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

# Amino acid properties for pocket descriptors
AA_HYDROPHOBICITY = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2,
}

AA_CHARGE = {
    'ALA': 0, 'ARG': 1, 'ASN': 0, 'ASP': -1, 'CYS': 0,
    'GLN': 0, 'GLU': -1, 'GLY': 0, 'HIS': 0.5, 'ILE': 0,
    'LEU': 0, 'LYS': 1, 'MET': 0, 'PHE': 0, 'PRO': 0,
    'SER': 0, 'THR': 0, 'TRP': 0, 'TYR': 0, 'VAL': 0,
}

AA_POLARITY = {
    'ALA': 0, 'ARG': 1, 'ASN': 1, 'ASP': 1, 'CYS': 0,
    'GLN': 1, 'GLU': 1, 'GLY': 0, 'HIS': 1, 'ILE': 0,
    'LEU': 0, 'LYS': 1, 'MET': 0, 'PHE': 0, 'PRO': 0,
    'SER': 1, 'THR': 1, 'TRP': 0, 'TYR': 1, 'VAL': 0,
}

# One-hot encoding for amino acids
AA_LIST = sorted(AA_HYDROPHOBICITY.keys())
AA_TO_IDX = {aa: i for i, aa in enumerate(AA_LIST)}


class GeometricPocketDetector:
    """Detect binding pockets using geometric analysis of protein structure.

    Uses a simplified approach based on:
    1. Identifying surface residues
    2. Finding concavities (potential pockets)
    3. Scoring pocket druggability based on composition
    """

    def __init__(self, probe_radius: float = 1.4, min_pocket_size: int = 5,
                 contact_distance: float = 6.0):
        self.probe_radius = probe_radius
        self.min_pocket_size = min_pocket_size
        self.contact_distance = contact_distance

    def detect_pockets(self, pdb_path: str) -> List[Dict]:
        """Detect binding pockets from PDB structure.

        Returns list of pocket dicts, each containing:
        - residue_indices: indices of residues in pocket
        - residue_names: amino acid names
        - center: pocket center coordinates
        - volume_estimate: approximate pocket volume
        - descriptors: computed pocket descriptors
        """
        # Parse structure
        coords, residues, bfactors = self._parse_pdb(pdb_path)
        if len(coords) < 10:
            logger.warning(f"Too few residues ({len(coords)}) to detect pockets")
            return []

        # Compute distance matrix
        dist_matrix = self._compute_distances(coords)

        # Identify surface residues (simplified: fewer neighbors within cutoff)
        surface_mask = self._identify_surface(dist_matrix, cutoff=10.0, max_neighbors=15)

        # Find pocket clusters (concavities among surface residues)
        pockets = self._find_concavities(coords, residues, bfactors, dist_matrix,
                                         surface_mask)

        # Score and rank pockets
        for pocket in pockets:
            pocket["descriptors"] = self._compute_pocket_descriptors(pocket, coords)
            pocket["druggability_score"] = self._estimate_druggability(pocket)

        # Sort by druggability score
        pockets.sort(key=lambda p: p["druggability_score"], reverse=True)

        return pockets

    def _parse_pdb(self, pdb_path: str) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Parse PDB file and extract CA coordinates, residue names, B-factors."""
        coords = []
        residues = []
        bfactors = []

        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    coords.append([
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    ])
                    residues.append(line[17:20].strip())
                    bfactors.append(float(line[60:66]))

        return np.array(coords), residues, np.array(bfactors)

    def _compute_distances(self, coords: np.ndarray) -> np.ndarray:
        """Compute pairwise distance matrix."""
        diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
        return np.sqrt(np.sum(diff ** 2, axis=-1))

    def _identify_surface(self, dist_matrix: np.ndarray, cutoff: float = 10.0,
                          max_neighbors: int = 15) -> np.ndarray:
        """Identify surface residues by counting neighbors."""
        n = len(dist_matrix)
        neighbor_counts = np.sum(dist_matrix < cutoff, axis=1) - 1  # subtract self
        # Surface residues have fewer close neighbors
        threshold = np.percentile(neighbor_counts, 40)
        return neighbor_counts <= threshold

    def _find_concavities(self, coords: np.ndarray, residues: List[str],
                          bfactors: np.ndarray, dist_matrix: np.ndarray,
                          surface_mask: np.ndarray) -> List[Dict]:
        """Find pocket-like concavities among surface residues."""
        n = len(coords)
        surface_indices = np.where(surface_mask)[0]

        if len(surface_indices) < self.min_pocket_size:
            return []

        # Cluster surface residues by spatial proximity
        visited = set()
        pockets = []

        for seed in surface_indices:
            if seed in visited:
                continue

            # BFS to find cluster
            cluster = []
            queue = [seed]
            while queue:
                curr = queue.pop(0)
                if curr in visited:
                    continue
                visited.add(curr)

                if surface_mask[curr]:
                    cluster.append(curr)
                    # Add neighbors within contact distance
                    neighbors = np.where(dist_matrix[curr] < self.contact_distance)[0]
                    for nb in neighbors:
                        if nb not in visited and surface_mask[nb]:
                            queue.append(nb)

            if len(cluster) >= self.min_pocket_size:
                cluster = np.array(cluster)
                pocket_coords = coords[cluster]
                center = pocket_coords.mean(axis=0)

                # Estimate volume (convex hull approximation)
                ranges = pocket_coords.max(axis=0) - pocket_coords.min(axis=0)
                volume_est = np.prod(ranges) * 0.5  # rough estimate

                pocket_residues = [residues[i] for i in cluster]
                pocket_bfactors = bfactors[cluster]

                pockets.append({
                    "residue_indices": cluster.tolist(),
                    "residue_names": pocket_residues,
                    "center": center.tolist(),
                    "volume_estimate": float(volume_est),
                    "num_residues": len(cluster),
                    "mean_plddt": float(pocket_bfactors.mean()),
                    "min_plddt": float(pocket_bfactors.min()),
                })

        return pockets

    def _compute_pocket_descriptors(self, pocket: Dict, all_coords: np.ndarray) -> Dict:
        """Compute physicochemical descriptors for a pocket."""
        residues = pocket["residue_names"]
        n = len(residues)

        if n == 0:
            return {}

        # Hydrophobicity
        hydro_values = [AA_HYDROPHOBICITY.get(r, 0.0) for r in residues]
        mean_hydro = np.mean(hydro_values)
        hydro_fraction = sum(1 for h in hydro_values if h > 0) / n

        # Charge
        charges = [AA_CHARGE.get(r, 0.0) for r in residues]
        net_charge = sum(charges)
        charged_fraction = sum(1 for c in charges if abs(c) > 0) / n

        # Polarity
        polarity = [AA_POLARITY.get(r, 0.0) for r in residues]
        polar_fraction = sum(polarity) / n

        # Amino acid composition vector
        aa_counts = np.zeros(20)
        for r in residues:
            if r in AA_TO_IDX:
                aa_counts[AA_TO_IDX[r]] += 1
        aa_composition = aa_counts / max(n, 1)

        # Pocket geometry
        indices = pocket["residue_indices"]
        pocket_coords = all_coords[indices]
        if len(pocket_coords) > 1:
            centroid = pocket_coords.mean(axis=0)
            distances_to_center = np.linalg.norm(pocket_coords - centroid, axis=1)
            radius = distances_to_center.max()
            depth = distances_to_center.mean()
        else:
            radius = 0.0
            depth = 0.0

        return {
            "hydrophobicity_mean": float(mean_hydro),
            "hydrophobic_fraction": float(hydro_fraction),
            "net_charge": float(net_charge),
            "charged_fraction": float(charged_fraction),
            "polar_fraction": float(polar_fraction),
            "aa_composition": aa_composition.tolist(),
            "radius": float(radius),
            "depth": float(depth),
            "volume": pocket["volume_estimate"],
        }

    def _estimate_druggability(self, pocket: Dict) -> float:
        """Estimate pocket druggability based on descriptors.

        Uses a simplified scoring function based on known druggability features:
        - Moderate hydrophobicity (not too high, not too low)
        - Appropriate size (500-2000 A³)
        - Mixed polar/hydrophobic character
        - Sufficient depth
        """
        desc = pocket.get("descriptors", {})
        if not desc:
            return 0.0

        score = 0.0

        # Volume component (optimal: 500-2000 A³)
        vol = desc.get("volume", 0)
        if 200 < vol < 5000:
            vol_score = 1.0 - abs(np.log10(vol / 1000)) / 1.5
            score += max(0, vol_score) * 0.3

        # Hydrophobicity component (moderate is best)
        hydro = desc.get("hydrophobicity_mean", 0)
        hydro_score = 1.0 - abs(hydro - 0.5) / 3.0
        score += max(0, hydro_score) * 0.2

        # Size component
        n_res = pocket.get("num_residues", 0)
        if 8 <= n_res <= 30:
            score += 0.2
        elif 5 <= n_res < 8 or 30 < n_res <= 50:
            score += 0.1

        # Mixed character (both polar and hydrophobic)
        polar = desc.get("polar_fraction", 0)
        if 0.3 <= polar <= 0.7:
            score += 0.15
        elif 0.2 <= polar <= 0.8:
            score += 0.07

        # Depth component
        depth = desc.get("depth", 0)
        if depth > 3.0:
            score += 0.15

        return min(1.0, max(0.0, score))


class PocketGraphBuilder:
    """Build graph representation of binding pockets for GNN input."""

    def __init__(self, contact_distance: float = 6.0,
                 node_feature_dim: int = 30):
        self.contact_distance = contact_distance
        self.node_feature_dim = node_feature_dim

    def build_pocket_graph(self, pocket: Dict, all_coords: np.ndarray,
                          esm_embeddings: Optional[np.ndarray] = None) -> Dict:
        """Build a graph representation of a pocket.

        Returns dict with:
        - node_features: (N, feature_dim) tensor
        - edge_index: (2, E) tensor of edges
        - edge_features: (E, edge_dim) tensor
        - pocket_label: druggability label (if available)
        """
        indices = pocket["residue_indices"]
        residues = pocket["residue_names"]
        n = len(indices)

        if n == 0:
            return None

        # Node features
        node_features = []
        for i, (idx, res) in enumerate(zip(indices, residues)):
            feat = self._residue_features(res)

            # Add pLDDT as feature
            plddt = pocket.get("mean_plddt", 70.0)
            feat.append(plddt / 100.0)

            # Add ESM embedding if available
            if esm_embeddings is not None and idx < len(esm_embeddings):
                # Use first 32 dims of ESM embedding as additional features
                esm_feat = esm_embeddings[idx, :32]
                feat.extend(esm_feat.tolist())

            node_features.append(feat)

        node_features = np.array(node_features, dtype=np.float32)

        # Edge construction: connect residues within contact distance
        pocket_coords = all_coords[indices]
        edge_index = []
        edge_features = []

        for i in range(n):
            for j in range(i + 1, n):
                dist = np.linalg.norm(pocket_coords[i] - pocket_coords[j])
                if dist < self.contact_distance:
                    edge_index.append([i, j])
                    edge_index.append([j, i])  # undirected

                    # Edge features: distance, sequence separation
                    seq_sep = abs(indices[i] - indices[j])
                    edge_feat = [dist / self.contact_distance, min(seq_sep / 50.0, 1.0)]
                    edge_features.append(edge_feat)
                    edge_features.append(edge_feat)

        if not edge_index:
            # Fully connect if no edges within cutoff
            for i in range(n):
                for j in range(n):
                    if i != j:
                        dist = np.linalg.norm(pocket_coords[i] - pocket_coords[j])
                        edge_index.append([i, j])
                        seq_sep = abs(indices[i] - indices[j])
                        edge_features.append([dist / 20.0, min(seq_sep / 50.0, 1.0)])

        edge_index = np.array(edge_index, dtype=np.int64).T if edge_index else np.zeros((2, 0), dtype=np.int64)
        edge_features = np.array(edge_features, dtype=np.float32) if edge_features else np.zeros((0, 2), dtype=np.float32)

        return {
            "node_features": torch.tensor(node_features),
            "edge_index": torch.tensor(edge_index),
            "edge_attr": torch.tensor(edge_features),
            "num_nodes": n,
            "pocket_center": pocket.get("center", [0, 0, 0]),
            "druggability_score": pocket.get("druggability_score", 0.0),
        }

    def _residue_features(self, residue_name: str) -> List[float]:
        """Compute feature vector for a residue."""
        features = []

        # One-hot amino acid type (20 dims)
        one_hot = [0.0] * 20
        if residue_name in AA_TO_IDX:
            one_hot[AA_TO_IDX[residue_name]] = 1.0
        features.extend(one_hot)

        # Physicochemical properties
        features.append(AA_HYDROPHOBICITY.get(residue_name, 0.0) / 4.5)  # normalized
        features.append(AA_CHARGE.get(residue_name, 0.0))
        features.append(AA_POLARITY.get(residue_name, 0.0))

        return features


class CrossIsoformPocketAligner:
    """Align and compare pockets between canonical and isoform structures."""

    def __init__(self, similarity_threshold: float = 0.5):
        self.similarity_threshold = similarity_threshold

    def align_pockets(self, canonical_pockets: List[Dict],
                     isoform_pockets: List[Dict]) -> List[Dict]:
        """Find corresponding pockets between canonical and isoform structures.

        Returns list of pocket alignment results.
        """
        alignments = []

        for i, can_pocket in enumerate(canonical_pockets):
            best_match = None
            best_similarity = -1

            for j, iso_pocket in enumerate(isoform_pockets):
                sim = self._compute_pocket_similarity(can_pocket, iso_pocket)
                if sim > best_similarity:
                    best_similarity = sim
                    best_match = j

            alignments.append({
                "canonical_pocket_idx": i,
                "isoform_pocket_idx": best_match,
                "similarity": best_similarity,
                "pocket_preserved": best_similarity >= self.similarity_threshold,
                "canonical_druggability": can_pocket.get("druggability_score", 0),
                "isoform_druggability": (
                    isoform_pockets[best_match].get("druggability_score", 0)
                    if best_match is not None else 0
                ),
                "delta_druggability": (
                    isoform_pockets[best_match].get("druggability_score", 0) -
                    can_pocket.get("druggability_score", 0)
                    if best_match is not None else -can_pocket.get("druggability_score", 0)
                ),
            })

        return alignments

    def _compute_pocket_similarity(self, pocket1: Dict, pocket2: Dict) -> float:
        """Compute similarity between two pockets based on descriptors."""
        desc1 = pocket1.get("descriptors", {})
        desc2 = pocket2.get("descriptors", {})

        if not desc1 or not desc2:
            return 0.0

        sim = 0.0
        n_features = 0

        # Amino acid composition similarity (cosine)
        comp1 = np.array(desc1.get("aa_composition", [0] * 20))
        comp2 = np.array(desc2.get("aa_composition", [0] * 20))
        norm1 = np.linalg.norm(comp1)
        norm2 = np.linalg.norm(comp2)
        if norm1 > 0 and norm2 > 0:
            cos_sim = np.dot(comp1, comp2) / (norm1 * norm2)
            sim += cos_sim
            n_features += 1

        # Volume similarity
        v1 = desc1.get("volume", 0)
        v2 = desc2.get("volume", 0)
        if v1 > 0 and v2 > 0:
            vol_sim = min(v1, v2) / max(v1, v2)
            sim += vol_sim
            n_features += 1

        # Hydrophobicity similarity
        h1 = desc1.get("hydrophobicity_mean", 0)
        h2 = desc2.get("hydrophobicity_mean", 0)
        h_sim = 1.0 - abs(h1 - h2) / 9.0  # range is -4.5 to 4.5
        sim += max(0, h_sim)
        n_features += 1

        # Size similarity (number of residues)
        n1 = pocket1.get("num_residues", 0)
        n2 = pocket2.get("num_residues", 0)
        if n1 > 0 and n2 > 0:
            size_sim = min(n1, n2) / max(n1, n2)
            sim += size_sim
            n_features += 1

        # Spatial overlap (center distance)
        c1 = np.array(pocket1.get("center", [0, 0, 0]))
        c2 = np.array(pocket2.get("center", [0, 0, 0]))
        center_dist = np.linalg.norm(c1 - c2)
        spatial_sim = max(0, 1.0 - center_dist / 20.0)
        sim += spatial_sim
        n_features += 1

        return sim / max(n_features, 1)

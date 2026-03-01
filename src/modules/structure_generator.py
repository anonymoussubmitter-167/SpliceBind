"""Module 1: Isoform Structure Generator

Generates protein structure representations using ESM-2 embeddings
and AlphaFold predicted structures. Tracks structural uncertainty
in spliced regions.
"""

import os
import logging
import numpy as np
import torch
import torch.nn as nn
from typing import Dict, List, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)


class ESMEmbedder:
    """Extract per-residue embeddings from ESM-2 model."""

    def __init__(self, model_name: str = "esm2_t33_650M_UR50D",
                 device: str = "cuda", cache_dir: str = None):
        self.device = torch.device(device if torch.cuda.is_available() else "cpu")
        self.model_name = model_name
        self.cache_dir = cache_dir
        self.model = None
        self.alphabet = None
        self.batch_converter = None

    def load_model(self):
        """Lazy-load ESM model."""
        if self.model is not None:
            return

        import esm
        logger.info(f"Loading ESM model: {self.model_name}")
        self.model, self.alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.batch_converter = self.alphabet.get_batch_converter()
        self.model = self.model.to(self.device)
        self.model.eval()
        logger.info(f"ESM model loaded on {self.device}")

    @torch.no_grad()
    def embed_sequence(self, sequence: str, label: str = "protein") -> np.ndarray:
        """Get per-residue embeddings for a protein sequence.

        Returns: numpy array of shape (seq_len, embed_dim)
        """
        self.load_model()

        # Check cache
        if self.cache_dir:
            cache_path = os.path.join(self.cache_dir, f"{label}_esm.npy")
            if os.path.exists(cache_path):
                return np.load(cache_path)

        # Truncate long sequences (ESM-2 max is 1022 tokens)
        max_len = 1022
        if len(sequence) > max_len:
            logger.warning(f"Sequence {label} truncated from {len(sequence)} to {max_len}")
            sequence = sequence[:max_len]

        data = [(label, sequence)]
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        batch_tokens = batch_tokens.to(self.device)

        results = self.model(batch_tokens, repr_layers=[33], return_contacts=False)

        # Shape: (1, seq_len+2, embed_dim) -> remove BOS/EOS tokens
        embeddings = results["representations"][33][0, 1:-1, :].cpu().numpy()

        # Cache
        if self.cache_dir:
            os.makedirs(self.cache_dir, exist_ok=True)
            cache_path = os.path.join(self.cache_dir, f"{label}_esm.npy")
            np.save(cache_path, embeddings)

        return embeddings

    @torch.no_grad()
    def embed_batch(self, sequences: List[Tuple[str, str]]) -> Dict[str, np.ndarray]:
        """Embed multiple sequences. Input: list of (label, sequence) tuples."""
        self.load_model()

        results_dict = {}
        # Process in small batches to fit in GPU memory
        batch_size = 4
        for i in range(0, len(sequences), batch_size):
            batch = sequences[i:i + batch_size]

            # Truncate
            batch = [(label, seq[:1022]) for label, seq in batch]

            batch_labels, batch_strs, batch_tokens = self.batch_converter(batch)
            batch_tokens = batch_tokens.to(self.device)

            results = self.model(batch_tokens, repr_layers=[33], return_contacts=False)
            representations = results["representations"][33]

            for j, (label, seq) in enumerate(batch):
                seq_len = len(seq)
                emb = representations[j, 1:seq_len + 1, :].cpu().numpy()
                results_dict[label] = emb

        return results_dict


class StructureFeatureExtractor:
    """Extract structural features from PDB files for pocket analysis."""

    def __init__(self):
        pass

    def extract_features(self, pdb_path: str) -> Dict:
        """Extract structural features from a PDB file.

        Returns dict with coordinates, contacts, secondary structure info, etc.
        """
        features = {
            "ca_coords": [],
            "residue_names": [],
            "residue_numbers": [],
            "bfactors": [],  # pLDDT in AlphaFold structures
        }

        try:
            with open(pdb_path) as f:
                for line in f:
                    if line.startswith("ATOM") and line[12:16].strip() == "CA":
                        features["ca_coords"].append([
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ])
                        features["residue_names"].append(line[17:20].strip())
                        features["residue_numbers"].append(int(line[22:26]))
                        features["bfactors"].append(float(line[60:66]))
        except Exception as e:
            logger.error(f"Error parsing PDB {pdb_path}: {e}")
            return features

        features["ca_coords"] = np.array(features["ca_coords"])
        features["bfactors"] = np.array(features["bfactors"])
        features["num_residues"] = len(features["residue_names"])

        # Compute contact map
        if len(features["ca_coords"]) > 0:
            coords = features["ca_coords"]
            diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
            dist_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))
            features["distance_matrix"] = dist_matrix
            features["contact_map_8A"] = (dist_matrix < 8.0).astype(np.float32)
            features["contact_map_12A"] = (dist_matrix < 12.0).astype(np.float32)

        return features

    def compute_splice_region_confidence(self, features: Dict,
                                          splice_positions: List[int],
                                          window: int = 10) -> Dict:
        """Assess structural confidence around splice junctions.

        splice_positions: list of residue positions at exon boundaries
        window: number of residues around junction to evaluate
        """
        if not splice_positions or len(features.get("bfactors", [])) == 0:
            return {"splice_confidence": None}

        bfactors = features["bfactors"]
        res_nums = features["residue_numbers"]

        splice_plddts = []
        non_splice_plddts = []

        for i, (res_num, plddt) in enumerate(zip(res_nums, bfactors)):
            near_splice = any(abs(res_num - sp) <= window for sp in splice_positions)
            if near_splice:
                splice_plddts.append(plddt)
            else:
                non_splice_plddts.append(plddt)

        return {
            "splice_region_mean_plddt": np.mean(splice_plddts) if splice_plddts else None,
            "non_splice_region_mean_plddt": np.mean(non_splice_plddts) if non_splice_plddts else None,
            "plddt_drop_at_splice": (
                np.mean(non_splice_plddts) - np.mean(splice_plddts)
                if splice_plddts and non_splice_plddts else None
            ),
            "num_low_confidence_splice_residues": sum(1 for p in splice_plddts if p < 70),
        }


class IsoformStructureComparator:
    """Compare structures between canonical and isoform proteins."""

    def __init__(self):
        self.feature_extractor = StructureFeatureExtractor()

    def compare_structures(self, canonical_pdb: str, isoform_pdb: str,
                          splice_positions: Optional[List[int]] = None) -> Dict:
        """Compare two protein structures and identify differences.

        Returns structural comparison metrics.
        """
        feat_can = self.feature_extractor.extract_features(canonical_pdb)
        feat_iso = self.feature_extractor.extract_features(isoform_pdb)

        if feat_can["num_residues"] == 0 or feat_iso["num_residues"] == 0:
            return {"error": "Could not extract features from one or both structures"}

        comparison = {
            "canonical_length": feat_can["num_residues"],
            "isoform_length": feat_iso["num_residues"],
            "length_difference": feat_iso["num_residues"] - feat_can["num_residues"],
            "canonical_mean_plddt": np.mean(feat_can["bfactors"]),
            "isoform_mean_plddt": np.mean(feat_iso["bfactors"]),
        }

        # Compute RMSD for aligned regions
        min_len = min(feat_can["num_residues"], feat_iso["num_residues"])
        coords_can = feat_can["ca_coords"][:min_len]
        coords_iso = feat_iso["ca_coords"][:min_len]

        from ..data_pipeline.fetch_structures import superimpose_structures
        _, rmsd = superimpose_structures(coords_can, coords_iso)
        comparison["global_rmsd"] = rmsd

        # Splice region analysis
        if splice_positions:
            comparison["splice_confidence_canonical"] = (
                self.feature_extractor.compute_splice_region_confidence(feat_can, splice_positions)
            )
            comparison["splice_confidence_isoform"] = (
                self.feature_extractor.compute_splice_region_confidence(feat_iso, splice_positions)
            )

        return comparison

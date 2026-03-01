"""ESM-2 protein language model embeddings for pocket node features.

Extracts per-residue embeddings from ESM-2 (esm2_t33_650M_UR50D)
to enrich pocket graph node features with evolutionary context.
"""

import os
import logging
import numpy as np
import torch
from typing import Dict, List, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)


class ESMEmbeddingExtractor:
    """Extract per-residue embeddings from ESM-2.

    Caches embeddings to disk for reuse.
    """

    def __init__(self, model_name: str = "esm2_t33_650M_UR50D",
                 cache_dir: Optional[str] = None,
                 device: Optional[str] = None):
        self.model_name = model_name
        self.cache_dir = cache_dir
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.model = None
        self.alphabet = None
        self.batch_converter = None
        self.repr_layer = 33  # Last layer for t33 model
        self.embedding_dim = 1280

        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)

    def _load_model(self):
        """Lazy-load ESM-2 model."""
        if self.model is not None:
            return

        logger.info(f"Loading ESM-2 model: {self.model_name}...")
        try:
            import esm
            self.model, self.alphabet = esm.pretrained.esm2_t33_650M_UR50D()
            self.batch_converter = self.alphabet.get_batch_converter()
            self.model = self.model.to(self.device)
            self.model.eval()
            logger.info(f"ESM-2 loaded on {self.device}")
        except Exception as e:
            logger.warning(f"Failed to load ESM-2 on {self.device}: {e}")
            if self.device != "cpu":
                logger.info("Retrying ESM-2 on CPU...")
                try:
                    import esm
                    if self.model is None:
                        self.model, self.alphabet = esm.pretrained.esm2_t33_650M_UR50D()
                        self.batch_converter = self.alphabet.get_batch_converter()
                    self.device = "cpu"
                    self.model = self.model.to("cpu")
                    self.model.eval()
                    torch.cuda.empty_cache()
                    logger.info("ESM-2 loaded on CPU (fallback)")
                except Exception as e2:
                    logger.warning(f"Failed to load ESM-2 on CPU: {e2}. Using random embeddings.")
                    self.model = None
            else:
                logger.warning("Using random embeddings.")
                self.model = None

    def extract_embeddings(self, sequence: str, protein_id: str = "") -> np.ndarray:
        """Extract per-residue embeddings for a protein sequence.

        Args:
            sequence: Amino acid sequence (single letter)
            protein_id: Optional ID for caching

        Returns:
            np.ndarray of shape (L, 1280) where L is sequence length
        """
        # Check cache
        if self.cache_dir and protein_id:
            cache_path = os.path.join(self.cache_dir, f"{protein_id}_esm2.npy")
            if os.path.exists(cache_path):
                return np.load(cache_path)

        self._load_model()

        if self.model is None:
            # Fallback: random embeddings (for when ESM-2 can't load)
            embeddings = np.random.randn(len(sequence), self.embedding_dim).astype(np.float32) * 0.01
        else:
            # Truncate if too long for ESM-2 (max ~1024 tokens recommended)
            max_len = 1022  # ESM-2 adds BOS/EOS tokens
            if len(sequence) > max_len:
                logger.warning(f"Sequence {protein_id} length {len(sequence)} > {max_len}, "
                              f"processing in chunks")
                embeddings = self._extract_chunked(sequence, max_len)
            else:
                embeddings = self._extract_single(sequence)

        # Cache
        if self.cache_dir and protein_id:
            cache_path = os.path.join(self.cache_dir, f"{protein_id}_esm2.npy")
            np.save(cache_path, embeddings)

        return embeddings

    def _extract_single(self, sequence: str) -> np.ndarray:
        """Extract embeddings for a single sequence."""
        data = [("protein", sequence)]
        _, _, tokens = self.batch_converter(data)
        tokens = tokens.to(self.device)

        with torch.no_grad():
            results = self.model(tokens, repr_layers=[self.repr_layer])

        # Remove BOS/EOS tokens: shape (1, L+2, 1280) -> (L, 1280)
        embeddings = results["representations"][self.repr_layer][0, 1:-1, :].cpu().numpy()
        return embeddings.astype(np.float32)

    def _extract_chunked(self, sequence: str, chunk_size: int = 1022,
                          overlap: int = 100) -> np.ndarray:
        """Extract embeddings in overlapping chunks for long sequences."""
        seq_len = len(sequence)
        embeddings = np.zeros((seq_len, self.embedding_dim), dtype=np.float32)
        counts = np.zeros(seq_len, dtype=np.float32)

        start = 0
        while start < seq_len:
            end = min(start + chunk_size, seq_len)
            chunk_seq = sequence[start:end]

            chunk_emb = self._extract_single(chunk_seq)

            embeddings[start:start + len(chunk_emb)] += chunk_emb
            counts[start:start + len(chunk_emb)] += 1

            if end >= seq_len:
                break
            start += chunk_size - overlap

        # Average overlapping regions
        mask = counts > 0
        embeddings[mask] /= counts[mask, np.newaxis]

        return embeddings

    def extract_for_pdb(self, pdb_path: str, protein_id: str = "") -> np.ndarray:
        """Extract embeddings for a PDB file by reading the sequence."""
        sequence = self._pdb_to_sequence(pdb_path)
        if not sequence:
            logger.warning(f"Could not extract sequence from {pdb_path}")
            return np.zeros((0, self.embedding_dim), dtype=np.float32)

        if not protein_id:
            protein_id = os.path.basename(pdb_path).replace(".pdb", "")

        return self.extract_embeddings(sequence, protein_id)

    def _pdb_to_sequence(self, pdb_path: str) -> str:
        """Extract single-letter amino acid sequence from PDB file."""
        three_to_one = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        }

        residues = []
        seen_res_ids = set()

        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    res_name = line[17:20].strip()
                    chain = line[21]
                    res_id = line[22:26].strip()
                    key = f"{chain}_{res_id}"

                    if key not in seen_res_ids:
                        seen_res_ids.add(key)
                        aa = three_to_one.get(res_name, 'X')
                        residues.append(aa)

        return "".join(residues)


def extract_all_embeddings(structures_dir: str, embeddings_dir: str,
                           device: str = "cpu") -> Dict[str, np.ndarray]:
    """Extract ESM-2 embeddings for all PDB structures.

    Returns dict mapping protein_id -> embeddings array.
    """
    extractor = ESMEmbeddingExtractor(cache_dir=embeddings_dir, device=device)
    results = {}

    pdb_files = [f for f in os.listdir(structures_dir) if f.endswith(".pdb")]
    logger.info(f"Extracting ESM-2 embeddings for {len(pdb_files)} structures...")

    for pdb_file in pdb_files:
        pdb_path = os.path.join(structures_dir, pdb_file)
        protein_id = pdb_file.replace(".pdb", "")

        try:
            embeddings = extractor.extract_for_pdb(pdb_path, protein_id)
            results[protein_id] = embeddings
            logger.info(f"  {protein_id}: {embeddings.shape}")
        except Exception as e:
            logger.warning(f"  Failed for {protein_id}: {e}")

    return results

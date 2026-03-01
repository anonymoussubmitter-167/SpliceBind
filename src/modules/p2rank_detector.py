"""P2Rank-based pocket detection.

Replaces the simplified GeometricPocketDetector with P2Rank,
a GNN-based pocket detection tool specifically benchmarked on AlphaFold structures.
Falls back to geometric detection if P2Rank is not available.
"""

import os
import re
import logging
import subprocess
import tempfile
import numpy as np
from typing import Dict, List, Optional, Tuple
from pathlib import Path

logger = logging.getLogger(__name__)

# Default P2Rank installation path
P2RANK_HOME = os.environ.get(
    "P2RANK_HOME",
    os.path.join(os.path.dirname(__file__), "..", "..", "tools", "p2rank_2.4.2")
)


class P2RankDetector:
    """Detect binding pockets using P2Rank.

    P2Rank is a machine-learning based tool for prediction of ligand-binding
    sites from protein structures. It outputs pocket residues, scores, and
    descriptors.
    """

    def __init__(self, p2rank_home: str = P2RANK_HOME, min_pocket_score: float = 0.1):
        self.p2rank_home = p2rank_home
        self.prank_cmd = os.path.join(p2rank_home, "prank")
        self.min_pocket_score = min_pocket_score

        # Verify P2Rank exists
        if not os.path.exists(self.prank_cmd):
            logger.warning(f"P2Rank not found at {self.prank_cmd}. "
                          f"Will fall back to geometric detection.")
            self.available = False
        else:
            self.available = True
            # Make executable
            os.chmod(self.prank_cmd, 0o755)

    def detect_pockets(self, pdb_path: str, output_dir: Optional[str] = None) -> List[Dict]:
        """Detect pockets using P2Rank.

        Returns list of pocket dicts compatible with PocketGraphBuilder.
        """
        if not self.available:
            logger.warning("P2Rank not available, falling back to geometric detection")
            from .pocket_analyzer import GeometricPocketDetector
            return GeometricPocketDetector().detect_pockets(pdb_path)

        # Parse PDB for coordinates (needed for graph building)
        coords, residues, bfactors = self._parse_pdb(pdb_path)

        if output_dir is None:
            output_dir = tempfile.mkdtemp(prefix="p2rank_")

        # Check if P2Rank output already exists (cache)
        pdb_name = os.path.basename(pdb_path).replace(".pdb", "")
        existing_csv = self._find_predictions_csv(output_dir)
        if existing_csv:
            logger.debug(f"Using cached P2Rank output: {existing_csv}")
            pockets = self._parse_p2rank_output(output_dir, pdb_name, coords, residues, bfactors)
            if pockets:
                pockets = [p for p in pockets if p.get("druggability_score", 0) >= self.min_pocket_score]
                pockets.sort(key=lambda p: p.get("druggability_score", 0), reverse=True)
                return pockets

        try:
            pockets = self._run_p2rank(pdb_path, output_dir, coords, residues, bfactors)
        except Exception as e:
            logger.warning(f"P2Rank failed: {e}. Falling back to geometric detection.")
            from .pocket_analyzer import GeometricPocketDetector
            return GeometricPocketDetector().detect_pockets(pdb_path)

        # Filter by score
        pockets = [p for p in pockets if p.get("druggability_score", 0) >= self.min_pocket_score]

        # Sort by score
        pockets.sort(key=lambda p: p.get("druggability_score", 0), reverse=True)

        return pockets

    def _run_p2rank(self, pdb_path: str, output_dir: str,
                     coords: np.ndarray, residues: List[str],
                     bfactors: np.ndarray) -> List[Dict]:
        """Run P2Rank on a single PDB file and parse results."""
        pdb_name = os.path.basename(pdb_path).replace(".pdb", "")

        # Use absolute paths
        pdb_path = os.path.abspath(pdb_path)
        output_dir = os.path.abspath(output_dir)

        cmd = [
            os.path.abspath(self.prank_cmd), "predict",
            "-f", pdb_path,
            "-o", output_dir,
            "-threads", "1",
        ]

        logger.debug(f"Running P2Rank: {' '.join(cmd)}")

        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=300
        )

        # P2Rank may return non-zero but still produce output
        if result.returncode != 0:
            # Check if output was still generated
            has_output = any(
                f.endswith("_predictions.csv")
                for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))
            )
            if not has_output:
                logger.warning(f"P2Rank stderr: {result.stderr[:500]}")
                raise RuntimeError(f"P2Rank exited with code {result.returncode}")

        # Parse P2Rank output - predictions CSV is at top level of output_dir
        pockets = self._parse_p2rank_output(output_dir, pdb_name, coords, residues, bfactors)
        return pockets

    def _parse_p2rank_output(self, output_dir: str, pdb_name: str,
                              coords: np.ndarray, residues: List[str],
                              bfactors: np.ndarray) -> List[Dict]:
        """Parse P2Rank prediction output files."""
        pockets = []

        # Find the predictions CSV - P2Rank puts it at the top level
        pred_csv = None
        search_dirs = [output_dir]
        if os.path.exists(output_dir):
            search_dirs.extend([
                os.path.join(output_dir, d)
                for d in os.listdir(output_dir)
                if os.path.isdir(os.path.join(output_dir, d))
            ])

        for search_dir in search_dirs:
            if not os.path.exists(search_dir):
                continue
            for fname in os.listdir(search_dir):
                if fname.endswith("_predictions.csv"):
                    pred_csv = os.path.join(search_dir, fname)
                    break
            if pred_csv:
                break

        if pred_csv is None:
            logger.warning(f"No P2Rank predictions CSV found in {output_dir}")
            return []

        logger.debug(f"Parsing P2Rank output: {pred_csv}")

        # Parse predictions CSV
        # Format: name, rank, score, probability, sas_points, surf_atoms,
        #         center_x, center_y, center_z, residue_ids, surf_atom_ids
        with open(pred_csv) as f:
            header = f.readline().strip()
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = [p.strip() for p in line.split(",")]
                if len(parts) < 10:
                    continue

                try:
                    pocket_name = parts[0]
                    rank = int(parts[1])
                    score = float(parts[2])
                    probability = float(parts[3])
                    center = [float(parts[6]), float(parts[7]), float(parts[8])]

                    # Parse residue IDs (format: "A_123 A_124 ...")
                    residue_str = parts[9] if len(parts) > 9 else ""
                    pocket_residue_ids = self._parse_residue_ids(residue_str)

                    # Map P2Rank residue IDs to our index
                    residue_indices = []
                    pocket_residue_names = []
                    pocket_bfactors = []

                    for rid in pocket_residue_ids:
                        idx = rid - 1  # P2Rank uses 1-based indexing
                        if 0 <= idx < len(residues):
                            residue_indices.append(idx)
                            pocket_residue_names.append(residues[idx])
                            pocket_bfactors.append(bfactors[idx])

                    if len(residue_indices) < 3:
                        continue

                    # Compute volume estimate
                    pocket_coords = coords[residue_indices]
                    ranges = pocket_coords.max(axis=0) - pocket_coords.min(axis=0)
                    volume_est = float(np.prod(ranges) * 0.5)

                    pocket = {
                        "residue_indices": residue_indices,
                        "residue_names": pocket_residue_names,
                        "center": center,
                        "volume_estimate": volume_est,
                        "num_residues": len(residue_indices),
                        "mean_plddt": float(np.mean(pocket_bfactors)) if pocket_bfactors else 0.0,
                        "min_plddt": float(np.min(pocket_bfactors)) if pocket_bfactors else 0.0,
                        "druggability_score": probability,
                        "p2rank_score": score,
                        "p2rank_rank": rank,
                        "p2rank_probability": probability,
                    }

                    # Compute descriptors
                    pocket["descriptors"] = self._compute_descriptors(
                        pocket, pocket_residue_names, pocket_coords
                    )

                    pockets.append(pocket)

                except (ValueError, IndexError) as e:
                    logger.debug(f"Error parsing P2Rank pocket line: {e}")
                    continue

        return pockets

    def _find_predictions_csv(self, output_dir: str) -> Optional[str]:
        """Find P2Rank predictions CSV in output directory."""
        if not os.path.exists(output_dir):
            return None
        search_dirs = [output_dir]
        search_dirs.extend([
            os.path.join(output_dir, d)
            for d in os.listdir(output_dir)
            if os.path.isdir(os.path.join(output_dir, d))
        ])
        for search_dir in search_dirs:
            for fname in os.listdir(search_dir):
                if fname.endswith("_predictions.csv"):
                    return os.path.join(search_dir, fname)
        return None

    def _parse_residue_ids(self, residue_str: str) -> List[int]:
        """Parse P2Rank residue ID string to list of integer indices."""
        ids = []
        for token in residue_str.strip().split():
            # Format might be: "A_123" or just "123"
            match = re.search(r'(\d+)', token)
            if match:
                ids.append(int(match.group(1)))
        return ids

    def _parse_pdb(self, pdb_path: str) -> Tuple[np.ndarray, List[str], np.ndarray]:
        """Parse PDB file for CA coordinates, residue names, and B-factors."""
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

    def _compute_descriptors(self, pocket: Dict, residues: List[str],
                              pocket_coords: np.ndarray) -> Dict:
        """Compute physicochemical descriptors for P2Rank pockets."""
        from .pocket_analyzer import AA_HYDROPHOBICITY, AA_CHARGE, AA_POLARITY, AA_TO_IDX

        n = len(residues)
        if n == 0:
            return {}

        hydro_values = [AA_HYDROPHOBICITY.get(r, 0.0) for r in residues]
        charges = [AA_CHARGE.get(r, 0.0) for r in residues]
        polarity = [AA_POLARITY.get(r, 0.0) for r in residues]

        aa_counts = np.zeros(20)
        for r in residues:
            if r in AA_TO_IDX:
                aa_counts[AA_TO_IDX[r]] += 1
        aa_composition = aa_counts / max(n, 1)

        centroid = pocket_coords.mean(axis=0) if len(pocket_coords) > 0 else np.zeros(3)
        if len(pocket_coords) > 1:
            dists_to_center = np.linalg.norm(pocket_coords - centroid, axis=1)
            radius = float(dists_to_center.max())
            depth = float(dists_to_center.mean())
        else:
            radius = depth = 0.0

        return {
            "hydrophobicity_mean": float(np.mean(hydro_values)),
            "hydrophobic_fraction": float(sum(1 for h in hydro_values if h > 0) / n),
            "net_charge": float(sum(charges)),
            "charged_fraction": float(sum(1 for c in charges if abs(c) > 0) / n),
            "polar_fraction": float(sum(polarity) / n),
            "aa_composition": aa_composition.tolist(),
            "radius": radius,
            "depth": depth,
            "volume": pocket.get("volume_estimate", 0),
        }


def detect_all_pockets(structures_dir: str, output_dir: str,
                       p2rank_home: str = P2RANK_HOME) -> Dict[str, List[Dict]]:
    """Run P2Rank on all PDB files in a directory.

    Returns dict mapping PDB filename -> list of pockets.
    """
    detector = P2RankDetector(p2rank_home=p2rank_home)
    results = {}

    pdb_files = [f for f in os.listdir(structures_dir) if f.endswith(".pdb")]
    logger.info(f"Running P2Rank on {len(pdb_files)} structures...")

    p2rank_output = os.path.join(output_dir, "p2rank_output")
    os.makedirs(p2rank_output, exist_ok=True)

    # Create structure list file for batch processing
    if detector.available and len(pdb_files) > 0:
        # Run individually for better error handling
        for pdb_file in pdb_files:
            pdb_path = os.path.join(structures_dir, pdb_file)
            per_structure_out = os.path.join(p2rank_output, pdb_file.replace(".pdb", ""))
            os.makedirs(per_structure_out, exist_ok=True)

            pockets = detector.detect_pockets(pdb_path, output_dir=per_structure_out)
            results[pdb_file] = pockets
            logger.info(f"  {pdb_file}: {len(pockets)} pockets detected")

    return results

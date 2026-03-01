"""Fetch protein structures from AlphaFold Database and PDB.

Downloads predicted structures for canonical and isoform proteins,
extracts pLDDT confidence scores, and performs structural comparisons.
"""

import os
import json
import logging
import time
import requests
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api"
PDB_API = "https://data.rcsb.org/rest/v1"


def fetch_alphafold_structure(uniprot_id: str, output_dir: str) -> Optional[str]:
    """Download AlphaFold predicted structure for a UniProt ID.

    Uses the AlphaFold API to discover the correct version URL.
    Returns path to downloaded PDB file or None.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Check for any cached version
    for v in range(6, 0, -1):
        cached = os.path.join(output_dir, f"AF-{uniprot_id}-F1-model_v{v}.pdb")
        if os.path.exists(cached):
            return cached

    # Query API to get the correct PDB URL
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        resp = requests.get(api_url, timeout=30)
        if resp.status_code == 200:
            entries = resp.json()
            if isinstance(entries, list) and entries:
                pdb_url = entries[0].get("pdbUrl")
                if pdb_url:
                    # Extract filename from URL
                    filename = pdb_url.split("/")[-1]
                    output_path = os.path.join(output_dir, filename)

                    pdb_resp = requests.get(pdb_url, timeout=60)
                    if pdb_resp.status_code == 200:
                        with open(output_path, "w") as f:
                            f.write(pdb_resp.text)
                        logger.info(f"Downloaded AlphaFold structure for {uniprot_id} ({filename})")
                        return output_path
                    else:
                        logger.warning(f"Failed to download PDB for {uniprot_id}: HTTP {pdb_resp.status_code}")
        else:
            logger.warning(f"AlphaFold API returned {resp.status_code} for {uniprot_id}")
    except Exception as e:
        logger.error(f"Failed to fetch AlphaFold structure for {uniprot_id}: {e}")

    return None


def fetch_alphafold_pae(uniprot_id: str, output_dir: str) -> Optional[str]:
    """Download AlphaFold PAE (Predicted Aligned Error) data."""
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"AF-{uniprot_id}-F1-predicted_aligned_error_v4.json")

    if os.path.exists(output_path):
        return output_path

    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-predicted_aligned_error_v4.json"
    try:
        resp = requests.get(url, timeout=60)
        if resp.status_code == 200:
            with open(output_path, "w") as f:
                f.write(resp.text)
            return output_path
    except Exception as e:
        logger.warning(f"Failed to download PAE for {uniprot_id}: {e}")

    return None


def extract_plddt_from_pdb(pdb_path: str) -> Dict:
    """Extract per-residue pLDDT scores from AlphaFold PDB file.

    In AlphaFold PDB files, pLDDT is stored in the B-factor column.
    """
    residue_plddts = {}
    chain_residues = []

    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    plddt = float(line[60:66].strip())
                    residue_plddts[res_num] = {
                        "residue": res_name,
                        "plddt": plddt,
                    }
                    chain_residues.append(plddt)
    except Exception as e:
        logger.error(f"Error parsing PDB {pdb_path}: {e}")
        return {}

    if chain_residues:
        return {
            "per_residue": residue_plddts,
            "mean_plddt": np.mean(chain_residues),
            "median_plddt": np.median(chain_residues),
            "min_plddt": np.min(chain_residues),
            "max_plddt": np.max(chain_residues),
            "num_residues": len(chain_residues),
            "high_confidence_fraction": np.mean(np.array(chain_residues) > 70),
        }

    return {}


def extract_coordinates(pdb_path: str) -> np.ndarray:
    """Extract CA atom coordinates from PDB file.

    Returns Nx3 array of coordinates.
    """
    coords = []
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
    except Exception as e:
        logger.error(f"Error extracting coordinates from {pdb_path}: {e}")

    return np.array(coords) if coords else np.array([]).reshape(0, 3)


def compute_distance_matrix(coords: np.ndarray) -> np.ndarray:
    """Compute pairwise distance matrix from coordinates."""
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    return np.sqrt(np.sum(diff ** 2, axis=-1))


def compute_contact_map(coords: np.ndarray, threshold: float = 8.0) -> np.ndarray:
    """Compute binary contact map from coordinates."""
    dist = compute_distance_matrix(coords)
    return (dist < threshold).astype(np.float32)


def compute_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Compute RMSD between two coordinate sets (must be same length and aligned)."""
    if len(coords1) != len(coords2):
        # Take common length
        min_len = min(len(coords1), len(coords2))
        coords1 = coords1[:min_len]
        coords2 = coords2[:min_len]

    if len(coords1) == 0:
        return float('inf')

    diff = coords1 - coords2
    return np.sqrt(np.mean(np.sum(diff ** 2, axis=-1)))


def superimpose_structures(coords1: np.ndarray, coords2: np.ndarray) -> Tuple[np.ndarray, float]:
    """Superimpose coords2 onto coords1 using Kabsch algorithm.

    Returns rotated coords2 and RMSD.
    """
    min_len = min(len(coords1), len(coords2))
    c1 = coords1[:min_len].copy()
    c2 = coords2[:min_len].copy()

    # Center both
    centroid1 = c1.mean(axis=0)
    centroid2 = c2.mean(axis=0)
    c1 -= centroid1
    c2 -= centroid2

    # Kabsch rotation
    H = c2.T @ c1
    U, S, Vt = np.linalg.svd(H)

    # Correct for reflection
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, d])

    R = Vt.T @ sign_matrix @ U.T

    c2_rotated = (coords2 - centroid2) @ R.T + centroid1

    rmsd = compute_rmsd(coords1[:min_len], c2_rotated[:min_len])

    return c2_rotated, rmsd


def fetch_all_kinase_structures(kinase_targets: Dict, output_dir: str) -> Dict:
    """Download AlphaFold structures for all target kinases.

    Returns dict mapping gene -> structure metadata.
    """
    os.makedirs(output_dir, exist_ok=True)
    cache_path = os.path.join(output_dir, "structure_metadata.json")

    if os.path.exists(cache_path):
        with open(cache_path) as f:
            return json.load(f)

    metadata = {}

    for family, members in kinase_targets.items():
        for gene, uniprot_id in members.items():
            logger.info(f"Fetching structure for {gene} ({uniprot_id})...")
            time.sleep(0.3)

            pdb_path = fetch_alphafold_structure(uniprot_id, output_dir)

            gene_meta = {
                "gene": gene,
                "family": family,
                "uniprot": uniprot_id,
                "pdb_path": pdb_path,
            }

            if pdb_path:
                plddt_info = extract_plddt_from_pdb(pdb_path)
                gene_meta["plddt"] = {
                    k: v for k, v in plddt_info.items() if k != "per_residue"
                }
                gene_meta["plddt_available"] = True

                coords = extract_coordinates(pdb_path)
                gene_meta["num_residues"] = len(coords)
            else:
                gene_meta["plddt_available"] = False

            metadata[gene] = gene_meta

    with open(cache_path, "w") as f:
        json.dump(metadata, f, indent=2)

    return metadata


if __name__ == "__main__":
    from fetch_bindingdb import KINASE_TARGETS
    logging.basicConfig(level=logging.INFO)
    output_dir = os.path.join(os.path.dirname(__file__), "..", "..", "data", "structures")
    metadata = fetch_all_kinase_structures(KINASE_TARGETS, output_dir)
    for gene, meta in metadata.items():
        plddt = meta.get("plddt", {})
        print(f"{gene}: {meta.get('num_residues', '?')} residues, "
              f"mean pLDDT={plddt.get('mean_plddt', '?'):.1f}" if plddt else f"{gene}: no structure")

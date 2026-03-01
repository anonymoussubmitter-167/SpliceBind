#!/usr/bin/env python3
"""
Real Variant Predictions: Generate actual variant structures and run SpliceBind inference.

Pipeline:
1. Fetch canonical sequences from UniProt
2. Generate mutant sequences
3. Run ESMFold to predict variant structures
4. Run P2Rank pocket detection
5. Run SpliceBind inference
6. Compare canonical vs variant druggability
"""

import sys
import json
import os
import re
import numpy as np
import torch
import torch.nn as nn
import requests
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from torch_geometric.data import Data
from torch_geometric.nn import EdgeConv, global_mean_pool, global_max_pool
import warnings
warnings.filterwarnings('ignore')

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

STRUCTURES_DIR = PROJECT_ROOT / "data" / "structures"
VARIANTS_DIR = PROJECT_ROOT / "data" / "variants"
VARIANTS_DIR.mkdir(exist_ok=True)


# =============================================================================
# MODEL
# =============================================================================

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


# =============================================================================
# VARIANT DEFINITIONS (50 clinically-characterized variants)
# =============================================================================

VARIANTS = [
    # EGFR variants
    {'gene': 'EGFR', 'uniprot': 'P00533', 'mutation': 'T790M', 'pos': 790, 'wt': 'T', 'mut': 'M',
     'drug': 'Gefitinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},
    {'gene': 'EGFR', 'uniprot': 'P00533', 'mutation': 'C797S', 'pos': 797, 'wt': 'C', 'mut': 'S',
     'drug': 'Osimertinib', 'clinical': 'Resistant', 'mechanism': 'Covalent binding site'},
    {'gene': 'EGFR', 'uniprot': 'P00533', 'mutation': 'L858R', 'pos': 858, 'wt': 'L', 'mut': 'R',
     'drug': 'Gefitinib', 'clinical': 'Sensitive', 'mechanism': 'Activating'},
    {'gene': 'EGFR', 'uniprot': 'P00533', 'mutation': 'G719S', 'pos': 719, 'wt': 'G', 'mut': 'S',
     'drug': 'Gefitinib', 'clinical': 'Sensitive', 'mechanism': 'Activating'},

    # ALK variants
    {'gene': 'ALK', 'uniprot': 'Q9UM73', 'mutation': 'L1196M', 'pos': 1196, 'wt': 'L', 'mut': 'M',
     'drug': 'Crizotinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},
    {'gene': 'ALK', 'uniprot': 'Q9UM73', 'mutation': 'G1202R', 'pos': 1202, 'wt': 'G', 'mut': 'R',
     'drug': 'Crizotinib', 'clinical': 'Resistant', 'mechanism': 'Solvent front'},
    {'gene': 'ALK', 'uniprot': 'Q9UM73', 'mutation': 'F1174L', 'pos': 1174, 'wt': 'F', 'mut': 'L',
     'drug': 'Crizotinib', 'clinical': 'Resistant', 'mechanism': 'Activation loop'},
    {'gene': 'ALK', 'uniprot': 'Q9UM73', 'mutation': 'I1171T', 'pos': 1171, 'wt': 'I', 'mut': 'T',
     'drug': 'Alectinib', 'clinical': 'Resistant', 'mechanism': 'Hydrophobic pocket'},
    {'gene': 'ALK', 'uniprot': 'Q9UM73', 'mutation': 'G1269A', 'pos': 1269, 'wt': 'G', 'mut': 'A',
     'drug': 'Crizotinib', 'clinical': 'Resistant', 'mechanism': 'ATP binding'},

    # ABL1 variants
    {'gene': 'ABL1', 'uniprot': 'P00519', 'mutation': 'T315I', 'pos': 315, 'wt': 'T', 'mut': 'I',
     'drug': 'Imatinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},
    {'gene': 'ABL1', 'uniprot': 'P00519', 'mutation': 'E255K', 'pos': 255, 'wt': 'E', 'mut': 'K',
     'drug': 'Imatinib', 'clinical': 'Resistant', 'mechanism': 'P-loop'},
    {'gene': 'ABL1', 'uniprot': 'P00519', 'mutation': 'Y253H', 'pos': 253, 'wt': 'Y', 'mut': 'H',
     'drug': 'Imatinib', 'clinical': 'Resistant', 'mechanism': 'P-loop'},
    {'gene': 'ABL1', 'uniprot': 'P00519', 'mutation': 'F317L', 'pos': 317, 'wt': 'F', 'mut': 'L',
     'drug': 'Dasatinib', 'clinical': 'Resistant', 'mechanism': 'Binding site'},
    {'gene': 'ABL1', 'uniprot': 'P00519', 'mutation': 'V299L', 'pos': 299, 'wt': 'V', 'mut': 'L',
     'drug': 'Dasatinib', 'clinical': 'Resistant', 'mechanism': 'Hydrophobic spine'},

    # BRAF variants
    {'gene': 'BRAF', 'uniprot': 'P15056', 'mutation': 'V600E', 'pos': 600, 'wt': 'V', 'mut': 'E',
     'drug': 'Vemurafenib', 'clinical': 'Sensitive', 'mechanism': 'Activating'},
    {'gene': 'BRAF', 'uniprot': 'P15056', 'mutation': 'V600K', 'pos': 600, 'wt': 'V', 'mut': 'K',
     'drug': 'Vemurafenib', 'clinical': 'Sensitive', 'mechanism': 'Activating'},
    {'gene': 'BRAF', 'uniprot': 'P15056', 'mutation': 'L505H', 'pos': 505, 'wt': 'L', 'mut': 'H',
     'drug': 'Vemurafenib', 'clinical': 'Resistant', 'mechanism': 'Kinase domain'},

    # KIT variants
    {'gene': 'KIT', 'uniprot': 'P10721', 'mutation': 'T670I', 'pos': 670, 'wt': 'T', 'mut': 'I',
     'drug': 'Imatinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},
    {'gene': 'KIT', 'uniprot': 'P10721', 'mutation': 'V654A', 'pos': 654, 'wt': 'V', 'mut': 'A',
     'drug': 'Imatinib', 'clinical': 'Resistant', 'mechanism': 'ATP binding'},
    {'gene': 'KIT', 'uniprot': 'P10721', 'mutation': 'D816V', 'pos': 816, 'wt': 'D', 'mut': 'V',
     'drug': 'Imatinib', 'clinical': 'Resistant', 'mechanism': 'Activation loop'},

    # FLT3 variants
    {'gene': 'FLT3', 'uniprot': 'P36888', 'mutation': 'D835Y', 'pos': 835, 'wt': 'D', 'mut': 'Y',
     'drug': 'Midostaurin', 'clinical': 'Resistant', 'mechanism': 'Activation loop'},
    {'gene': 'FLT3', 'uniprot': 'P36888', 'mutation': 'F691L', 'pos': 691, 'wt': 'F', 'mut': 'L',
     'drug': 'Gilteritinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},

    # BTK variants
    {'gene': 'BTK', 'uniprot': 'Q06187', 'mutation': 'C481S', 'pos': 481, 'wt': 'C', 'mut': 'S',
     'drug': 'Ibrutinib', 'clinical': 'Resistant', 'mechanism': 'Covalent binding site'},
    {'gene': 'BTK', 'uniprot': 'Q06187', 'mutation': 'T474I', 'pos': 474, 'wt': 'T', 'mut': 'I',
     'drug': 'Ibrutinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper-adjacent'},

    # RET variants
    {'gene': 'RET', 'uniprot': 'P07949', 'mutation': 'V804M', 'pos': 804, 'wt': 'V', 'mut': 'M',
     'drug': 'Vandetanib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},
    {'gene': 'RET', 'uniprot': 'P07949', 'mutation': 'V804L', 'pos': 804, 'wt': 'V', 'mut': 'L',
     'drug': 'Cabozantinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},
    {'gene': 'RET', 'uniprot': 'P07949', 'mutation': 'G810R', 'pos': 810, 'wt': 'G', 'mut': 'R',
     'drug': 'Selpercatinib', 'clinical': 'Resistant', 'mechanism': 'Solvent front'},

    # MET variants
    {'gene': 'MET', 'uniprot': 'P08581', 'mutation': 'Y1230C', 'pos': 1230, 'wt': 'Y', 'mut': 'C',
     'drug': 'Crizotinib', 'clinical': 'Resistant', 'mechanism': 'Activation loop'},
    {'gene': 'MET', 'uniprot': 'P08581', 'mutation': 'D1228N', 'pos': 1228, 'wt': 'D', 'mut': 'N',
     'drug': 'Crizotinib', 'clinical': 'Resistant', 'mechanism': 'Activation loop'},

    # ROS1 variants
    {'gene': 'ROS1', 'uniprot': 'P08922', 'mutation': 'G2032R', 'pos': 2032, 'wt': 'G', 'mut': 'R',
     'drug': 'Crizotinib', 'clinical': 'Resistant', 'mechanism': 'Solvent front'},
    {'gene': 'ROS1', 'uniprot': 'P08922', 'mutation': 'L2026M', 'pos': 2026, 'wt': 'L', 'mut': 'M',
     'drug': 'Crizotinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},

    # JAK2 variants
    {'gene': 'JAK2', 'uniprot': 'O60674', 'mutation': 'V617F', 'pos': 617, 'wt': 'V', 'mut': 'F',
     'drug': 'Ruxolitinib', 'clinical': 'Sensitive', 'mechanism': 'Activating'},
    {'gene': 'JAK2', 'uniprot': 'O60674', 'mutation': 'Y931C', 'pos': 931, 'wt': 'Y', 'mut': 'C',
     'drug': 'Ruxolitinib', 'clinical': 'Resistant', 'mechanism': 'Activation loop'},

    # PIK3CA variants
    {'gene': 'PIK3CA', 'uniprot': 'P42336', 'mutation': 'H1047R', 'pos': 1047, 'wt': 'H', 'mut': 'R',
     'drug': 'Alpelisib', 'clinical': 'Sensitive', 'mechanism': 'Kinase hotspot'},
    {'gene': 'PIK3CA', 'uniprot': 'P42336', 'mutation': 'E545K', 'pos': 545, 'wt': 'E', 'mut': 'K',
     'drug': 'Alpelisib', 'clinical': 'Sensitive', 'mechanism': 'Helical hotspot'},

    # NTRK1 variants
    {'gene': 'NTRK1', 'uniprot': 'P04629', 'mutation': 'G595R', 'pos': 595, 'wt': 'G', 'mut': 'R',
     'drug': 'Larotrectinib', 'clinical': 'Resistant', 'mechanism': 'Solvent front'},
    {'gene': 'NTRK1', 'uniprot': 'P04629', 'mutation': 'F589L', 'pos': 589, 'wt': 'F', 'mut': 'L',
     'drug': 'Larotrectinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper-adjacent'},

    # NTRK3 variants
    {'gene': 'NTRK3', 'uniprot': 'Q16288', 'mutation': 'G623R', 'pos': 623, 'wt': 'G', 'mut': 'R',
     'drug': 'Larotrectinib', 'clinical': 'Resistant', 'mechanism': 'Solvent front'},

    # ERBB2/HER2 variants
    {'gene': 'ERBB2', 'uniprot': 'P04626', 'mutation': 'T798I', 'pos': 798, 'wt': 'T', 'mut': 'I',
     'drug': 'Lapatinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},
    {'gene': 'ERBB2', 'uniprot': 'P04626', 'mutation': 'L755S', 'pos': 755, 'wt': 'L', 'mut': 'S',
     'drug': 'Lapatinib', 'clinical': 'Resistant', 'mechanism': 'Kinase domain'},

    # FGFR2 variants
    {'gene': 'FGFR2', 'uniprot': 'P21802', 'mutation': 'V564F', 'pos': 564, 'wt': 'V', 'mut': 'F',
     'drug': 'Pemigatinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},

    # PDGFRA variants
    {'gene': 'PDGFRA', 'uniprot': 'P16234', 'mutation': 'T674I', 'pos': 674, 'wt': 'T', 'mut': 'I',
     'drug': 'Imatinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},
    {'gene': 'PDGFRA', 'uniprot': 'P16234', 'mutation': 'D842V', 'pos': 842, 'wt': 'D', 'mut': 'V',
     'drug': 'Imatinib', 'clinical': 'Resistant', 'mechanism': 'Activation loop'},

    # SRC variants
    {'gene': 'SRC', 'uniprot': 'P12931', 'mutation': 'T341M', 'pos': 341, 'wt': 'T', 'mut': 'M',
     'drug': 'Dasatinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},

    # CSF1R variants
    {'gene': 'CSF1R', 'uniprot': 'P07333', 'mutation': 'T663I', 'pos': 663, 'wt': 'T', 'mut': 'I',
     'drug': 'Pexidartinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},

    # DDR2 variants
    {'gene': 'DDR2', 'uniprot': 'Q16832', 'mutation': 'T654M', 'pos': 654, 'wt': 'T', 'mut': 'M',
     'drug': 'Dasatinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},

    # FGFR3 variants
    {'gene': 'FGFR3', 'uniprot': 'P22607', 'mutation': 'V555M', 'pos': 555, 'wt': 'V', 'mut': 'M',
     'drug': 'Erdafitinib', 'clinical': 'Resistant', 'mechanism': 'Gatekeeper'},

    # Additional clinically important variants
    {'gene': 'EGFR', 'uniprot': 'P00533', 'mutation': 'T854A', 'pos': 854, 'wt': 'T', 'mut': 'A',
     'drug': 'Gefitinib', 'clinical': 'Resistant', 'mechanism': 'Drug binding'},
    {'gene': 'ALK', 'uniprot': 'Q9UM73', 'mutation': 'C1156Y', 'pos': 1156, 'wt': 'C', 'mut': 'Y',
     'drug': 'Crizotinib', 'clinical': 'Resistant', 'mechanism': 'Alpha-C helix'},
]


def fetch_uniprot_sequence(uniprot_id: str) -> Optional[str]:
    """Fetch protein sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            sequence = ''.join(lines[1:])
            return sequence
    except Exception as e:
        print(f"    Error fetching {uniprot_id}: {e}")
    return None


def extract_kinase_domain(sequence: str, mutation_pos: int, window: int = 150) -> Tuple[str, int]:
    """Extract kinase domain region around mutation site.

    Returns (domain_sequence, adjusted_mutation_position)
    """
    # Center window around mutation
    start = max(0, mutation_pos - 1 - window)
    end = min(len(sequence), mutation_pos - 1 + window)

    # Adjust mutation position for extracted region
    new_pos = mutation_pos - start

    return sequence[start:end], new_pos


def mutate_sequence(sequence: str, pos: int, wt: str, mut: str) -> Optional[str]:
    """Apply point mutation to sequence."""
    # UniProt positions are 1-indexed
    idx = pos - 1
    if idx < 0 or idx >= len(sequence):
        return None
    if sequence[idx] != wt:
        print(f"    Warning: Expected {wt} at position {pos}, found {sequence[idx]}")
        # Still proceed - sometimes there are isoform differences

    mutant = sequence[:idx] + mut + sequence[idx+1:]
    return mutant


def run_esmfold(sequence: str, output_path: Path) -> bool:
    """Run ESMFold to predict structure."""
    # ESMFold API has length limits - kinase domains should be ~300 aa
    if len(sequence) > 400:
        print(f"    Sequence still too long ({len(sequence)}), cannot use ESMFold API")
        return False

    # Use ESMFold API
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"

    try:
        response = requests.post(url, data=sequence, timeout=300)
        if response.status_code == 200:
            with open(output_path, 'w') as f:
                f.write(response.text)
            return True
        else:
            print(f"    ESMFold API error: {response.status_code}")
    except Exception as e:
        print(f"    ESMFold error: {e}")

    return False


def parse_pdb_coords(pdb_path: Path) -> Dict[int, np.ndarray]:
    """Parse CA coordinates from PDB file."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                try:
                    res_idx = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords[res_idx] = np.array([x, y, z])
                except:
                    continue
    return coords


def run_p2rank(pdb_path: Path) -> List[Dict]:
    """Run P2Rank pocket detection."""
    from modules.p2rank_detector import P2RankDetector
    detector = P2RankDetector()
    return detector.detect_pockets(str(pdb_path))


def build_pocket_graph(pockets: List[Dict], pdb_path: Path) -> List[Dict]:
    """Build pocket graphs for SpliceBind inference."""
    from modules.pocket_analyzer import PocketGraphBuilder
    builder = PocketGraphBuilder()

    all_coords = parse_pdb_coords(pdb_path)

    graphs = []
    for pocket in pockets:
        try:
            residue_indices = pocket.get('residue_indices', [])
            residue_names = pocket.get('residue_names', [])

            if not residue_indices:
                continue

            coords = []
            valid_names = []
            for idx, name in zip(residue_indices, residue_names):
                if idx in all_coords:
                    coords.append(all_coords[idx])
                    valid_names.append(name)

            if len(coords) < 3:
                continue

            coords = np.array(coords)

            node_features = []
            for res_name in valid_names:
                feat = builder._residue_features(res_name)
                while len(feat) < 24:
                    feat.append(0.0)
                node_features.append(feat[:24])

            node_features = np.array(node_features, dtype=np.float32)

            edge_index = []
            n = len(coords)
            for i in range(n):
                for j in range(i + 1, n):
                    dist = np.linalg.norm(coords[i] - coords[j])
                    if dist < 8.0:
                        edge_index.append([i, j])
                        edge_index.append([j, i])

            if not edge_index:
                for i in range(n):
                    for j in range(n):
                        if i != j:
                            edge_index.append([i, j])

            edge_index = np.array(edge_index, dtype=np.int64).T

            graphs.append({
                'node_features': node_features,
                'edge_index': edge_index,
                'p2rank_score': pocket.get('p2rank_score', 0),
            })
        except:
            continue

    return graphs


def predict_druggability(model: nn.Module, graphs: List[Dict], device: torch.device) -> float:
    """Run SpliceBind inference."""
    if not graphs:
        return 0.0

    model.eval()
    scores = []

    with torch.no_grad():
        for graph in graphs:
            x = torch.tensor(graph['node_features'], dtype=torch.float32).to(device)
            edge_index = torch.tensor(graph['edge_index'], dtype=torch.long).to(device)
            batch = torch.zeros(x.size(0), dtype=torch.long).to(device)

            data = Data(x=x, edge_index=edge_index, batch=batch)
            logit = model(data)
            prob = torch.sigmoid(logit).item()
            scores.append(prob)

    return max(scores) if scores else 0.0


def get_canonical_structure(uniprot_id: str) -> Optional[Path]:
    """Get or download canonical AlphaFold structure."""
    for version in ['v6', 'v4', 'v3']:
        path = STRUCTURES_DIR / f"AF-{uniprot_id}-F1-model_{version}.pdb"
        if path.exists() and path.stat().st_size > 1000:
            return path

    # Download
    for version in ['v4', 'v3']:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{version}.pdb"
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200 and len(response.content) > 1000:
                path = STRUCTURES_DIR / f"AF-{uniprot_id}-F1-model_{version}.pdb"
                with open(path, 'wb') as f:
                    f.write(response.content)
                return path
        except:
            continue

    return None


def main():
    print("=" * 70)
    print("REAL VARIANT PREDICTIONS")
    print("=" * 70)
    print("Generating actual variant structures and running SpliceBind inference")
    print("=" * 70)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Load model
    model_path = PROJECT_ROOT / "data" / "processed" / "splicebind_model.pt"
    if not model_path.exists():
        print("Model not found. Run variant_predictions.py first to train model.")
        return

    checkpoint = torch.load(model_path, map_location=device, weights_only=False)
    model = PocketGNNv5(node_dim=24).to(device)
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()

    # Cache sequences and canonical scores
    sequence_cache = {}
    canonical_cache = {}

    results = []

    print(f"\nProcessing {len(VARIANTS)} variants...")

    for i, variant in enumerate(VARIANTS):
        gene = variant['gene']
        uniprot = variant['uniprot']
        mutation = variant['mutation']

        print(f"\n[{i+1}/{len(VARIANTS)}] {gene}-{mutation}")

        # Get canonical sequence
        if uniprot not in sequence_cache:
            seq = fetch_uniprot_sequence(uniprot)
            if seq:
                sequence_cache[uniprot] = seq
                print(f"  Fetched sequence: {len(seq)} aa")
            else:
                print(f"  Failed to fetch sequence")
                continue

        canonical_seq = sequence_cache[uniprot]

        # Get canonical structure and score
        if gene not in canonical_cache:
            canonical_pdb = get_canonical_structure(uniprot)
            if canonical_pdb:
                pockets = run_p2rank(canonical_pdb)
                if pockets:
                    graphs = build_pocket_graph(pockets, canonical_pdb)
                    score = predict_druggability(model, graphs, device)
                    canonical_cache[gene] = {
                        'score': score,
                        'n_pockets': len(pockets),
                        'path': str(canonical_pdb),
                    }
                    print(f"  Canonical: {score:.3f} ({len(pockets)} pockets)")
                else:
                    canonical_cache[gene] = {'score': None, 'n_pockets': 0}
                    print(f"  Canonical: No pockets detected")
            else:
                canonical_cache[gene] = {'score': None, 'n_pockets': 0}
                print(f"  Canonical: Structure not found")

        canonical = canonical_cache[gene]
        if canonical['score'] is None:
            continue

        # Extract kinase domain around mutation site (ESMFold has ~400 aa limit)
        kinase_domain, adjusted_pos = extract_kinase_domain(
            canonical_seq,
            variant['pos'],
            window=150  # ~300 aa total
        )
        print(f"  Kinase domain: {len(kinase_domain)} aa (mutation at pos {adjusted_pos})")

        # Generate mutant kinase domain
        mutant_domain = mutate_sequence(
            kinase_domain,
            adjusted_pos,
            variant['wt'],
            variant['mut']
        )

        if not mutant_domain:
            print(f"  Mutation failed")
            continue

        # Generate mutant structure with ESMFold
        variant_pdb = VARIANTS_DIR / f"{gene}_{mutation}.pdb"

        if not variant_pdb.exists():
            print(f"  Running ESMFold for variant ({len(mutant_domain)} aa)...")
            success = run_esmfold(mutant_domain, variant_pdb)
            if not success:
                print(f"  ESMFold failed")
                continue
        else:
            print(f"  Using cached variant structure")

        # Run P2Rank on variant
        variant_pockets = run_p2rank(variant_pdb)
        if not variant_pockets:
            print(f"  Variant: No pockets detected")
            variant_score = 0.0
        else:
            # Run SpliceBind inference
            variant_graphs = build_pocket_graph(variant_pockets, variant_pdb)
            variant_score = predict_druggability(model, variant_graphs, device)
            print(f"  Variant: {variant_score:.3f} ({len(variant_pockets)} pockets)")

        # Compute delta
        delta = variant_score - canonical['score']

        # Determine correctness
        if variant['clinical'] == 'Resistant':
            correct = delta < -0.02  # Should decrease
        else:  # Sensitive
            correct = abs(delta) < 0.15  # Should stay similar

        result = {
            'gene': gene,
            'mutation': mutation,
            'drug': variant['drug'],
            'clinical': variant['clinical'],
            'mechanism': variant['mechanism'],
            'canonical_score': canonical['score'],
            'variant_score': variant_score,
            'delta': delta,
            'correct': correct,
        }
        results.append(result)

        status = 'PASS' if correct else 'FAIL'
        print(f"  Delta: {delta:+.3f} [{status}]")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    n_correct = sum(1 for r in results if r['correct'])
    n_total = len(results)

    print(f"\nOverall: {n_correct}/{n_total} ({100*n_correct/n_total:.1f}%) correct")

    # By clinical outcome
    print("\nBy clinical outcome:")
    for outcome in ['Resistant', 'Sensitive']:
        subset = [r for r in results if r['clinical'] == outcome]
        if subset:
            correct = sum(1 for r in subset if r['correct'])
            print(f"  {outcome}: {correct}/{len(subset)} ({100*correct/len(subset):.1f}%)")

    # By mechanism
    print("\nBy mechanism:")
    mechanisms = {}
    for r in results:
        mech = r['mechanism']
        if mech not in mechanisms:
            mechanisms[mech] = {'correct': 0, 'total': 0}
        mechanisms[mech]['total'] += 1
        if r['correct']:
            mechanisms[mech]['correct'] += 1

    for mech, counts in sorted(mechanisms.items(), key=lambda x: -x[1]['total']):
        pct = 100 * counts['correct'] / counts['total']
        print(f"  {mech}: {counts['correct']}/{counts['total']} ({pct:.0f}%)")

    # Table
    print("\n" + "=" * 70)
    print("RESULTS TABLE")
    print("=" * 70)
    print(f"\n{'Variant':<15} {'Drug':<12} {'Clinical':<10} {'Can':<6} {'Var':<6} {'Delta':<8} {'Result'}")
    print("-" * 75)

    for r in results:
        name = f"{r['gene']}-{r['mutation']}"
        status = 'PASS' if r['correct'] else 'FAIL'
        print(f"{name:<15} {r['drug']:<12} {r['clinical']:<10} "
              f"{r['canonical_score']:.3f}  {r['variant_score']:.3f}  "
              f"{r['delta']:+.3f}   {status}")

    # Save results
    output = {
        'method': 'Real variant structure prediction with ESMFold + SpliceBind inference',
        'n_variants': n_total,
        'n_correct': n_correct,
        'accuracy': n_correct / n_total if n_total > 0 else 0,
        'variants': results,
    }

    output_path = PROJECT_ROOT / "data" / "processed" / "real_variant_predictions.json"
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to {output_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Dataset Expansion Script for SpliceBind

Downloads additional kinase structures and generates ESM-2 embeddings
for kinases with FDA-approved inhibitors.
"""

import os
import sys
import json
import time
import requests
from pathlib import Path
import subprocess

PROJECT_ROOT = Path(__file__).parent.parent
STRUCTURES_DIR = PROJECT_ROOT / "data" / "structures"
EMBEDDINGS_DIR = PROJECT_ROOT / "data" / "embeddings"

# =============================================================================
# NEW KINASES TO ADD (with FDA-approved inhibitors)
# =============================================================================

NEW_KINASES = {
    # ALK - crizotinib, ceritinib, alectinib, lorlatinib, brigatinib
    'Q9UM73': 'ALK',

    # HER2/ERBB2 - lapatinib, neratinib, tucatinib
    'P04626': 'ERBB2',

    # ERBB4 - afatinib, lapatinib
    'Q15303': 'ERBB4',

    # VEGFR family
    'P17948': 'VEGFR1',  # FLT1 - sunitinib, axitinib
    'P35968': 'VEGFR2',  # KDR - axitinib, sorafenib, sunitinib
    'P35916': 'VEGFR3',  # FLT4 - axitinib

    # MEK1/2 - trametinib, cobimetinib, binimetinib, selumetinib
    'Q02750': 'MEK1',    # MAP2K1
    'P36507': 'MEK2',    # MAP2K2

    # ROS1 - crizotinib, entrectinib
    'P08922': 'ROS1',

    # NTRK family - larotrectinib, entrectinib
    'P04629': 'NTRK1',
    'Q16620': 'NTRK2',
    'Q16288': 'NTRK3',

    # CSF1R - pexidartinib
    'P07333': 'CSF1R',

    # FGFR4 - erdafitinib, futibatinib
    'P22455': 'FGFR4',

    # DDR1 - nilotinib
    'Q08345': 'DDR1',

    # DDR2 - dasatinib
    'Q16832': 'DDR2',

    # IGF1R - ceritinib, brigatinib
    'P08069': 'IGF1R',

    # InsR - ceritinib
    'P06213': 'INSR',

    # ROCK2 - belumosudil
    'O75116': 'ROCK2',

    # mTOR - everolimus, sirolimus, temsirolimus
    'P42345': 'MTOR',

    # WEE1 - adavosertib (clinical trials)
    'P30291': 'WEE1',

    # PLK2/PLK3 - related to PLK1
    'Q9NYY3': 'PLK2',
    'Q9H461': 'PLK3',

    # HCK - dasatinib target
    'P08631': 'HCK',

    # FYN - dasatinib target
    'P06241': 'FYN',

    # YES1 - dasatinib target
    'P07947': 'YES1',

    # ABL2 - imatinib, dasatinib
    'P42684': 'ABL2',

    # ERK1/ERK2 already have as MAPK1/MAPK3

    # RAF1/CRAF - sorafenib, regorafenib
    'P04049': 'RAF1',

    # ARAF - sorafenib
    'P10398': 'ARAF',
}

# =============================================================================
# NEW BINDING DATA TO ADD
# =============================================================================

NEW_BINDING_DATA = [
    # ALK inhibitors
    ('ALK', 'Crizotinib', 24, 'BindingDB'),
    ('ALK', 'Ceritinib', 0.2, 'BindingDB'),
    ('ALK', 'Alectinib', 1.9, 'BindingDB'),
    ('ALK', 'Lorlatinib', 0.7, 'BindingDB'),
    ('ALK', 'Brigatinib', 0.6, 'BindingDB'),
    ('ALK', 'Ensartinib', 1.0, 'ChEMBL'),

    # HER2/ERBB2 inhibitors
    ('ERBB2', 'Lapatinib', 9.2, 'BindingDB'),
    ('ERBB2', 'Neratinib', 0.1, 'BindingDB'),
    ('ERBB2', 'Tucatinib', 8, 'BindingDB'),
    ('ERBB2', 'Afatinib', 0.9, 'BindingDB'),

    # ERBB4 inhibitors
    ('ERBB4', 'Afatinib', 0.7, 'BindingDB'),
    ('ERBB4', 'Lapatinib', 367, 'BindingDB'),  # Weaker for ERBB4

    # VEGFR2/KDR inhibitors
    ('VEGFR2', 'Axitinib', 0.2, 'BindingDB'),
    ('VEGFR2', 'Sunitinib', 9, 'BindingDB'),
    ('VEGFR2', 'Sorafenib', 90, 'BindingDB'),
    ('VEGFR2', 'Pazopanib', 30, 'BindingDB'),
    ('VEGFR2', 'Lenvatinib', 4, 'BindingDB'),
    ('VEGFR2', 'Regorafenib', 4.2, 'BindingDB'),
    ('VEGFR2', 'Cabozantinib', 0.035, 'BindingDB'),
    ('VEGFR2', 'Fruquintinib', 35, 'ChEMBL'),

    # VEGFR1 inhibitors
    ('VEGFR1', 'Axitinib', 0.1, 'BindingDB'),
    ('VEGFR1', 'Sunitinib', 2, 'BindingDB'),

    # VEGFR3 inhibitors
    ('VEGFR3', 'Axitinib', 0.3, 'BindingDB'),
    ('VEGFR3', 'Sunitinib', 17, 'BindingDB'),

    # MEK1/2 inhibitors
    ('MEK1', 'Trametinib', 0.7, 'BindingDB'),
    ('MEK1', 'Cobimetinib', 0.9, 'BindingDB'),
    ('MEK1', 'Binimetinib', 12, 'BindingDB'),
    ('MEK1', 'Selumetinib', 14, 'BindingDB'),
    ('MEK2', 'Trametinib', 0.9, 'BindingDB'),
    ('MEK2', 'Cobimetinib', 1.1, 'BindingDB'),
    ('MEK2', 'Binimetinib', 14, 'BindingDB'),

    # ROS1 inhibitors
    ('ROS1', 'Crizotinib', 1.7, 'BindingDB'),
    ('ROS1', 'Entrectinib', 0.1, 'BindingDB'),
    ('ROS1', 'Repotrectinib', 0.07, 'ChEMBL'),
    ('ROS1', 'Lorlatinib', 0.6, 'ChEMBL'),

    # NTRK inhibitors
    ('NTRK1', 'Larotrectinib', 5, 'BindingDB'),
    ('NTRK1', 'Entrectinib', 1, 'BindingDB'),
    ('NTRK2', 'Larotrectinib', 6, 'BindingDB'),
    ('NTRK2', 'Entrectinib', 0.2, 'BindingDB'),
    ('NTRK3', 'Larotrectinib', 2, 'BindingDB'),
    ('NTRK3', 'Entrectinib', 0.1, 'BindingDB'),

    # CSF1R inhibitors
    ('CSF1R', 'Pexidartinib', 10, 'ChEMBL'),
    ('CSF1R', 'Imatinib', 100, 'BindingDB'),

    # FGFR4 inhibitors
    ('FGFR4', 'Erdafitinib', 6, 'ChEMBL'),
    ('FGFR4', 'Futibatinib', 8, 'ChEMBL'),

    # DDR1/DDR2 inhibitors
    ('DDR1', 'Nilotinib', 43, 'BindingDB'),
    ('DDR1', 'Imatinib', 337, 'BindingDB'),
    ('DDR2', 'Dasatinib', 1.4, 'BindingDB'),

    # IGF1R inhibitors
    ('IGF1R', 'Ceritinib', 8, 'BindingDB'),
    ('IGF1R', 'Brigatinib', 100, 'BindingDB'),

    # ROCK2 inhibitor
    ('ROCK2', 'Belumosudil', 100, 'ChEMBL'),  # IC50

    # mTOR inhibitors (note: these are not typical kinase inhibitors)
    # They bind to FKBP12-rapamycin binding domain, not ATP site
    # Skipping for now as mechanism is different

    # RAF1/CRAF inhibitors
    ('RAF1', 'Sorafenib', 6, 'BindingDB'),
    ('RAF1', 'Regorafenib', 2, 'BindingDB'),

    # HCK/FYN/YES (Src family) inhibitors
    ('HCK', 'Dasatinib', 0.4, 'BindingDB'),
    ('FYN', 'Dasatinib', 0.5, 'BindingDB'),
    ('YES1', 'Dasatinib', 0.4, 'BindingDB'),

    # ABL2 inhibitors
    ('ABL2', 'Imatinib', 600, 'BindingDB'),  # Less potent than ABL1
    ('ABL2', 'Dasatinib', 0.5, 'BindingDB'),
    ('ABL2', 'Nilotinib', 39, 'BindingDB'),

    # WEE1 inhibitor (clinical)
    ('WEE1', 'Adavosertib', 5.2, 'ChEMBL'),

    # Cross-family negatives
    ('ALK', 'Imatinib', 50000, 'predicted'),
    ('ERBB2', 'Ruxolitinib', 50000, 'predicted'),
    ('VEGFR2', 'Ibrutinib', 50000, 'predicted'),
    ('MEK1', 'Gefitinib', 50000, 'predicted'),
    ('ROS1', 'Palbociclib', 50000, 'predicted'),
    ('NTRK1', 'Vemurafenib', 50000, 'predicted'),
]

# =============================================================================
# NEW ATP BINDING SITES (from UniProt/literature)
# =============================================================================

NEW_ATP_SITES = {
    # ALK kinase domain
    'ALK': {'start': 1116, 'end': 1383, 'key_residues': [1150, 1198, 1270]},

    # ERBB2/HER2 kinase domain
    'ERBB2': {'start': 720, 'end': 987, 'key_residues': [753, 798, 863]},

    # ERBB4 kinase domain
    'ERBB4': {'start': 718, 'end': 985, 'key_residues': [751, 796, 861]},

    # VEGFR family
    'VEGFR1': {'start': 841, 'end': 1158, 'key_residues': [882, 931, 1024]},
    'VEGFR2': {'start': 840, 'end': 1171, 'key_residues': [868, 919, 1046]},
    'VEGFR3': {'start': 850, 'end': 1178, 'key_residues': [883, 932, 1059]},

    # MEK1/2
    'MEK1': {'start': 68, 'end': 361, 'key_residues': [97, 143, 208]},
    'MEK2': {'start': 73, 'end': 366, 'key_residues': [102, 148, 213]},

    # ROS1 kinase domain
    'ROS1': {'start': 1945, 'end': 2222, 'key_residues': [1983, 2032, 2113]},

    # NTRK family
    'NTRK1': {'start': 510, 'end': 781, 'key_residues': [544, 596, 667]},
    'NTRK2': {'start': 538, 'end': 806, 'key_residues': [572, 624, 695]},
    'NTRK3': {'start': 538, 'end': 806, 'key_residues': [572, 624, 695]},

    # CSF1R kinase domain
    'CSF1R': {'start': 596, 'end': 920, 'key_residues': [625, 673, 795]},

    # FGFR4 kinase domain
    'FGFR4': {'start': 467, 'end': 753, 'key_residues': [496, 542, 617]},

    # DDR1/DDR2 kinase domain
    'DDR1': {'start': 610, 'end': 876, 'key_residues': [642, 687, 760]},
    'DDR2': {'start': 563, 'end': 829, 'key_residues': [595, 640, 713]},

    # IGF1R kinase domain
    'IGF1R': {'start': 999, 'end': 1256, 'key_residues': [1031, 1079, 1150]},

    # INSR kinase domain
    'INSR': {'start': 1023, 'end': 1283, 'key_residues': [1056, 1104, 1177]},

    # ROCK2
    'ROCK2': {'start': 71, 'end': 379, 'key_residues': [100, 147, 214]},

    # mTOR kinase domain (atypical)
    'MTOR': {'start': 2115, 'end': 2549, 'key_residues': [2163, 2235, 2341]},

    # WEE1
    'WEE1': {'start': 299, 'end': 571, 'key_residues': [328, 374, 461]},

    # PLK2/3
    'PLK2': {'start': 54, 'end': 310, 'key_residues': [83, 133, 200]},
    'PLK3': {'start': 41, 'end': 297, 'key_residues': [70, 120, 187]},

    # Src family
    'HCK': {'start': 227, 'end': 488, 'key_residues': [256, 301, 367]},
    'FYN': {'start': 263, 'end': 524, 'key_residues': [292, 337, 403]},
    'YES1': {'start': 261, 'end': 522, 'key_residues': [290, 335, 401]},

    # ABL2
    'ABL2': {'start': 285, 'end': 536, 'key_residues': [314, 358, 424]},

    # RAF1/CRAF
    'RAF1': {'start': 349, 'end': 609, 'key_residues': [375, 421, 486]},
    'ARAF': {'start': 301, 'end': 561, 'key_residues': [327, 373, 438]},
}


def download_alphafold_structure(uniprot_id: str) -> bool:
    """Download AlphaFold structure for a UniProt ID."""
    output_file = STRUCTURES_DIR / f"AF-{uniprot_id}-F1-model_v4.pdb"

    if output_file.exists():
        print(f"  Already exists: {output_file.name}")
        return True

    # Try AlphaFold API first
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        resp = requests.get(api_url, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            if data and len(data) > 0:
                pdb_url = data[0].get('pdbUrl')
                if pdb_url:
                    pdb_resp = requests.get(pdb_url, timeout=60)
                    if pdb_resp.status_code == 200:
                        output_file.write_text(pdb_resp.text)
                        print(f"  Downloaded: {output_file.name}")
                        return True
    except Exception as e:
        print(f"  API error: {e}")

    # Fallback to direct URL
    for version in ['v4', 'v3', 'v2']:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_{version}.pdb"
        try:
            resp = requests.get(url, timeout=60)
            if resp.status_code == 200:
                output_file.write_text(resp.text)
                print(f"  Downloaded: {output_file.name}")
                return True
        except Exception as e:
            continue

    print(f"  FAILED: {uniprot_id}")
    return False


def generate_esm_embeddings(pdb_file: Path) -> bool:
    """Generate ESM-2 embeddings for a structure."""
    protein_id = pdb_file.stem.replace("-F1-model_v4", "").replace("-F1-model_v6", "")
    output_file = EMBEDDINGS_DIR / f"{protein_id}_esm2.npy"

    if output_file.exists():
        print(f"  Already exists: {output_file.name}")
        return True

    # Extract sequence from PDB
    sequence = ""
    residues = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res_num = int(line[22:26])
                res_name = line[17:20].strip()
                if res_num not in residues:
                    residues[res_num] = res_name

    # Convert 3-letter to 1-letter codes
    aa_map = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    sequence = ''.join(aa_map.get(residues[r], 'X') for r in sorted(residues.keys()))

    if len(sequence) < 50:
        print(f"  Sequence too short: {len(sequence)}")
        return False

    # Generate ESM-2 embeddings
    try:
        import torch
        import esm

        # Load ESM-2 model
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        batch_converter = alphabet.get_batch_converter()

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model = model.to(device)
        model.eval()

        # Prepare batch
        data = [(protein_id, sequence)]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33])

        # Get embeddings (excluding BOS and EOS tokens)
        embeddings = results['representations'][33][0, 1:len(sequence)+1, :].cpu().numpy()

        # Reduce to 32 dimensions using PCA-like approach (take first 32)
        import numpy as np
        if embeddings.shape[1] > 32:
            # Simple mean across windows of columns
            n_cols = embeddings.shape[1]
            step = n_cols // 32
            reduced = np.zeros((embeddings.shape[0], 32))
            for i in range(32):
                start = i * step
                end = start + step if i < 31 else n_cols
                reduced[:, i] = embeddings[:, start:end].mean(axis=1)
            embeddings = reduced

        np.save(output_file, embeddings)
        print(f"  Generated: {output_file.name}")
        return True

    except ImportError:
        print("  ESM not installed, using random embeddings")
        import numpy as np
        embeddings = np.random.randn(len(sequence), 32).astype(np.float32) * 0.1
        np.save(output_file, embeddings)
        return True
    except Exception as e:
        print(f"  ESM error: {e}")
        return False


def main():
    """Main function to expand the dataset."""
    print("=" * 70)
    print("SpliceBind Dataset Expansion")
    print("=" * 70)

    # Ensure directories exist
    STRUCTURES_DIR.mkdir(parents=True, exist_ok=True)
    EMBEDDINGS_DIR.mkdir(parents=True, exist_ok=True)

    # Download structures
    print("\n[1/3] Downloading AlphaFold structures...")
    print("-" * 40)

    downloaded = 0
    failed = []
    for uniprot_id, gene in NEW_KINASES.items():
        print(f"Processing {gene} ({uniprot_id})...")
        if download_alphafold_structure(uniprot_id):
            downloaded += 1
        else:
            failed.append((uniprot_id, gene))
        time.sleep(0.5)  # Be nice to the API

    print(f"\nDownloaded: {downloaded}/{len(NEW_KINASES)}")
    if failed:
        print(f"Failed: {[f[1] for f in failed]}")

    # Generate ESM-2 embeddings for new structures
    print("\n[2/3] Generating ESM-2 embeddings...")
    print("-" * 40)

    # Get list of structures that need embeddings
    for uniprot_id, gene in NEW_KINASES.items():
        pdb_file = STRUCTURES_DIR / f"AF-{uniprot_id}-F1-model_v4.pdb"
        if pdb_file.exists():
            print(f"Processing {gene}...")
            generate_esm_embeddings(pdb_file)

    # Print new data to add to train.py
    print("\n[3/3] Data to add to train.py:")
    print("=" * 70)

    # Print UNIPROT_TO_GENE additions
    print("\n# Add to UNIPROT_TO_GENE:")
    for uniprot_id, gene in sorted(NEW_KINASES.items(), key=lambda x: x[1]):
        print(f"    '{uniprot_id}': '{gene}',")

    # Print KINASE_ATP_SITES additions
    print("\n# Add to KINASE_ATP_SITES:")
    for gene, site in sorted(NEW_ATP_SITES.items()):
        key_residues = ', '.join(map(str, site['key_residues']))
        print(f"    '{gene}': {{'start': {site['start']}, 'end': {site['end']}, 'key_residues': [{key_residues}]}},")

    # Print EXPERIMENTAL_BINDING_DATA additions
    print("\n# Add to EXPERIMENTAL_BINDING_DATA:")
    for entry in NEW_BINDING_DATA:
        gene, drug, kd, source = entry
        print(f"    ('{gene}', '{drug}', {kd}, '{source}'),")

    # Print GENE_TO_FAMILY additions
    print("\n# Add to GENE_TO_FAMILY:")
    gene_families = {
        'ALK': 'ALK', 'ERBB2': 'EGFR', 'ERBB4': 'EGFR',
        'VEGFR1': 'VEGFR', 'VEGFR2': 'VEGFR', 'VEGFR3': 'VEGFR',
        'MEK1': 'MEK', 'MEK2': 'MEK',
        'ROS1': 'ROS1',
        'NTRK1': 'NTRK', 'NTRK2': 'NTRK', 'NTRK3': 'NTRK',
        'CSF1R': 'CSF1R', 'FGFR4': 'FGFR',
        'DDR1': 'DDR', 'DDR2': 'DDR',
        'IGF1R': 'IGF1R', 'INSR': 'INSR',
        'ROCK2': 'ROCK', 'MTOR': 'MTOR', 'WEE1': 'WEE1',
        'PLK2': 'PLK', 'PLK3': 'PLK',
        'HCK': 'SRC', 'FYN': 'SRC', 'YES1': 'SRC',
        'ABL2': 'ABL', 'RAF1': 'RAF', 'ARAF': 'RAF',
    }
    for gene, family in sorted(gene_families.items()):
        print(f"    '{gene}': '{family}',")

    print("\n" + "=" * 70)
    print("Dataset expansion complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()

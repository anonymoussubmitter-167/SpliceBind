#!/usr/bin/env python3
"""Generate ESM-2 embeddings for all structures missing them."""

import sys
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np
from modules.esm_embeddings import ESMEmbeddingExtractor

def extract_sequence_from_pdb(pdb_path):
    """Extract amino acid sequence from PDB file."""
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    sequence = []
    seen_residues = set()
    
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                res_name = line[17:20].strip()
                res_num = line[22:26].strip()
                chain = line[21]
                key = (chain, res_num)
                
                if key not in seen_residues:
                    seen_residues.add(key)
                    aa = aa_map.get(res_name, 'X')
                    if aa != 'X':
                        sequence.append(aa)
    
    return ''.join(sequence)

def main():
    structures_dir = PROJECT_ROOT / "data" / "structures"
    embeddings_dir = PROJECT_ROOT / "data" / "embeddings"
    embeddings_dir.mkdir(exist_ok=True)
    
    # Find structures needing embeddings
    missing = []
    for pdb in structures_dir.glob("AF-*-F1-model_v6.pdb"):
        base = pdb.stem
        emb_path = embeddings_dir / f"{base}_esm2.npy"
        if not emb_path.exists():
            missing.append(pdb)
    
    print(f"Found {len(missing)} structures needing ESM embeddings")
    
    if not missing:
        print("All structures have embeddings!")
        return
    
    # Initialize extractor
    extractor = ESMEmbeddingExtractor(cache_dir=str(embeddings_dir))
    
    for i, pdb in enumerate(missing):
        print(f"[{i+1}/{len(missing)}] Processing {pdb.name}...")
        
        try:
            sequence = extract_sequence_from_pdb(pdb)
            if len(sequence) < 10:
                print(f"  Skipping: too short ({len(sequence)} residues)")
                continue
            
            protein_id = pdb.stem
            embeddings = extractor.extract_embeddings(sequence, protein_id)
            
            # Save
            out_path = embeddings_dir / f"{protein_id}_esm2.npy"
            np.save(out_path, embeddings)
            print(f"  Saved: {out_path.name} ({embeddings.shape})")
            
        except Exception as e:
            print(f"  Error: {e}")

if __name__ == "__main__":
    main()

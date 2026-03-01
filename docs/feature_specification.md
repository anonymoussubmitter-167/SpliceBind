# SpliceBind Node Feature Specification

Complete specification of all 56 node features used in PocketGNN.

## Overview

| Category | Features | Dimensions |
|----------|----------|------------|
| One-hot amino acid | 20 standard amino acids | 20 |
| Physicochemical | Hydrophobicity, charge, etc. | 4 |
| ESM-2 embeddings | Protein language model | 32 |
| **Total** | | **56** |

---

## 1. One-Hot Amino Acid Encoding (20 dimensions)

Standard amino acid identity, ordered alphabetically:

| Index | Amino Acid | 3-Letter | 1-Letter |
|-------|------------|----------|----------|
| 0 | Alanine | ALA | A |
| 1 | Arginine | ARG | R |
| 2 | Asparagine | ASN | N |
| 3 | Aspartic acid | ASP | D |
| 4 | Cysteine | CYS | C |
| 5 | Glutamic acid | GLU | E |
| 6 | Glutamine | GLN | Q |
| 7 | Glycine | GLY | G |
| 8 | Histidine | HIS | H |
| 9 | Isoleucine | ILE | I |
| 10 | Leucine | LEU | L |
| 11 | Lysine | LYS | K |
| 12 | Methionine | MET | M |
| 13 | Phenylalanine | PHE | F |
| 14 | Proline | PRO | P |
| 15 | Serine | SER | S |
| 16 | Threonine | THR | T |
| 17 | Tryptophan | TRP | W |
| 18 | Tyrosine | TYR | Y |
| 19 | Valine | VAL | V |

**Encoding**: Binary (0 or 1), exactly one position is 1 per residue.

---

## 2. Physicochemical Properties (4 dimensions)

### 2.1 Hydrophobicity Score (Index 20)
- **Scale**: Kyte-Doolittle hydropathy index
- **Range**: [-4.5, 4.5] (normalized to [0, 1])
- **Source**: Kyte & Doolittle, 1982

| Amino Acid | Raw Value | Normalized |
|------------|-----------|------------|
| I | 4.5 | 1.000 |
| V | 4.2 | 0.967 |
| L | 3.8 | 0.922 |
| F | 2.8 | 0.811 |
| C | 2.5 | 0.778 |
| M | 1.9 | 0.711 |
| A | 1.8 | 0.700 |
| G | -0.4 | 0.456 |
| T | -0.7 | 0.422 |
| S | -0.8 | 0.411 |
| W | -0.9 | 0.400 |
| Y | -1.3 | 0.356 |
| P | -1.6 | 0.322 |
| H | -3.2 | 0.144 |
| E | -3.5 | 0.111 |
| Q | -3.5 | 0.111 |
| D | -3.5 | 0.111 |
| N | -3.5 | 0.111 |
| K | -3.9 | 0.067 |
| R | -4.5 | 0.000 |

### 2.2 Net Charge (Index 21)
- **Scale**: Net charge at physiological pH (7.4)
- **Values**: {-1, 0, +1}
- **Encoding**: Mapped to [0, 0.5, 1]

| Category | Amino Acids | Charge | Encoded |
|----------|-------------|--------|---------|
| Negative | D, E | -1 | 0.0 |
| Neutral | A,C,F,G,H,I,L,M,N,P,Q,S,T,V,W,Y | 0 | 0.5 |
| Positive | K, R | +1 | 1.0 |

### 2.3 Polar/Non-polar (Index 22)
- **Binary classification**
- **Polar (1.0)**: S, T, N, Q, Y, H, K, R, D, E, C
- **Non-polar (0.0)**: A, V, L, I, M, F, W, P, G

### 2.4 Aromatic (Index 23)
- **Binary classification**
- **Aromatic (1.0)**: F, Y, W, H
- **Non-aromatic (0.0)**: All others

---

## 3. ESM-2 Embeddings (32 dimensions)

### Source
- **Model**: ESM-2 (esm2_t33_650M_UR50D)
- **Original dimension**: 1280
- **Reduction method**: Binned mean pooling

### Reduction Process
1. Extract per-residue embeddings from ESM-2 layer 33
2. Partition 1280 dimensions into 32 bins of 40 dimensions each
3. Average each bin to produce 32-dimensional embedding

```python
def reduce_esm_embedding(embedding_1280):
    """Reduce 1280-dim ESM-2 to 32-dim via binned pooling."""
    n_bins = 32
    bin_size = 40  # 1280 / 32
    reduced = np.zeros(n_bins)
    for i in range(n_bins):
        start = i * bin_size
        end = start + bin_size
        reduced[i] = embedding_1280[start:end].mean()
    return reduced
```

### Validation
- Cosine similarity preservation: 82.5% correlation with original pairwise similarities
- Comparison with PCA (32 components): PCA preserves 79.8% (binned pooling is better)

---

## 4. Feature Indices Summary

| Index Range | Feature Type | Description |
|-------------|--------------|-------------|
| 0-19 | One-hot AA | Amino acid identity |
| 20 | Hydrophobicity | Kyte-Doolittle scale, normalized |
| 21 | Charge | Net charge at pH 7.4 |
| 22 | Polarity | Binary polar/non-polar |
| 23 | Aromatic | Binary aromatic/non-aromatic |
| 24-55 | ESM-2 | Binned mean pooling of 1280-dim |

---

## 5. Edge Features (2 dimensions)

Edges connect residues within 6.0 Å (Cα-Cα distance).

| Index | Feature | Description |
|-------|---------|-------------|
| 0 | Distance | Euclidean distance, normalized to [0,1] by dividing by cutoff |
| 1 | Contact type | 1.0 if < 4Å (strong), 0.5 if 4-6Å (weak) |

---

## References

1. Kyte J, Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein. J Mol Biol 157:105-132.
2. Lin Z, et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379:1123-1130.

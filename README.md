# SpliceBind

A graph neural network for predicting kinase binding pocket druggability, with applications to splice variant drug resistance prediction.

## Overview

SpliceBind uses a GNN architecture to predict whether kinase binding pockets are druggable, trained on experimental binding data from BindingDB and ChEMBL. The model incorporates ESM-2 protein language model embeddings alongside physicochemical features.

**Key Result**: SpliceBind achieves AUROC 0.667, outperforming P2Rank (0.621) on experimental binding labels.

## Features

- **GNN-based pocket druggability prediction** using EdgeConv architecture
- **ESM-2 embeddings** for protein sequence features (+0.09 AUROC improvement)
- **P2Rank integration** for pocket detection
- **Splice variant analysis** for drug resistance prediction
- **78 kinase structures** across 35+ families

## Installation

```bash
# Clone repository
git clone https://github.com/anonymoussubmitter-167/SpliceBind.git
cd SpliceBind

# Install dependencies
pip install torch torch-geometric esm biopython requests

# Download P2Rank (included in tools/)
```

## Usage

### Training

```bash
python scripts/train.py
```

Results saved to `data/processed/results.json`

### Splice Variant Validation

```bash
# Generate kinase domain structures
python scripts/generate_variant_structures.py

# Validate predictions
python scripts/validate_kinase_domains.py
```

## Results

| Method | AUROC | Notes |
|--------|-------|-------|
| Random | 0.505 | Baseline |
| P2Rank | 0.621 | State-of-the-art pocket detector |
| **SpliceBind** | **0.667** | Our model |

### Splice Variant Case Studies

| Variant | Drug | Prediction | Status |
|---------|------|------------|--------|
| PIK3CD-S | Idelalisib | Reduced druggability | Correct |
| AR-V7 | Enzalutamide | No binding domain | Correct |
| MET-ex14skip | Capmatinib | Druggable | Correct |

See [results.md](results.md) for detailed results.

## Architecture

```
Input: Pocket residue graph
  |
  v
Node Features (56-dim)
  - 24 physicochemical properties
  - 32 ESM-2 embeddings
  |
  v
3x EdgeConv layers (256 hidden)
  |
  v
Mean + Max pooling
  |
  v
MLP classifier -> Druggability score
```

## Dataset

- **78 kinase structures** from AlphaFold
- **178 binding pockets** detected by P2Rank
- **150+ experimental binding records** from BindingDB/ChEMBL
- **Labels**: Kd < 100nM = potent, Kd > 10µM = weak

## Limitations

- Validated only on kinases
- Pocket-based methods cannot detect allosteric resistance mechanisms
- No wet-lab validation

## License

MIT License

## Citation

If you use SpliceBind, please cite this repository.

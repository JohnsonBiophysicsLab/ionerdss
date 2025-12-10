# Affinity Prediction with ProAffinity-GNN

This document describes how to use the ProAffinity-GNN integration for predicting protein-protein binding affinities in ionerdss.

## Overview

The `coarse_grain()` method now supports optional binding affinity prediction using ProAffinity-GNN, a graph neural network model trained on protein-protein complex structures.

## Quick Start

```python
from ionerdss.model.pdb_model import PDBModel

# Basic usage with affinity prediction
model = PDBModel(pdb_id='8erq', save_dir='./output')
model.coarse_grain(
    predict_affinity=True,
    adfr_path='/path/to/ADFRsuite/bin/prepare_receptor'
)
```

## Parameters

### `coarse_grain()` method

- **`predict_affinity`** (bool, default=False): Enable ProAffinity-GNN prediction
- **`adfr_path`** (str, optional): Path to ADFR `prepare_receptor` tool (required if `predict_affinity=True`)
- **`distance_cutoff`** (float, default=0.35): Max distance (nm) for interface detection
- **`residue_cutoff`** (int, default=3): Min residue pairs for valid interface
- **`standard_output`** (bool, default=False): Print detailed output

## Requirements

### ADFR Suite Installation

Download from: https://ccsb.scripps.edu/adfr/downloads/

```bash
# Example installation
wget https://ccsb.scripps.edu/adfr/download/1038/
tar -xzvf ADFRsuite_x86_64Linux_1.0.tar.gz
cd ADFRsuite_x86_64Linux_1.0
./install.sh

# Set ADFR_PATH environment variable
export ADFR_PATH="/path/to/ADFRsuite/bin/prepare_receptor"
```

### Python Dependencies

```bash
pip install torch torch_geometric transformers
```

See `ionerdss/model/proaffinity_requirements.txt` for specific versions.

## Energy Values

- **With ProAffinity**: Predicted binding energy in kJ/mol
- **Without ProAffinity** (default): -39.5 kJ/mol (-16 RT at 298K)
- **Fallback**: Uses default value if prediction fails

## Example Output

```
Chain-Chain Interaction: A-B
    Interface 1: [1, 2, 3, 4, 5]
    Interface 2: [10, 11, 12, 13, 14]
    Interface Energy: -45.23 kJ/mol
    Predicted binding energy for chains A-B: -45.23 kJ/mol
```

## Error Handling

The system automatically falls back to default energy if:
- ProAffinity prediction fails
- ADFR tools are not available
- PDB conversion errors occur

## Performance Notes

- ProAffinity prediction adds ~10-30 seconds per interface
- Only runs for valid interfaces (residue_cutoff threshold met)
- Predictions are cached during the same `coarse_grain()` call

## Advanced Usage

### Using Pre-downloaded PDB Files

```python
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy_from_file

energy = predict_proaffinity_binding_energy_from_file(
    pdb_file='./structures/8erq.pdb',
    chains='A,B',
    adfr_path='/path/to/prepare_receptor',
    verbose=True
)
print(f"Binding energy: {energy:.2f} kJ/mol")
```

### Custom Temperature

ProAffinity model returns K_d, converted to ΔG using:

```
ΔG = -RT ln(K_d)
```

Default temperature is 298.15 K.

## Troubleshooting

### "ADFR path not provided"
Set the `adfr_path` parameter or `ADFR_PATH` environment variable.

### "ProAffinity prediction failed"
Check that:
1. ADFR tools are installed and accessible
2. PDB file contains both specified chains
3. Chains have sufficient interface residues

### Import errors
Ensure ProAffinity dependencies are installed:
```bash
pip install -r ionerdss/model/proaffinity_requirements.txt
```

## References

ProAffinity-GNN uses ESM2 protein language model embeddings and graph neural networks to predict binding affinities. See the original ProAffinity paper for methodology details.

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

## If you are on a mac, you can use the following command to install ADFR to bypass macOS marking python2 as "untrusted developer":

chmod +x ./examples/install_ADFR_mac.sh
./examples/install_ADFR_mac.sh

# Set ADFR_PATH environment variable
# A script cannot permanently modify your shell’s PATH just by echoing export PATH=... inside itself
# Because each script runs in its own subshell, and environment changes do not propagate back to your interactive terminal.
export ADFR_PATH="/path/to/ADFRsuite/bin/prepare_receptor"
```

### Python Dependencies

It is suggested that the user create a new environment separately because the version with proaffinity enabled uses an earlier version of `numpy` and `scipy`.

```bash
pip install "ionerdss[proaffinity]"
```

See `pyproject.toml` for specific information about dependencies and their versions.

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

- ProAffinity prediction adds ~30-90 seconds per interface (hardware dependent)
- Only runs for valid interfaces (residue_cutoff threshold met)
- Predictions are cached during the same `coarse_grain()` call

## Custom Temperature

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
2. PDB file contains the specified chains
3. Chains have sufficient interface residues

## References

Zhiyuan Zhou, Yueming Yin, Hao Han, Yiping Jia, Jun Hong Koh, Adams Wai-Kin Kong, Yuguang Mu. ProAffinity-GNN: A Novel Approach to Structure-Based Protein–Protein Binding Affinity Prediction via a Curated Data Set and Graph Neural Networks. J. Chem. Inf. Model. 2024, 64, 23, 8796–8808. https://doi.org/10.1021/acs.jcim.4c01850


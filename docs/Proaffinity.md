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

ProAffinity pins numpy 1.x and torch 2.2. Those cannot share an environment with anything that needs numpy 2 -- the OVITO renderer, in particular -- so it belongs in an environment of its own:

```bash
pip install "ionerdss[proaffinity]"
```

See `pyproject.toml` for specific information about dependencies and their versions.

### Running ProAffinity from another environment

You do not have to work inside the ProAffinity environment. ioNERDSS can call into it as a sidecar: only a PDB path, chain pairs, and the resulting energies cross the boundary as JSON, so the two dependency stacks never meet.

Create the sidecar once:

```bash
python -m venv ~/.ionerdss-proaffinity
~/.ionerdss-proaffinity/bin/pip install "ioNERDSS[proaffinity]"
export IONERDSS_PROAFFINITY_PYTHON=~/.ionerdss-proaffinity/bin/python
```

With that variable set, `predict_affinity=True` works from your main environment unchanged. To point at it explicitly instead of using the environment variable:

```python
build_system_from_pdb(
    source='8erq',
    predict_affinity=True,
    proaffinity_python='~/.ionerdss-proaffinity/bin/python',
)
```

Two hyperparameters control this:

- **`proaffinity_backend`** (str, default=`'auto'`): `'auto'` uses the sidecar when one is configured and runs in-process otherwise; `'sidecar'` requires one; `'in_process'` never spawns one.
- **`proaffinity_python`** (str, optional): the sidecar interpreter. Defaults to `$IONERDSS_PROAFFINITY_PYTHON`.

The sidecar needs ProAffinity's dependencies, not a second ioNERDSS install -- the worker imports ioNERDSS from the source tree it ships with. ADFR still has to be reachable from the sidecar, so set `ADFR_PATH` in that environment or pass `adfr_path`.

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


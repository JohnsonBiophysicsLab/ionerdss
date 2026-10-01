# Affinity Prediction with ProAffinity-GNN

This document describes how to use the ProAffinity-GNN integration for predicting protein-protein binding affinities in ionerdss.

## Overview

`build_system_from_pdb()` supports optional binding affinity prediction using ProAffinity-GNN, a graph neural network model trained on protein-protein complex structures. The predicted energy of each detected interface replaces the default binding energy.

## Quick Start

```python
from ionerdss import build_system_from_pdb

# Basic usage with affinity prediction
system = build_system_from_pdb(
    source='8erq',
    workspace_path='./output',
    predict_affinity=True,
    adfr_path='/path/to/ADFRsuite/bin/prepare_receptor'
)
```

## Parameters

### `build_system_from_pdb()` hyperparameters

- **`predict_affinity`** (bool, default=False): Enable ProAffinity-GNN prediction
- **`adfr_path`** (str, optional): Path to ADFR `prepare_receptor` tool, or to the ADFR install directory or its `bin/`. Defaults to `$ADFR_PATH`; one of the two is required if `predict_affinity=True`
- **`interface_detect_distance_cutoff`** (float, default=0.9): Max distance (nm) for interface detection
- **`interface_detect_n_residue_cutoff`** (int, default=2): Min contacting residues on each chain for a valid interface

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

chmod +x ./tutorials/install_ADFR_mac.sh
./tutorials/install_ADFR_mac.sh

# Set ADFR_PATH environment variable
# A script cannot permanently modify your shell’s PATH just by echoing export PATH=... inside itself
# Because each script runs in its own subshell, and environment changes do not propagate back to your interactive terminal.
export ADFR_PATH="/path/to/ADFRsuite/bin/prepare_receptor"
```

### Python Dependencies

ProAffinity pins numpy 1.x and torch 2.2. Those cannot share an environment with anything that needs numpy 2 -- OVITO 3.16 and later, in particular -- so it belongs in an environment of its own:

```bash
pip install "ionerdss[proaffinity]"
```

See `pyproject.toml` for specific information about dependencies and their versions.

### Running ProAffinity from another environment

You do not have to work inside the ProAffinity environment. ioNERDSS can call into it as a sidecar: only a PDB path, chain pairs, and the resulting energies cross the boundary as JSON, so the two dependency stacks never meet.

Create the sidecar once. With conda (recommended -- it picks the Python version for you, and torch 2.2.2 publishes cp38-cp312 wheels only):

```bash
conda env create -f env/proaffinity/environment.yml
export IONERDSS_PROAFFINITY_PYTHON="$(conda info --base)/envs/ionerdss-proaffinity/bin/python"
```

Or without a checkout of this repository:

```bash
conda create -y -n ionerdss-proaffinity python=3.10
conda run -n ionerdss-proaffinity pip install "ioNERDSS[proaffinity]"
export IONERDSS_PROAFFINITY_PYTHON="$(conda info --base)/envs/ionerdss-proaffinity/bin/python"
```

Or with `venv`, which inherits the Python version of the interpreter you run it from -- fine on 3.10-3.12, impossible on 3.13+:

```bash
python -m venv ~/.ionerdss-proaffinity
~/.ionerdss-proaffinity/bin/pip install "ioNERDSS[proaffinity]"
export IONERDSS_PROAFFINITY_PYTHON=~/.ionerdss-proaffinity/bin/python
```

On a cluster, create the environment somewhere with room for several GB of torch and model weights, and where a compute node can read it:

```bash
conda create -y -p /scratch/$USER/envs/ionerdss-proaffinity python=3.10
conda run -p /scratch/$USER/envs/ionerdss-proaffinity pip install "ioNERDSS[proaffinity]"
export IONERDSS_PROAFFINITY_PYTHON=/scratch/$USER/envs/ionerdss-proaffinity/bin/python
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

The sidecar needs ProAffinity's dependencies and ioNERDSS's own (`ioNERDSS[proaffinity]` brings both), but not a matching ioNERDSS version -- the worker imports ioNERDSS from the source tree it ships with. ADFR still has to be reachable from the sidecar, which inherits the environment variables of the process that launches it, so set `ADFR_PATH` where you run ioNERDSS or pass `adfr_path`.

## Energy Values

- **With ProAffinity**: Predicted binding energy in kJ/mol
- **Without ProAffinity** (default): -39.6 kJ/mol (-16 RT at 298K)
- **Fallback**: Uses default value if prediction fails

## Example Output

Logged at INFO level for each interface:

```
Predicted energy for A-B: -45.23 kJ/mol
```

A pair that could not be predicted logs `ProAffinity prediction failed for A-B, using default energy` instead.

## Error Handling

The system automatically falls back to default energy if:
- ProAffinity prediction fails
- ADFR tools are not available
- PDB conversion errors occur

## Performance Notes

- ProAffinity prediction adds ~30-90 seconds per interface (hardware dependent)
- Only runs for valid interfaces (`interface_detect_n_residue_cutoff` threshold met)
- The ESM-2 model is loaded once and reused for every interface, and the structure's PDBQT conversion is reused once it exists

## Energy Conversion

ProAffinity model returns pK_d (= -log10 K_d), converted to ΔG using:

```
ΔG = -RT ln(10^pK_d) = RT ln(K_d)
```

The temperature is fixed at 298.15 K.

## Troubleshooting

### "ADFR_PATH environment variable not set" or "prepare_receptor not found"
Set the `adfr_path` parameter or `ADFR_PATH` environment variable.

### "ProAffinity prediction failed"
Check that:
1. ADFR tools are installed and accessible
2. PDB file contains the specified chains
3. Chains have sufficient interface residues

## References

Zhiyuan Zhou, Yueming Yin, Hao Han, Yiping Jia, Jun Hong Koh, Adams Wai-Kin Kong, Yuguang Mu. ProAffinity-GNN: A Novel Approach to Structure-Based Protein–Protein Binding Affinity Prediction via a Curated Data Set and Graph Neural Networks. J. Chem. Inf. Model. 2024, 64, 23, 8796–8808. https://doi.org/10.1021/acs.jcim.4c01850


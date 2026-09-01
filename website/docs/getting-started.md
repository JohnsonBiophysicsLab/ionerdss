# Getting Started

## Install from PyPI

```bash
pip install ioNERDSS
```

Optional extras:

```bash
pip install "ioNERDSS[jupyter]"
pip install "ioNERDSS[ovito_rendering]"
```

On clusters running RHEL/Rocky 8 -- Rockfish among them -- pip installs OVITO 3.15 instead of 3.16: 3.16 pins a PySide6 whose Linux wheels need glibc 2.34, newer than those systems ship. Rendering still works there; ioNERDSS falls back to OVITO's software ray tracer, which is slower than the OpenGL path but needs no display. If the install instead fails with `No matching distribution found for PySide6`, you are on ioNERDSS 2.2.3 or earlier -- `pip install -U ioNERDSS` and retry.

ProAffinity-GNN, which predicts protein-protein binding free energies, is installed differently -- see below.

## Development install

```bash
git clone https://github.com/JohnsonBiophysicsLab/ionerdss.git
cd ionerdss
pip install -e ".[test,jupyter]"
pytest
```

## Binding affinity prediction (ProAffinity-GNN)

ProAffinity-GNN pins numpy 1.x and torch 2.2. Those cannot share an environment with the numpy 2 that OVITO rendering requires, and installing them together leaves one of the two quietly broken. So ProAffinity gets an environment of its own, and ioNERDSS calls into it when it needs a binding energy: only a PDB path, the chain pairs, and the resulting energies cross between the two.

Set it up once:

```bash
conda env create -f env/proaffinity/environment.yml
export IONERDSS_PROAFFINITY_PYTHON="$(conda info --base)/envs/ionerdss-proaffinity/bin/python"
```

Without a checkout of the repository:

```bash
conda create -y -n ionerdss-proaffinity python=3.10
conda run -n ionerdss-proaffinity pip install "ioNERDSS[proaffinity]"
export IONERDSS_PROAFFINITY_PYTHON="$(conda info --base)/envs/ionerdss-proaffinity/bin/python"
```

Python is pinned because torch 2.2.2 publishes wheels for CPython 3.8-3.12 only. Add the `export` to your shell profile to make it persist.

You never activate that environment yourself. With the variable set, affinity prediction works from your normal environment:

```python
system = ion.build_system_from_pdb(
    source="8y7s",
    predict_affinity=True,
    adfr_path="~/Documents/ADFR",
)
```

Prediction also needs the ADFR suite, which has to be reachable from the ProAffinity environment. Two hyperparameters control where prediction runs: `proaffinity_backend` (`auto`, `sidecar`, or `in_process`) and `proaffinity_python`, which overrides the environment variable.

## First model build

```python
import ionerdss as ion

system = ion.build_system_from_pdb(
    source="5l93",
    workspace_path="5l93_dir",
    generate_visualizations=True,
    generate_nerdss_files=True,
)
```

The generated workspace can include:

- `structures/` for fetched or normalized structures
- `visualizations/` for coarse-grained outputs
- `nerdss_files/` for `.mol` files and `parms.inp`
- `ode_results/` when the ODE pipeline is enabled
- `logs/` for pipeline execution details

## Local docs preview

```bash
pip install -e ".[docs]"
pip install -r website/requirements.txt
mkdocs serve -f website/mkdocs.yml
```

Open `http://127.0.0.1:8000/`.

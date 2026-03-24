# Getting Started

## Install from PyPI

```bash
pip install ioNERDSS
```

Optional extras:

```bash
pip install "ioNERDSS[ode,ovito_rendering]"
pip install "ioNERDSS[jupyter]"
```

## Development install

```bash
git clone https://github.com/JohnsonBiophysicsLab/ionerdss.git
cd ionerdss
pip install -e ".[test,jupyter]"
pytest
```

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

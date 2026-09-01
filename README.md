[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://johnsonbiophysicslab.github.io/ionerdss/)
[![Run Unit Tests](https://github.com/JohnsonBiophysicsLab/ionerdss/actions/workflows/unittest.yml/badge.svg?branch=main&event=push)](https://github.com/JohnsonBiophysicsLab/ionerdss/actions/workflows/unittest.yml)
![PyPI](https://img.shields.io/pypi/v/ioNERDSS.svg)
![PyPI - Downloads](https://img.shields.io/pypi/dm/ioNERDSS.svg)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/ioNERDSS.svg)


![ioNERDSS Banner](https://raw.githubusercontent.com/JohnsonBiophysicsLab/ionerdss/main/website/assets/banner.png)


> ### Try NERDSS / ioNERDSS online  
> Try NERDSS / ioNERDSS on the webserver at **nerdssdemo.org** without local installation.  
>[![Webserver](https://img.shields.io/badge/Webserver-nerdssdemo.org-blue?style=for-the-badge&logo=internet-explorer)](https://nerdssdemo.org)


[![Documentation](https://img.shields.io/badge/Documentation-ioNERDSS%20site-lightgrey?style=for-the-badge&logo=readthedocs)](https://johnsonbiophysicslab.github.io/ionerdss/)
[![Article](https://img.shields.io/badge/Article-bioRxiv-orange?style=for-the-badge&logo=readme)](https://www.biorxiv.org/content/10.64898/2026.01.27.702082v1)
[![PyPI](https://img.shields.io/badge/PyPI-ioNERDSS-blue?style=for-the-badge&logo=pypi)](https://pypi.org/project/ioNERDSS/)
[![GitHub](https://img.shields.io/badge/GitHub-Source-black?style=for-the-badge&logo=github)](https://github.com/JohnsonBiophysicsLab/ionerdss)
[![Learn More](https://img.shields.io/badge/Learn%20More-Website-9cf?style=for-the-badge&logo=info)](https://johnsonbiophysicslab.github.io/NERDSS/)

**ionerdss** is a Python library for building NERDSS-ready models from structures, running simulation workflows, and analyzing simulation outputs.

## Installation
Recommended: Install in an isolated environment
We strongly recommend installing ioNERDSS in an isolated environment using a package manager like conda or mamba to avoid conflicts with existing Python libraries.

### Using Conda
```bash
# Create a new environment with Python 3.10 or later
conda create -n ionerdss python=3.10
conda activate ionerdss
```

### 1. From PyPI

Requires Python 3.10 or later.

```bash
pip install ioNERDSS
```

To include optional features:

```bash
pip install "ioNERDSS[jupyter]"
pip install "ioNERDSS[ovito_rendering]"
```

On an HPC cluster, read step 4 before installing: the OVITO extras resolve differently there, and the environment belongs on scratch rather than in your home directory.

The ODE pipeline does not require a separate extra; it is included in the main package install.

`proaffinity` is the exception: it needs an environment of its own, set up in step 3 below.

### 2. For development

```bash
git clone https://github.com/JohnsonBiophysicsLab/ionerdss.git
cd ionerdss
pip install -e ".[test,jupyter]"
pytest
```

For the full contributor environment:

```bash
pip install -e ".[all]"
```

`[all]` covers everything except `proaffinity`, which needs its own environment (step 3).

### 3. Binding affinity prediction (ProAffinity-GNN)

ProAffinity-GNN pins numpy 1.x and torch 2.2, which cannot share an environment with the numpy 2 that OVITO rendering needs. Installing it alongside the rest quietly breaks one of the two, so it goes in an environment of its own and ioNERDSS calls into that environment when it needs a binding energy. Only a PDB path and the resulting energies cross between them.

Set it up once:

```bash
conda env create -f env/proaffinity/environment.yml
export IONERDSS_PROAFFINITY_PYTHON="$(conda info --base)/envs/ionerdss-proaffinity/bin/python"
```

Without a checkout of this repository:

```bash
conda create -y -n ionerdss-proaffinity python=3.10
conda run -n ionerdss-proaffinity pip install "ioNERDSS[proaffinity]"
export IONERDSS_PROAFFINITY_PYTHON="$(conda info --base)/envs/ionerdss-proaffinity/bin/python"
```

Add that `export` to your shell profile to make it stick. You never activate the environment yourself -- with the variable set, `predict_affinity=True` works from your normal environment:

```python
system = build_system_from_pdb(source="8erq", predict_affinity=True)
```

Affinity prediction also needs the ADFR suite. See [docs/Proaffinity.md](docs/Proaffinity.md) for that, for the `venv` and cluster variants, and for pointing at the environment explicitly instead of using the variable.

### 4. On an HPC cluster (Rockfish and other Slurm systems)

Nothing below needs root or a GPU, and all of it runs on a login node.

**Build the environment on scratch, not in your home directory.** `[all]` is several GB once OVITO and Qt land, and the ProAffinity sidecar adds several more in torch and model weights -- more than most home quotas allow. Send pip's cache and temporary directory there too, or the build will fill your home quota even when the environment itself does not live there.

```bash
module load anaconda                      # whatever `module avail conda` offers
export PIP_CACHE_DIR=/scratch/$USER/.cache/pip
export TMPDIR=/scratch/$USER/tmp && mkdir -p "$TMPDIR"

conda create -y -p /scratch/$USER/envs/ionerdss python=3.12
conda activate /scratch/$USER/envs/ionerdss
pip install "ioNERDSS[all]"
```

Substitute your cluster's scratch path -- on Rockfish that is your `~/scratch4-<PI>` or `~/data-<PI>` space. Python 3.10 through 3.13 all work for this environment; the ProAffinity sidecar below is the one piece that needs 3.10-3.12.

**OVITO resolves to 3.15 on RHEL/Rocky 8 clusters, by design.** OVITO 3.16 pins a PySide6 whose Linux wheels require glibc 2.34, and RHEL/Rocky 8 -- Rockfish included -- ships glibc 2.28. ioNERDSS therefore asks for `ovito>=3.15` rather than `>=3.16`, so pip installs OVITO 3.15.5 with PySide6 6.9.3 there and takes 3.16 on newer systems. To see which side of that line you are on:

```bash
python -c "import platform; print(platform.libc_ver())"
```

If that prints 2.34 or newer you get OVITO 3.16 and its OpenGL renderer; below that you get 3.15 and the software ray tracer. Both render without a display.

If the install fails with `No matching distribution found for PySide6`, you are on ioNERDSS 2.2.3 or earlier, which floors OVITO at 3.16 and so cannot resolve on glibc 2.28 at all:

```bash
pip install -U "ioNERDSS[all]"
```

If you have to stay on an older ioNERDSS, install the rendering extras by hand instead: `pip install "ovito<3.16" imageio Pillow`.

**Rendering needs no display, but does want cores.** `visualize_trajectory_ovito` renders in a child process with no X11 connection, so it works on a compute node as it stands. Where OVITO 3.15 is what resolved, the frames come from a CPU ray tracer, so give the job a few cores and expect it to be slower than a desktop OpenGL render:

```bash
srun -n 1 -c 8 --pty bash        # or the equivalent in your sbatch script
```

**Notebooks.** Register the environment as a Jupyter kernel so your cluster's notebook portal can find it:

```bash
python -m ipykernel install --user --name ionerdss --display-name "ioNERDSS"
```

**ProAffinity.** The sidecar environment belongs on scratch as well, for the same quota reason:

```bash
conda create -y -p /scratch/$USER/envs/ionerdss-proaffinity python=3.10
conda run -p /scratch/$USER/envs/ionerdss-proaffinity pip install "ioNERDSS[proaffinity]"
export IONERDSS_PROAFFINITY_PYTHON=/scratch/$USER/envs/ionerdss-proaffinity/bin/python
```

Put that `export` in your `~/.bashrc` and in any sbatch script, since batch jobs do not inherit your login shell's environment. See [docs/Proaffinity.md](docs/Proaffinity.md) for the ADFR suite, which prediction also requires.

To run simulations for structures generated from ioNERDSS you will also need to install NERDSS.

NERDSS installation instructions - [NERDSS Github](https://github.com/JohnsonBiophysicsLab/NERDSS)
## Usage

```python
import ionerdss as ion

system = ion.build_system_from_pdb(
    source="6bno",
    workspace_path="6bno_dir",
    ode_enabled=True,
)

analyzer = ion.Analyzer("path/to/simulation/root")
analyzer.plot.free_energy()
```

## Tutorials

Current examples live under `tutorials/` and are the supported starting point for notebook-based workflows:

- `tutorials/quick_start_6bno.ipynb`
- `tutorials/ionerdss_tutorial_5l93.ipynb`
- `tutorials/ionerdss_tutorial_6bno.ipynb`
- `tutorials/ionerdss_tutorial_8y7s.ipynb`
- `tutorials/additional_examples/4yd9.ipynb`

Open them locally with Jupyter:

```bash
git clone https://github.com/JohnsonBiophysicsLab/ionerdss.git
cd ionerdss
pip install -e ".[jupyter]"
jupyter lab tutorials/
```

Hosted documentation and tutorial index:

- Docs: [johnsonbiophysicslab.github.io/ionerdss](https://johnsonbiophysicslab.github.io/ionerdss/)
- Tutorials page: [Tutorials](https://johnsonbiophysicslab.github.io/ionerdss/tutorials/)

### Run a quick trial with our server

Go to the [NERDSS server](http://52.15.142.249:5000/).

## Documentation

The documentation site is published with GitHub Pages from the `website/` folder.

To preview it locally:

```bash
pip install -e ".[docs]"
mkdocs serve -f website/mkdocs.yml
```

Then open `http://127.0.0.1:8000/`.

## Running tests

```bash
pytest
```

## Docker development environment

```bash
docker build --no-cache -t ionerdss_dev .
docker run -it --rm -v $(pwd):/app -p 8888:8888 ionerdss_dev
```

This creates a containerized environment with Jupyter Lab accessible at `http://localhost:8888`.

## License

This project is licensed under the GPL-3.0 License.

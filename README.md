![ioNERDSS Banner](https://raw.githubusercontent.com/JohnsonBiophysicsLab/ionerdss/main/website/assets/banner.png)

---

[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://johnsonbiophysicslab.github.io/ionerdss/)
[![Run Unit Tests](https://github.com/JohnsonBiophysicsLab/ionerdss/actions/workflows/unittest.yml/badge.svg?branch=main&event=push)](https://github.com/JohnsonBiophysicsLab/ionerdss/actions/workflows/unittest.yml)
![PyPI](https://img.shields.io/pypi/v/ioNERDSS.svg)
![PyPI - Downloads](https://img.shields.io/pypi/dm/ioNERDSS.svg)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/ioNERDSS.svg)

---

> ### Try NERDSS / ioNERDSS online  
> Try NERDSS / ioNERDSS on the webserver at **nerdssdemo.org** without local installation.  
>
> [![Webserver](https://img.shields.io/badge/Webserver-nerdssdemo.org-blue?style=for-the-badge&logo=internet-explorer)](https://nerdssdemo.org)
> [![Documentation](https://img.shields.io/badge/Documentation-NERDSS%20site-lightgrey?style=for-the-badge&logo=readthedocs)](https://johnsonbiophysicslab.github.io/NERDSS/)
> [![Article](https://img.shields.io/badge/Article-bioRxiv-orange?style=for-the-badge&logo=readme)](https://www.biorxiv.org/content/10.64898/2026.01.27.702082v1)
> [![GitHub](https://img.shields.io/badge/GitHub-Repository-black?style=for-the-badge&logo=github)](https://github.com/JohnsonBiophysicsLab/NERDSS)
> [![PyPI](https://img.shields.io/badge/PyPI-ioNERDSS-blue?style=for-the-badge&logo=pypi)](https://pypi.org/project/ioNERDSS/)
> [![Learn More](https://img.shields.io/badge/Learn%20More-Website-9cf?style=for-the-badge&logo=info)](https://johnsonbiophysicslab.github.io/NERDSS/)

**ionerdss** is a Python library for building NERDSS-ready models from structures, running simulation workflows, and analyzing simulation outputs.

## Installation

### 1. From PyPI

Requires Python 3.10 or later.

```bash
pip install ioNERDSS
```

To include optional features:

```bash
pip install "ioNERDSS[ode,ovito_rendering]"
```

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
pip install -r website/requirements.txt
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

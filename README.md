# ionerdss
[![Documentation Status](https://readthedocs.org/projects/ionerdss/badge/?version=latest)](https://ionerdss.readthedocs.io/en/latest/?badge=latest)
[![Run Unit Tests](https://github.com/JohnsonBiophysicsLab/ionerdss/actions/workflows/unittest.yml/badge.svg?branch=main&event=push)](https://github.com/JohnsonBiophysicsLab/ionerdss/actions/workflows/unittest.yml)
![PyPI](https://img.shields.io/pypi/v/ioNERDSS.svg)
![PyPI - Downloads](https://img.shields.io/pypi/dm/ioNERDSS.svg)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/ioNERDSS.svg)
[![codecov](https://codecov.io/gh/JohnsonBiophysicsLab/ionerdss/graph/badge.svg?token=IUUUOZT0VJ)](https://codecov.io/gh/JohnsonBiophysicsLab/ionerdss)

**ionerdss** is a Python library that provides user‐friendly tools for setting up and analyzing output from the [NERDSS](https://github.com/JohnsonBiophysicsLab/NERDSS) reaction‐diffusion simulator. Its goal is to streamline model building (from PDB files or from scratch), data analysis, and visualization for simulation workflows.

---

## Installation

1. From PyPI (recommended):
   - **Python version:** 3.9 or later

Create a conda environment (optional but recommended):

Download and install Anaconda or Miniconda, then create a new conda environment for `ionerdss`:

```bash
conda create -n ionerdss python=3.9
conda activate ionerdss
```

Install the latest release directly from PyPI:

```bash
pip install ioNERDSS
```

2. From GitHub (for the latest development version):
   - If you want to use the latest features or contribute to the development, you can install directly from the GitHub repository:

To install from source (e.g., if you've cloned this repo and want the most recent changes):

```bash
git clone https://github.com/JohnsonBiophysicsLab/ionerdss.git
cd ionerdss

# Modern approach (recommended)
pip install -e .  # Editable mode: updates reflect immediately in the environment

# Alternative: Legacy approach using environment files
pip install -r env/requirements.txt
pip install -e .
```

---

## Quick Start

```python
import ionerdss as ion
ion.some_function()  # Replace with actual function calls to set up model and analyze results
```

For extended examples, see the [tutorials](https://ionerdss.readthedocs.io/en/latest/ionerdss_tutorials.html).

### Run a quick trial with our server

Go to the [NERDSS server](http://52.15.142.249:5000/).

---

## Documentation
- **User Guide:** [ionerdss user guide](https://ionerdss.readthedocs.io/en/latest/ionerdss_documentation_v1_1.html).

- **API Reference:** [API](https://ionerdss.readthedocs.io/en/latest/ionerdss.html). You can also build the docs locally using Sphinx:
```bash
sphinx-apidoc -o docs/source ionerdss
cd docs
make html
```
Then open docs/build/html/index.html in your browser.

---

## Repository Structure
```
ionerdss/
├── .github/workflows/     # Continuous Integration workflows
├── docs/                  # Documentation
│   ├── source/            # Sphinx source files
│   ├── make.bat           # Windows build script
│   └── Makefile           # Unix build script
├── ionerdss/              # Main Python package
│   ├── model/      # Model building tools (v1.2.0)
│   ├── nerdss_simulation/ # Simulation tools (v1.2.0)
│   ├── nerdss_analysis/   # Data analysis tools (v1.2.0)
│   └── __init__.py 
├── tests/                 # Unit tests
├── data/                  # Test and tutorial data
└── setup.py               # Installation & packaging
```

---

## Develop using docker container:  
```bash
docker build --no-cache -t ionerdss_dev . 
docker run -it --rm -v $(pwd):/app -p 8888:8888 ionerdss_dev
```

---

## Best Practices

1. **Docstrings & Sphinx**  
   - Write clear docstrings in Google‐style to help auto‐generate documentation.

2. **Code Organization**  
   - Keep related functionality grouped in submodules.

3. **Development Setup and Tests**  
   - Add or update unit tests in `tests/` for any new function. We use [unittest](https://docs.python.org/3/library/unittest.html).

   - **Recommended setup for development and testing:**
     ```bash
     # Create conda environment (recommended for isolation)
     conda create -n ionerdss python=3.9
     conda activate ionerdss
     
     # Install in development mode with test dependencies
     pip install -e ".[tests]"
     
     # Run tests
     pytest
     ```

   - **Alternative setup using legacy environment files:**
     ```bash
     # Using conda
     conda env create -f env/environment-full.yml
     conda activate ionerdss_env
     
     # OR using pip
     pip install -r env/requirements-full.txt
     
     # Install package in development mode
     pip install -e .
     ```

4. **Versioning & Releases**  
   - Update `setup.py` with a new version number. A GitHub release will auto‐update the PyPI package.

5. **Contributions**  
   - Fork the repo, create a feature branch, and open a pull request.

---

## License
This project is licensed under the GPL‐3.0 License.

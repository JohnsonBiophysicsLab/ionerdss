# ionerdss
[![Documentation Status](https://readthedocs.org/projects/ionerdss/badge/?version=latest)](https://ionerdss.readthedocs.io/en/latest/?badge=latest)
[![Run Unit Tests](https://github.com/JohnsonBiophysicsLab/ionerdss/actions/workflows/unittest.yml/badge.svg?branch=main&event=push)](https://github.com/JohnsonBiophysicsLab/ionerdss/actions/workflows/unittest.yml)
![PyPI](https://img.shields.io/pypi/v/ioNERDSS.svg)
![PyPI - Downloads](https://img.shields.io/pypi/dm/ioNERDSS.svg)
![PyPI - Wheel](https://img.shields.io/pypi/wheel/ioNERDSS.svg)
[![codecov](https://codecov.io/gh/JohnsonBiophysicsLab/ionerdss/graph/badge.svg?token=IUUUOZT0VJ)](https://codecov.io/gh/JohnsonBiophysicsLab/ionerdss)

**ionerdss** is a Python library that provides user‐friendly tools for setting up and analyzing output from the [NERDSS](https://github.com/JohnsonBiophysicsLab/NERDSS) reaction‐diffusion simulator.

---

## Installation

### 1. From PyPI (Recommended)

Requires Python 3.9 or later.

```bash
pip install ioNERDSS
```
To include optional features, you can specify them during installation:
```bash
# Example: Install with ODE solver and rendering tools
pip install "ioNERDSS[ode,ovito_rendering]"
```

### 2. For Development (from GitHub)

**Quick Start (Clone and Test):**
```bash
# Requires Python 3.9+
git clone https://github.com/JohnsonBiophysicsLab/ionerdss.git
cd ionerdss
pip install -e ".[full]"
pytest
```

**Detailed Setup for Contributors:**

If you want to contribute to development, work with examples, or need a specific environment setup:

**Prerequisites:**
*   [Git](https://git-scm.com/)
*   Python 3.9+
*   Your choice of environment manager: `conda`, `venv`, etc.
*   Optionally, [uv](https://github.com/astral-sh/uv) for faster performance

**Setup Instructions:**

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/JohnsonBiophysicsLab/ionerdss.git
    cd ionerdss
    ```

2.  **Create and activate an environment:**

    *   **Using `conda` (Recommended for full environment):**
        ```bash
        conda env create -f env/environment.yml
        conda activate ionerdss-dev
        ```
    *   **Using `uv` and `venv` (Fastest):**
        ```bash
        python -m venv .venv
        source .venv/bin/activate  # On Windows: .venv\Scripts\activate
        uv pip install -e ".[full]"
        ```
    *   **Using `pip` and `venv`:**
        ```bash
        python -m venv .venv
        source .venv/bin/activate  # On Windows: .venv\Scripts\activate
        pip install -e ".[all]"
        ```

This installs `ionerdss` in editable mode with all dependencies for development, testing, and running examples.

---

## Running Tests

To run the unit tests locally, ensure you have installed the development environment (which includes `pytest` and `pytest-cov`).

Then, from the project root folder, run:
```bash
pytest
```

---

## Usage

```python
import ionerdss as ion

# Create a PDB model from structure
model = ion.PDBModel(pdb_id="1abc")
model.generate_model()

# Run analysis on simulation data
analysis = ion.Analysis("path/to/simulation/data")
analysis.plot_time_evolution()
```

## Examples

The `examples/` directory contains hands-on Jupyter notebooks demonstrating real molecular systems:

- **`Homo-3mer-5VA4.ipynb`** - 3-component homogeneous assembly
- **`Hetero-30mer-4YD9.ipynb`** - 30-component heterogeneous system  
- **`Homo-720mer-6MX4.ipynb`** - Large 720-component viral capsid
- **`book_chapter_example_system_1.ipynb`** - Comprehensive tutorial example

To run the examples locally:
```bash
git clone https://github.com/JohnsonBiophysicsLab/ionerdss.git
cd ionerdss
pip install -e ".[jupyter]"  # Install with Jupyter support
jupyter lab examples/
```

For additional tutorials, see the [online documentation](https://ionerdss.readthedocs.io/en/latest/ionerdss_tutorials.html).

### Run a quick trial with our server
Go to the [NERDSS server](http://52.15.142.249:5000/).

---

## Documentation
- **User Guide:** [ionerdss user guide](https://ionerdss.readthedocs.io/en/latest/ionerdss_documentation_v1_1.html).
- **API Reference:** [API](https://ionerdss.readthedocs.io/en/latest/ionerdss.html).

You can also build the docs locally using Sphinx:
```bash
# Ensure you are in your activated environment
pip install -e ".[docs]"  # Install documentation dependencies
sphinx-apidoc -o website/source ionerdss
cd website
make html
```
Then open `website/build/html/index.html` in your browser.

---

## Docker Development Environment

For isolated development with Jupyter Lab:
```bash
docker build --no-cache -t ionerdss_dev . 
docker run -it --rm -v $(pwd):/app -p 8888:8888 ionerdss_dev
```
This creates a containerized environment with Jupyter Lab accessible at `http://localhost:8888`.

---

## License
This project is licensed under the GPL‐3.0 License.

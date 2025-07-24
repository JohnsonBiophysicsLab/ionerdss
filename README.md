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

If you want to contribute to the development, you can install the package in editable mode from a local clone of this repository.

**Prerequisites:**
*   [Git](https://git-scm.com/)
*   Python 3.9+
*   Your choice of environment manager: `conda`, `venv`, etc.
*   Optionally, [uv](https://github.com/astral-sh/uv) for faster performance.

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
        uv pip install -r env/requirements-dev.txt
        ```
    *   **Using `pip` and `venv`:**
        ```bash
        python -m venv .venv
        source .venv/bin/activate  # On Windows: .venv\Scripts\activate
        pip install -r env/requirements-dev.txt
        ```

This will install `ionerdss` in editable mode (`-e`) with all dependencies required for development and testing.

---

## Running Tests

To run the unit tests locally, ensure you have installed the development environment (which includes `pytest` and `pytest-cov`).

Then, from the project root folder, run:
```bash
pytest
```

---

## Quick Start

```python
import ionerdss as ion
# Example usage:
# ion.some_function()
```
For extended examples, see the [tutorials](https://ionerdss.readthedocs.io/en/latest/ionerdss_tutorials.html).

### Run a quick trial with our server
Go to the [NERDSS server](http://52.15.142.249:5000/).

---

## Documentation
- **User Guide:** [ionerdss user guide](https://ionerdss.readthedocs.io/en/latest/ionerdss_documentation_v1_1.html).
- **API Reference:** [API](https://ionerdss.readthedocs.io/en/latest/ionerdss.html).

You can also build the docs locally using Sphinx:
```bash
# Ensure you are in your activated environment
pip install ".[dev]"  # Install dev dependencies if you haven't already
pip install sphinx
sphinx-apidoc -o docs/source ionerdss
cd docs
make html
```
Then open `docs/build/html/index.html` in your browser.

---

## Develop using docker container:  
```bash
docker build --no-cache -t ionerdss_dev . 
docker run -it --rm -v $(pwd):/app -p 8888:8888 ionerdss_dev
```

---

## License
This project is licensed under the GPL‐3.0 License.

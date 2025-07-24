# Environment Setup

This folder contains configuration files for setting up the Python environment for developing `ionerdss`.

Our project uses `pyproject.toml` as the single source of truth for dependencies. The files here are for convenience.

---

### Option 1: Using Conda (`environment.yml`)

This is the recommended method if you have Conda installed. It sets up a `conda` environment and then uses `pip` to install `ionerdss` in editable mode with all dependencies.

```bash
# From the project root
conda env create -f env/environment.yml
conda activate ionerdss-dev
```

To update the environment later:
```bash
conda env update -f env/environment.yml --prune
```

---

### Option 2: Using pip or uv (`requirements-dev.txt`)

This method is for those who prefer standard Python virtual environments.

**With `uv` (Recommended for speed):**
```bash
# Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate

# Install with uv
uv pip install -r env/requirements-dev.txt
```

**With `pip`:**
```bash
# Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate

# Install with pip
pip install -r env/requirements-dev.txt
```

---

### Jupyter Kernel

After setting up your environment, you can make it available in Jupyter:
```bash
python -m ipykernel install --user --name ionerdss-dev --display-name "Python (ionerdss-dev)"
```

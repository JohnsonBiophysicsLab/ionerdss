
# Environment Setup for Development

This folder contains configuration files for setting up the Python environment required for **developing this project**. For users, the dependencies are already installed during `pip install` (encoded in `pyproject.toml`).

----

##  Using Conda (`environment.yml`)

This is the recommended method if you have Conda (via [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)) installed.

```bash
# From the project root
conda env create -f environment.yml
conda activate ionerdss_env
````

To update the environment later if the file changes:

```bash
conda env update -f env/environment.yml --prune
```


---

## ProAffinity (`proaffinity/environment.yml`)

ProAffinity-GNN pins numpy 1.x and torch 2.2, which cannot coexist with numpy 2 -- required by `ovito>=3.16`. It therefore gets an environment of its own rather than a place in the dev environment above:

```bash
conda env create -f env/proaffinity/environment.yml
export IONERDSS_PROAFFINITY_PYTHON="$(conda info --base)/envs/ionerdss-proaffinity/bin/python"
```

With that variable set, `predict_affinity=True` works from your normal environment -- ioNERDSS runs the prediction in this one and passes JSON across. There is no need to activate it. See `docs/Proaffinity.md`.

If conda stops with a Terms of Service error for `repo.anaconda.com`, your conda is configured to use Anaconda's `defaults` channels. This environment does not need them; accept the terms conda names, or build it from conda-forge alone:

```bash
conda create -y --override-channels -c conda-forge -n ionerdss-proaffinity python=3.10 pip
conda run -n ionerdss-proaffinity pip install "ioNERDSS[proaffinity]"
```

---

## Questions

If you're unsure which setup to use, start with **Conda** if available.

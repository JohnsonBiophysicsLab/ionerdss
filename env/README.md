
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

## Questions

If you're unsure which setup to use, start with **Conda** if available.

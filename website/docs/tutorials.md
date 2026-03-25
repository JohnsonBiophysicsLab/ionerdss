# Tutorials

The maintained examples now live in the repository-level `tutorials/` folder. The legacy notebooks that previously lived under `website/source/` are no longer the supported examples.

## Primary notebooks

- [`quick_start_6bno.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/quick_start_6bno.ipynb): shortest path from structure input to generated NERDSS files.
- [`ionerdss_tutorial_5l93.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_5l93.ipynb): structure processing and exported visualization outputs for `5l93`.
- [`ionerdss_tutorial_6bno.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_6bno.ipynb): PDB-to-model workflow with ODE outputs and NERDSS-ready files.
- [`ionerdss_tutorial_with_proaffinity_8y7s.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_with_proaffinity_8y7s.ipynb): Tutorial using ProAffinity-GNN to predict binding affinity with larger assembly workflow with generated workspace artifacts.

## Additional examples

- [`additional_examples/4yd9.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/additional_examples/4yd9.ipynb)

## Open locally

```bash
pip install -e ".[jupyter]"
jupyter lab tutorials/
```

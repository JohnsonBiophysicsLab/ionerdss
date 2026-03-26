# Tutorials

The maintained examples now live in the repository-level `tutorials/` folder. The legacy notebooks that previously lived under `website/source/` are no longer the supported examples.

## Primary notebooks

- [`quick_start_6bno.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/quick_start_6bno.ipynb): shortest path from structure input to generated NERDSS files.
- [`ionerdss_tutorial_5l93.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_5l93.ipynb): structure processing and exported visualization outputs for `5l93`.
- [`ionerdss_tutorial_6bno.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_6bno.ipynb): PDB-to-model workflow with ODE outputs and NERDSS-ready files.
- [`ionerdss_tutorial_with_proaffinity_8y7s.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_with_proaffinity_8y7s.ipynb): Tutorial using ProAffinity-GNN to predict binding affinity with larger assembly workflow with generated workspace artifacts.
- [`ionerdss_tutorial_simularium.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_simularium.ipynb): Tutorial on converting NERDSS trajectories to the Simularium 3D viewer format. (Assuming that you have already run one of the above simulations and have the NERDSS output files ready)
- [`ionerdss_tutorial_ovito.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_ovito.ipynb): Tutorial on converting NERDSS trajectories into high-quality GIFs using the OVITO rendering engine.

## Additional examples

- [`additional_examples/4yd9.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/additional_examples/4yd9.ipynb)

## Open locally

```bash
pip install -e ".[jupyter]"
jupyter lab tutorials/
```

# Tutorials

The maintained examples now live in the repository-level `tutorials/` folder. The legacy notebooks that previously lived under `website/source/` are no longer the supported examples.

## Core modeling tutorials

- [`quick_start_6bno.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/quick_start_6bno.ipynb): shortest path from structure input to generated NERDSS files.
- [`ionerdss_tutorial_5l93.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_5l93.ipynb): structure processing and exported visualization outputs for `5l93`.
- [`ionerdss_tutorial_6bno.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_6bno.ipynb): PDB-to-model workflow with ODE outputs and NERDSS-ready files.
- [`ionerdss_dodecahedron_tutorial.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_dodecahedron_tutorial.ipynb): platonic-solid workflow for generating a dodecahedron-based assembly.

## Validation and advanced workflow tutorials

- [`ionerdss_tutorial_pdb_validation.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_pdb_validation.ipynb): structure-validation workflow for comparing an assembled result against the designed coarse-grained target.
- [`ionerdss_tutorial_with_proaffinity_8y7s.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_with_proaffinity_8y7s.ipynb): ProAffinity-enabled workflow for binding-affinity prediction on a larger assembly.
- [`ionerdss_tutorial_simularium.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_simularium.ipynb): convert NERDSS trajectories into the Simularium 3D viewer format.
- [`ionerdss_tutorial_ovito_gif.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/ionerdss_tutorial_ovito_gif.ipynb): render NERDSS trajectories into GIF animations using OVITO.

## Additional examples

- [`additional_examples/4yd9.ipynb`](https://github.com/JohnsonBiophysicsLab/ionerdss/blob/main/tutorials/additional_examples/4yd9.ipynb)

## Supplemental files in `tutorials/`

The tutorials folder also contains supporting outputs and helper files that are referenced by some notebooks:

- `install_ADFR_mac.sh`: helper script for installing the ADFR toolchain used by some ProAffinity workflows

## Open locally

```bash
pip install -e ".[jupyter]"
jupyter lab tutorials/
```

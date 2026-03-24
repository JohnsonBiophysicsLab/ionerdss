# PDB Modeling

The `ionerdss.model.pdb` package is the main structure-to-model pipeline. It handles structure parsing, interface detection, repeated-chain grouping, template generation, visualization, and NERDSS export.

## `PDBModelBuilder`

```python
from ionerdss.model.pdb.main import PDBModelBuilder

builder = PDBModelBuilder("6bno")
system = builder.build_system(workspace_path="6bno_dir")
```

### Responsibilities

- Accept a PDB ID or local structure path.
- Create a managed workspace for logs, structures, outputs, and exported files.
- Parse the structure and detect interfaces.
- Group repeated chains into molecule templates.
- Assemble the final `System`.
- Optionally generate visualization artifacts and NERDSS input files.
- Optionally run the ODE pipeline when enabled in hyperparameters.

## `PDBModelHyperparameters`

This dataclass controls the builder behavior. Common fields include:

- `interface_detect_distance_cutoff`
- `interface_detect_n_residue_cutoff`
- `chain_grouping_matching_mode`
- `generate_visualizations`
- `generate_nerdss_files`
- `nerdss_water_box`
- `nerdss_total_molecule_count`
- `predict_affinity`
- `ode_enabled`
- `ode_time_span`
- `ode_solver_method`
- `count_transition`
- `transition_matrix_size`

Example:

```python
from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

hyperparams = PDBModelHyperparameters(
    interface_detect_distance_cutoff=0.8,
    chain_grouping_matching_mode="sequence",
    ode_enabled=True,
)

builder = PDBModelBuilder("1ABC", hyperparams=hyperparams)
system = builder.build_system(workspace_path="workspace")
```

## Hyperparameter helpers

The helper functions in `ionerdss.model.pdb.api` support exporting, importing, printing, and updating hyperparameters without manually rebuilding the dataclass each time.

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
- Optionally export the one-copy structure-validation workflow when requested through the public API.

## `PDBModelHyperparameters`

`PDBModelHyperparameters` controls the behavior of the full PDB-to-NERDSS pipeline. It is documented on a separate page because the parameter surface is large and the fields affect different stages of the workflow.

See [PDBModelHyperparameters](pdb-hyperparameters.md) for a full parameter-by-parameter reference.

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

The top-level wrapper `ionerdss.build_system_from_pdb(...)` also supports:

- `structure_validation=True`
- `structure_validation_options={...}`

to export the validation-ready NERDSS files alongside the main workspace.

## Hyperparameter helpers

The helper functions in `ionerdss.model.pdb.api` support exporting, importing, printing, and updating hyperparameters without manually rebuilding the dataclass each time.

Useful helpers include:

- `set_hyperparameters(builder, **kwargs)`: create or update the builder's hyperparameter object.
- `export_hyperparameters(builder, filepath)`: save a builder configuration to JSON.
- `import_hyperparameters(builder, filepath)`: load a saved configuration into a builder.
- `print_hyperparameters(builder)`: inspect the current values attached to a builder.

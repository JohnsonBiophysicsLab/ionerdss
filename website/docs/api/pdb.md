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

The top-level wrapper `ionerdss.build_system_from_pdb(...)` and `PDBModelBuilder.build_system(...)` also support:

- `structure_validation=True`
- `structure_validation_options={...}`

to also export the [validation deck](top-level.md#validation-deck) into `structure_validation/` in the workspace, beside the regular NERDSS files in `nerdss_files/`. The deck is exported with the build's hyperparameters, and `PDBModelBuilder` keeps the result on `builder.structure_validation_artifacts`.

## Hyperparameter helpers

The helper functions in `ionerdss.model.pdb.api` support exporting, importing, printing, and updating hyperparameters without manually rebuilding the dataclass each time.

Useful helpers include:

- `set_hyperparameters(builder, **kwargs)`: create or update the builder's hyperparameter object.
- `export_hyperparameters(builder, filepath)`: save a builder configuration to JSON.
- `import_hyperparameters(builder, filepath)`: load a saved configuration into a builder.
- `print_hyperparameters(builder)`: inspect the current values attached to a builder.

## Interface and site names

Every interface type is named after the molecule type it sits on, the molecule type it binds, an index that separates several interface types between the same pair, and, for the two halves of a homodimeric heterotypic pair, an `f`/`b` tag. NERDSS site names may only contain letters and digits, so the parts are concatenated without a separator:

| interaction | interface types | NERDSS sites |
|---|---|---|
| heterodimeric, A binds B | `AB1` on A, `BA1` on B | `ab1`, `ba1` |
| homodimeric heterotypic (head-to-tail) | `AA1f` and `AA1b` on A | `aa1f`, `aa1b` |
| homodimeric homotypic (self-binding) | `AA1` on A | `aa1` |

A single-letter molecule name is written as is. Any other molecule name is written with its length in front, so that names of different lengths can still be told apart once concatenated: the A side of an A–AA interface is `A2AA1` (site `a2aa1`) and the AA side is `2AAA1` (site `2aaa1`), a chain renamed `AA0` by the parser gives `3AA0...`, and names of ten or more characters use a `0` followed by a two-digit length (`012Dodecahedron`). The site label is the interface type name with the molecule names in lower case.

The scheme lives in `ionerdss.model.components.interface_naming`, whose `make_interface_name` and `parse_interface_name` are the only places that spell or split these names.

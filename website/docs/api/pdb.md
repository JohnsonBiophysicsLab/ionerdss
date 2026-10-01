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

### `build_system` arguments

- `workspace_path`: output directory for logs, structures, reports, the system JSON and `nerdss_files/`.
- `hyperparams`: a `PDBModelHyperparameters`. Defaults to the one attached to the builder, else the defaults.
- `molecule_counts`: copies per molecule type for the NERDSS export. When omitted, `nerdss_total_molecule_count` is split across the molecule types by stoichiometry.
- `box_nm`: default `(100.0, 100.0, 100.0)`. The NERDSS export uses it only when `nerdss_water_box` is empty; it is also the validation box when `structure_validation_options` gives none.
- `structure_validation`, `structure_validation_options`: also export the one-copy validation deck, which replaces the regular files in `nerdss_files/` (see [Top-Level API](top-level.md#parameters)); the result is kept on `builder.structure_validation_artifacts`.
- `nerdss_params`: extra `parms.inp` parameters written over the generated ones. Set the time step with the `nerdss_time_step` hyperparameter instead, because the automatic time step, whenever it can be computed, replaces a `timestep` given here.
- `**kwargs`: hyperparameter fields merged over `hyperparams`. Names that are not fields are ignored.

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

to export the validation-ready NERDSS files into the workspace's `nerdss_files/`, where they replace the regular `parms.inp` and `.mol` files.

## Hyperparameter helpers

The helper functions in `ionerdss.model.pdb.api` support exporting, importing, printing, and updating hyperparameters without manually rebuilding the dataclass each time. They are also importable from `ionerdss.model.pdb`, and `PDBModelBuilder` has methods of the same names that take the same arguments without `builder`.

Useful helpers include:

- `set_hyperparameters(builder, **kwargs)`: create or update the builder's hyperparameter object.
- `export_hyperparameters(builder, filepath)`: save a builder configuration to JSON.
- `import_hyperparameters(builder, filepath)`: load a saved configuration into a builder.
- `print_hyperparameters(builder)`: print and return a grouped summary of the values attached to a builder. It leaves out some fields, the NERDSS export settings among them; `builder.hyperparams.to_dict()` holds every field.

## Structure validation workflow

`ionerdss.model.pdb.validation` runs the one-copy validation end to end: export a deck holding one copy of the designed assembly, with irreversible binding and titrated subunits, run NERDSS on it, and compare the assembly that forms with the design.

- `setup_simulation(system, *, workspace_manager=None, box_nm=(100.0, 100.0, 100.0), initial_molecule_count=1, titration_on_rate=1e-5, target_filename="structure_validation_target.json", titration_parms_filename="parms_titrate.inp", parms_overrides=None, designed_coordinates=None, interface_com_proximity_threshold_nm=None)`: write the validation deck and the design target, and return `StructureValidationArtifacts`. `initial_molecule_count` multiplies the designed copy numbers. `prepare(...)` does the same without `initial_molecule_count` and `titration_parms_filename`.
- `run_simulation(artifacts, nerdss_dir, *, sim_index=1, sim_dir_name="validation_output", env=None)`: run the NERDSS executable at `<nerdss_dir>/bin/nerdss` on the titration deck in `nerdss_files/<sim_dir_name>/<sim_index>/`, and return a `StructureValidationSimulationResult` saying whether and when the full assembly formed, the largest assembly seen, and the observed coordinates. `env` adds environment variables, such as `LD_LIBRARY_PATH`, for the executable. A NERDSS process that exits with an error raises `RuntimeError` instead of being reported as a run in which nothing assembled.
- `collect_results(artifacts, simulation_dir)`: read the result of a validation run launched outside ioNERDSS, for example as a cluster job. `simulation_dir` is the directory the run wrote `DATA/` into.
- `align_structure(designed_coordinates, observed_coordinates, *, backend="kabsch", plot=False)`, also available as `compare`: rigid alignment and RMSD, returning a `StructureAlignmentResult`. `backend` may also be `"biopython"`.
- `get_designed_structure(system)`: the designed assembly keyed by molecule instance, with each instance's type, centre of mass and bound interfaces (interface type, binding partner and coordinate), for checking what the design connects.

The site-geometry and box-fit preflight helpers described under [Top-Level API](top-level.md#structure-validation-helpers) (`get_near_com_interface_sites`, `get_degenerate_site_layouts`, `get_interface_com_proximity_message`, `get_designed_assembly_extent`, `get_box_fit_message`) can be imported from this module too.

## Interface and site names

Every interface type is named after the molecule type it sits on, the molecule type it binds, an index that separates several interface types between the same pair, and, for the two halves of a homodimeric heterotypic pair, an `f`/`b` tag. NERDSS site names may only contain letters and digits, so the parts are concatenated without a separator:

| interaction | interface types | NERDSS sites |
|---|---|---|
| heterodimeric, A binds B | `AB1` on A, `BA1` on B | `ab1`, `ba1` |
| homodimeric heterotypic (head-to-tail) | `AA1f` and `AA1b` on A | `aa1f`, `aa1b` |
| homodimeric homotypic (self-binding) | `AA1` on A | `aa1` |

A single-letter molecule name is written as is. Any other molecule name is written with its length in front, so that names of different lengths can still be told apart once concatenated: the A side of an A–AA interface is `A2AA1` (site `a2aa1`) and the AA side is `2AAA1` (site `2aaa1`), a chain renamed `AA0` by the parser gives `3AA0...`, and names of ten or more characters use a `0` followed by a two-digit length (`012Dodecahedron`). The site label is the interface type name with the molecule names in lower case.

The scheme lives in `ionerdss.model.components.interface_naming`, whose `make_interface_name` and `parse_interface_name` are the only places that spell or split these names.

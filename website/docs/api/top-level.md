# Top-Level API

These are the main objects and helpers exposed from `import ionerdss as ion`.

## Public entry points

- `build_system_from_pdb`: convenience wrapper for the structure-to-system pipeline.
- `System`: top-level molecular system container with registries and JSON serialization.
- `Simulation`: helper for editing and running a NERDSS workspace.
- `Analyzer`: entry point for post-processing simulation outputs.
- `ODEPipelineConfig`: configuration object for ODE solves.
- `run_ode_pipeline`: execute the ODE workflow.
- `platonic_solid_generator`: build a platonic-solid-based model.
- `build_system_from_plat`: convenience helper for platonic solid workflows.
- `prepare_structure_validation_for_system`: export the irreversible one-copy validation setup.
- `compare_structure_to_design`: align observed coordinates against the design target and report RMSD.
- `StructureValidationConfig`: configuration dataclass for the validation workflow.
- `StructureValidationArtifacts`: generated files and metadata from validation export.
- `StructureAlignmentResult`: output object returned by rigid alignment.
- `visualize_trajectory_ovito`: render XYZ trajectories to GIFs with OVITO.
- `convert_simularium`: export simulation outputs to Simularium.

## `build_system_from_pdb`

```python
import ionerdss as ion

system = ion.build_system_from_pdb(
    source="4v6x",
    workspace_path="4v6x_dir",
    ode_enabled=True,
    count_transition=True,
)
```

### Parameters

- `source`: PDB identifier such as `"4v6x"` or a local PDB/mmCIF path.
- `workspace_path`: output directory for generated files. Defaults to `<source>_dir`.
- `fetch_format`: optional remote structure format.
- `molecule_counts`: optional explicit counts for NERDSS export.
- `structure_validation`: if `True`, export a one-copy-per-type validation setup during the build.
- `structure_validation_options`: optional settings forwarded to the validation export workflow.
- `**hyperparams_kwargs`: any field accepted by `PDBModelHyperparameters`.

### Returns

A populated `System` object ready for export, simulation setup, or analysis.

## Structure validation helpers

### `prepare_structure_validation_for_system`

Prepare the special validation deck that exports one representative copy per molecule type, forces off-rates to zero, adds titration behavior, and writes the designed target coordinates for later comparison.

Typical return values are packaged in `StructureValidationArtifacts`, including:

- chosen molecule counts
- designed coarse-grained coordinates
- target JSON file path
- generated NERDSS input files
- `preflight_warning_message`: set when the designed assembly graph is disconnected
- `free_interface_warning_message`: set when the design leaves binding capacity unused

### Over-assembly preflight warning

A coarse-grained design can leave interfaces unbound. A molecule *type* declares every
interface it can bind through, but a given copy in the deposited assembly only realizes
the contacts observed in the PDB; anything declared but unrealized is a **free slot**.
Because a fresh copy arrives with all of its own interfaces free, that slot can bind
another subunit and grow the assembly past the deposited stoichiometry — over-assembly.

`get_free_interface_capacity(system)` returns
`{molecule type: {interface type: number of copies leaving it free}}`, and
`get_free_interface_message(system, prefix=...)` formats it. Both live in
`ionerdss.model.pdb.structure_validation`. The message is also raised as a
`RuntimeWarning` during validation export and carried on the artifacts.

Free interfaces are not necessarily a defect. A filament such as actin genuinely
nucleates beyond the deposited asymmetric unit, so the warning says a larger structure
*may* form, not that the model is wrong.

Two caveats worth knowing:

- The check is **static**, over the finished design. Over-assembly can also arise
  kinetically: during assembly, partial complexes transiently expose interfaces, and
  with irreversible binding they can fuse before completing. A saturated design is not
  a guarantee.
- Over-assembly is **unobservable at one copy** of the deposited stoichiometry, since
  the largest possible assembly is then the target itself. Supply more copies to test
  for it.

### `compare_structure_to_design`

Rigidly align observed coarse-grained coordinates onto the design target and compute RMSD.

This function accepts either:

- dictionaries keyed by molecule type
- ordered coordinate arrays

It returns a `StructureAlignmentResult` containing labels, RMSD, the rotation and translation, and the aligned coordinates.

## Rendering and export helpers

### `visualize_trajectory_ovito`

Render an XYZ trajectory using OVITO and optionally save it as a GIF. This requires the `ovito_rendering` optional extra.

### `convert_simularium`

Convert supported simulation outputs into Simularium files for interactive 3D viewing.

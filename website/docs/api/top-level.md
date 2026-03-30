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

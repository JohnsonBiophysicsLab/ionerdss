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
- `**hyperparams_kwargs`: any field accepted by `PDBModelHyperparameters`.

### Returns

A populated `System` object ready for export, simulation setup, or analysis.

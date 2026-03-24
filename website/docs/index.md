# ioNERDSS

`ioNERDSS` provides a Python-first interface for assembling NERDSS models from structures, managing simulation inputs, and analyzing simulation outputs.

## What is documented here

- The supported top-level API exposed by `import ionerdss as ion`
- The current structure-to-model pipeline under `ionerdss.model.pdb`
- Simulation utilities in `ionerdss.nerdss_simulation`
- Analysis tools in `ionerdss.analysis`
- The ODE pipeline for reaction-network kinetics

## Core workflow

1. Build a `System` from a PDB identifier or a local structure file.
2. Export NERDSS files and optional visualizations into a workspace.
3. Run or modify a NERDSS simulation from that workspace.
4. Analyze outputs with the `Analyzer` API.

## Quick start

```python
import ionerdss as ion

system = ion.build_system_from_pdb(
    source="6bno",
    workspace_path="6bno_dir",
    ode_enabled=True,
    count_transition=True,
)

simulation = ion.Simulation("6bno_dir/nerdss_files")
simulation.print_inp_file()

analyzer = ion.Analyzer("path/to/simulation/root")
analyzer.plot.free_energy()
```

## Where to go next

- Start with [Getting Started](getting-started.md) for install and local preview steps.
- Use [Tutorials](tutorials.md) for the maintained notebooks in `tutorials/`.
- Browse [API Reference](api/index.md) for module-level details pulled from the package docstrings.

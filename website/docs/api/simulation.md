# Simulation

`ionerdss.nerdss_simulation.Simulation` works with generated NERDSS workspaces.

## Typical usage

```python
from ionerdss.nerdss_simulation import Simulation

sim = Simulation("6bno_dir/nerdss_files")
sim.print_inp_file()
sim.modify_inp_file({"nItr": 500000})
```

## Supported tasks

- Create or reuse a simulation working directory.
- Modify `.mol` files with `modify_mol_file()`.
- Modify `parms.inp` or related input files with `modify_inp_file()`.
- Add interface states with `add_interface_state()`.
- Print molecule or input file contents for inspection.
- Install or invoke NERDSS-related tooling from a workspace.

This API is aimed at workspace editing and execution rather than model generation.

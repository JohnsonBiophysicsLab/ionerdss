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
- Clone and compile NERDSS with `install_nerdss()`.
- Run NERDSS with `run_new_simulations()`. Each replicate runs in `<sim_dir>/<index>/` (`sim_dir` defaults to `nerdss_output/` in the working directory) with its own `output.log`, `env` adds environment variables (for example `LD_LIBRARY_PATH`) for the executable, and a run that exits with an error raises `RuntimeError` quoting the end of its log.
- Continue finished runs from their `DATA/restart.dat` with `run_restart_simulations()`.

This API is aimed at workspace editing and execution rather than model generation.

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
- Clone and compile NERDSS with `install_nerdss(install_dir)`, which returns the compiled executable.
- Run NERDSS with `run_new_simulations(nerdss_path=...)`. Each replicate runs in `<sim_dir>/<index>/` (`sim_dir` defaults to `nerdss_output/` in the working directory) with its own `output.log`, `env` adds environment variables (for example `LD_LIBRARY_PATH`) for the executable, and a run that exits with an error raises `RuntimeError` quoting the end of its log.
- Continue finished runs from their `DATA/restart.dat` with `run_restart_simulations()`.

This API is aimed at workspace editing and execution rather than model generation.

## Locating NERDSS

Every ioNERDSS function that launches NERDSS (`run_new_simulations()`, `run_restart_simulations()`, `ionerdss.model.pdb.validation.run_simulation()` and the benchmark scripts' `--nerdss_path`) takes a `nerdss_path` and resolves it with `ionerdss.nerdss_simulation.resolve_nerdss_executable()`:

- a path to an executable file is run as-is, whatever it is called;
- a directory is searched for `nerdss`, `bin/nerdss`, `nerdss_mpi` and `bin/nerdss_mpi`, in that order, so a NERDSS checkout and its `bin/` folder both work;
- a bare command name such as `"nerdss"` is looked up on `PATH`;
- `None` (the default) searches `<work_dir>/NERDSS` (where `install_nerdss()` puts it when `install_dir` is the working directory) and then `PATH`.

`~` is expanded. The executable runs in place from each replicate directory; it is not copied there. A `nerdss_mpi` build launched this way runs as a single MPI process. The old `nerdss_dir` keyword, which had to be the checkout containing `bin/nerdss`, is still accepted with a `DeprecationWarning`.

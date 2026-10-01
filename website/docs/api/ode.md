# ODE Pipeline

The ODE pipeline computes reaction-network kinetics for generated systems before or alongside NERDSS simulations.

## `ODEPipelineConfig`

Key configuration fields include:

- `t_span`: integration interval, default `(0.0, 10.0)`.
- `initial_concentrations`: optional species-to-concentration mapping. If omitted, the first species starts at 1.0 and all others at 0.
- `solver_method`: solver such as `"BDF"`.
- `atol`: absolute tolerance.
- `plot`: whether to generate plots.
- `plot_species_indices`: optional subset of species for plotting.
- `plot_sample_points`: number of sampled points for plots.
- `save_csv`: whether to save CSV output.
- `species_labels`: optional custom legend labels.

## Main helpers

These live in `ionerdss.ode.ode_pipeline`; `ODEPipelineConfig` and `run_ode_pipeline` are also exported from `ionerdss`.

- `calculate_ode_solution()`: solve concentrations over time for a reaction system.
- `save_ode_results()`: write CSV and plot outputs.
- `run_ode_pipeline()`: high-level pipeline entry point used by model-building workflows. Returns `(time, concentrations, species_names, saved_files)`.

The reaction system for a built `System` comes from `generate_ode_model_from_system(system)` in `ionerdss.ode.system_ode_generator`, which returns `(complex_names, reaction_system)` and raises `ValueError` when the assembly has more molecules than `max_complex_size` (default 12).

Example:

```python
from ionerdss import ODEPipelineConfig, run_ode_pipeline
from ionerdss.ode.system_ode_generator import generate_ode_model_from_system

_, reaction_system = generate_ode_model_from_system(system)

config = ODEPipelineConfig(
    t_span=(0.0, 100.0),
    solver_method="BDF",
    plot=True,
)

run_ode_pipeline(reaction_system, output_dir="ode_results", config=config)
```

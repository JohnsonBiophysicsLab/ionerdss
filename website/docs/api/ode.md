# ODE Pipeline

The ODE pipeline computes reaction-network kinetics for generated systems before or alongside NERDSS simulations.

## `ODEPipelineConfig`

Key configuration fields include:

- `t_span`: integration interval.
- `initial_concentrations`: optional species-to-concentration mapping.
- `solver_method`: solver such as `"BDF"`.
- `atol`: absolute tolerance.
- `plot`: whether to generate plots.
- `plot_species_indices`: optional subset of species for plotting.
- `plot_sample_points`: number of sampled points for plots.
- `save_csv`: whether to save CSV output.
- `species_labels`: optional custom legend labels.

## Main helpers

- `calculate_ode_solution()`: solve concentrations over time for a reaction system.
- `save_ode_results()`: write CSV and plot outputs.
- `run_ode_pipeline()`: high-level pipeline entry point used by model-building workflows.

Example:

```python
from ionerdss import ODEPipelineConfig, run_ode_pipeline

config = ODEPipelineConfig(
    t_span=(0.0, 100.0),
    solver_method="BDF",
    plot=True,
)

run_ode_pipeline(reaction_system, output_dir="ode_results", config=config)
```

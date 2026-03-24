# Analysis

`ionerdss.analysis` provides the current post-processing API for simulation outputs.

## `Analyzer`

```python
from ionerdss.analysis import Analyzer

analyzer = Analyzer("./simulation_root")
print(len(analyzer.simulations))
analyzer.plot.free_energy()
```

### Constructor

- `Analyzer(root_dir)`: create an analysis controller rooted at a directory containing one or more simulation folders.

On initialization, the analyzer:

- stores `root_dir` as a `Path`
- creates a `DataLoader`
- discovers simulations recursively
- exposes the plot namespace as `analyzer.plot`

### Main attributes

- `root_dir`: filesystem root searched for simulations.
- `loader`: `DataLoader` instance used for discovery.
- `simulations`: list of discovered `Simulation` objects.
- `plot`: `Plotter` namespace bound to this analyzer.

### Analyzer methods

#### `get_simulation(index_or_id)`

Retrieve a simulation either by integer index or by simulation ID.

- accepts: `int | str`
- returns: a single `Simulation`
- raises: `IndexError` for an invalid index, `KeyError` for an unknown ID

Use this method when you want an explicit simulation object before calling lower-level helpers.

#### `load_simulations(simulations=None, time_frame=None)`

Compatibility helper that returns a subset of the discovered simulations.

- if `simulations` is `None`, it returns the full `analyzer.simulations` list
- if `simulations` is a list of indices or IDs, it resolves each one through `get_simulation`
- invalid identifiers are skipped rather than raising
- `time_frame` is currently accepted for compatibility but not applied in this implementation

This is mainly useful when adapting older analysis code to the newer API.

#### `compute_size_distribution(sim)`

Compute a cluster-size distribution from a simulation transition matrix.

- accepts: a `Simulation` object
- internally calls: `sim.get_transition_matrix()`
- returns: a `pandas.DataFrame` with columns:
  - `size`
  - `count`
  - `probability`

If no transition matrix is present, the method logs an error and returns an empty-size distribution through the processing layer.

#### `compute_free_energy(sim, temperature=1.0)`

Compute a free-energy profile from the transition-matrix-derived size distribution.

- accepts: a `Simulation` object and an optional temperature scale
- returns: a `pandas.DataFrame` containing:
  - `size`
  - `count`
  - `probability`
  - `free_energy`

The result is cached in `sim.data.df_free_energy`, so repeated calls on the same simulation avoid recomputing the DataFrame.

## `Plotter`

`analyzer.plot` is a thin namespace that turns processed data into Matplotlib plots. Each method resolves one simulation, computes the needed data, and forwards plotting kwargs to the visualization layer.

### `plot.free_energy(simulation_index=0, ax=None, **kwargs)`

Plots free energy vs. cluster size for a specific simulation.

- resolves the simulation with `get_simulation`
- computes the DataFrame with `compute_free_energy`
- calls `plots.plot_free_energy`
- returns: a Matplotlib `Axes`

Typical kwargs include line styling options such as `color`, `linewidth`, and `linestyle`.

### `plot.size_distribution(simulation_index=0, ax=None, **kwargs)`

Plots the cluster size probability distribution for a specific simulation.

- resolves the simulation with `get_simulation`
- computes the DataFrame with `compute_size_distribution`
- calls `plots.plot_size_distribution`
- returns: a Matplotlib `Axes`

The underlying plot function also supports `log_scale=True` by default.

### `plot.transitions(simulation_index=0, ax=None, **kwargs)`

Plots growth and shrinkage probabilities derived from the simulation transition matrix.

- resolves the simulation
- calls `sim.get_transition_matrix()`
- computes probabilities with `transitions.compute_transition_probabilities`
- renders the result with `plots.plot_growth_probabilities`
- returns: a Matplotlib `Axes`

This is the main convenience function for viewing association vs. dissociation trends as a function of cluster size.

### `plot.heatmap(simulation_index=0, ax=None, **kwargs)`

Plots the raw aggregated transition matrix as a heatmap.

- resolves the simulation
- calls `sim.get_transition_matrix()`
- renders it with `plots.plot_heatmap`
- returns: a Matplotlib `Axes`

Useful kwargs include `log_scale`, `cmap`, and `title`.

## `Simulation` objects exposed through `Analyzer`

Each entry in `analyzer.simulations` is an instance of `ionerdss.analysis.core.simulation.Simulation`. These are often the next layer of API you call after selecting a simulation from the analyzer.

### Main attributes

- `path`: root directory of the simulation
- `id`: simulation identifier, usually derived from the directory name
- `data`: lazily loaded `SimulationData`

### Important `Simulation` methods

#### `load()`

Load transition matrices, lifetimes, copy numbers, and histogram data from the simulation directory.

Expected files are searched under `DATA/`, including:

- `transition_matrix_time.dat`
- `copy_numbers_time.dat`
- `histogram_complexes_time.dat`

#### `get_transition_matrix(time_range=None)`

Aggregate the transition matrices across all recorded time points, or only within a selected `(start, end)` interval.

- returns: a single summed NumPy matrix
- pads smaller matrices if needed before summing
- returns an empty array if no transition data is available

This is the main input to `compute_size_distribution`, `compute_free_energy`, `plot.transitions`, and `plot.heatmap`.

#### `get_lifetimes(cluster_size)`

Return all recorded lifetimes for complexes of one cluster size.

- accepts: integer cluster size
- returns: `list[float]`

#### `get_time_series(complex_name)`

Return histogram-based time series for one or more target complexes.

Accepted complex selectors:

- a string such as `"A: 84. c1: 75. L: 84."`
- a composition dictionary such as `{"A": 84, "c1": 75, "L": 84}`
- a list of either form

Returns:

- `time`: 1D NumPy array
- `counts`: 1D or 2D NumPy array depending on whether one or many complexes were requested

#### `get_largest_size_time_series(include=None, exclude=None, only_count_these=None)`

Compute the time series of the largest matching complex size from histogram data.

This is useful for assembly tracking, for example:

- only complexes containing selected monomers
- excluding contaminants or helper species
- counting only a subset of monomer types in the size definition

#### `get_average_size_time_series(include=None, exclude=None, only_count_these=None)`

Compute the mass-weighted average complex size over time using histogram data.

This is useful when you want a smoother summary statistic than the single largest complex.

## Supporting modules used by `Analyzer`

`Analyzer` is intentionally thin and delegates the numerical work to a few internal modules:

- `analysis.io.loader.DataLoader`: discovers simulation directories containing `DATA/`
- `analysis.processing.transitions`: computes size distributions, free energies, and transition probabilities
- `analysis.visualization.plots`: turns processed tables and matrices into Matplotlib plots

The callable functions currently used by `Analyzer` and `Plotter` are:

- `compute_size_distribution_transition_matrix(transition_matrix)`
- `compute_free_energy(size_dist, temperature=1.0)`
- `compute_transition_probabilities(transition_matrix, symmetric=True)`
- `plot_free_energy(df, ax=None, label=None, **kwargs)`
- `plot_size_distribution(df, ax=None, log_scale=True, label=None, **kwargs)`
- `plot_growth_probabilities(df, ax=None, **kwargs)`
- `plot_heatmap(matrix, ax=None, log_scale=True, cmap="viridis", title="Transition Matrix")`

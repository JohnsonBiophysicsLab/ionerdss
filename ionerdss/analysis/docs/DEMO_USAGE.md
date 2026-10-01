# Demo Usage Cases

This document demonstrates how to use the refactored `ionerdss.analysis` library.

## 1. Basic Setup

```python
from ionerdss.analysis import Analyzer

# Initialize the analyzer with the root directory of your simulations
analyzer = Analyzer(root_dir="./my_simulations")

# The analyzer automatically discovers simulation subdirectories
print(f"Found {len(analyzer.simulations)} simulations.")
```

## 2. Computing Free Energy

```python
# Compute Free Energy profile for the first simulation, from its transition matrix
# (F = -kT ln P(size), temperature=1.0 by default)
# Returns a Pandas DataFrame
df_fe = analyzer.compute_free_energy(analyzer.simulations[0])

print(df_fe.head())
#    size     count  probability  free_energy
# 0     1   5793964     0.157809     1.846368
# 1     2   7769490     0.211616     1.552980
```

## 3. Plotting

### 3.1. Modern API (Recommended)

```python
import matplotlib.pyplot as plt

# Plot Free Energy of one simulation (its index in analyzer.simulations, or its ID)
analyzer.plot.free_energy(simulation_index=0)
plt.show()

# Compare simulations by drawing them on one Axes
fig, ax = plt.subplots()
for i in [0, 1]:
    analyzer.plot.free_energy(simulation_index=i, ax=ax)

# Plot Cluster Size Distribution (log-scale y axis unless log_scale=False)
analyzer.plot.size_distribution(
    simulation_index=0,
    log_scale=False
)

# Growth vs shrinkage probabilities, and the transition-matrix heatmap
analyzer.plot.transitions(simulation_index=0)
analyzer.plot.heatmap(simulation_index=0)
```

### 3.2. Legacy API (Backward Compatibility)

`Analyzer` has no `plot_figure` method. `LegacyPlotInterface` wraps an analyzer and maps
old `plot_figure(figure_type, x=..., y=...)` calls onto `analyzer.plot`. It issues no
DeprecationWarning; for a combination it does not map, it prints a warning and returns `None`.

```python
from ionerdss.analysis import LegacyPlotInterface

legacy = LegacyPlotInterface(analyzer)

# Each call draws the same plot as the analyzer.plot call in its comment
legacy.plot_figure(figure_type="line", x="size", y="free_energy")                     # plot.free_energy(simulation_index=0)
legacy.plot_figure(figure_type="line", x="size", y="growth_probability")              # plot.transitions(simulation_index=0)
legacy.plot_figure(figure_type="hist", x="size", y="complex_count", simulations=[1])  # plot.size_distribution(simulation_index=1)
legacy.plot_figure(figure_type="heatmap")                                             # plot.heatmap(simulation_index=0)
```

`"line"` also takes any other `y` containing `probability`, and `"hist"` any `y` containing
`count`. `simulations`, `x` and `y` only choose the plot, and only the first entry of
`simulations` is drawn. Every other keyword goes on to the `analyzer.plot` method, so the call
accepts what that method accepts; options only the old `plot_figure` had, such as `legend`
or `time_frame`, raise an error.

## 4. Advanced Analysis

### 4.1. Custom Processing

You can access the raw NumPy arrays for custom analysis.

```python
import numpy as np

# Get the transition matrix for Simulation 0, summed over every time point
# (pass time_range=(start, end) to sum only that window;
#  the per-time-point matrices are in analyzer.simulations[0].data.transitions)
# Shape: (N_sizes, N_sizes)
T_matrix = analyzer.simulations[0].get_transition_matrix()

# Calculate custom metric: e.g., Eigenvalues
eigenvals = np.linalg.eigvals(T_matrix)
```

### 4.2. Filtering Data

```python
sim = analyzer.simulations[0]

# Copy numbers are a DataFrame: "Time (s)" plus one column per species NERDSS writes
copies = sim.data.copy_numbers
filtered = copies[copies["A(A1)"] > 5]  # time points where species A(A1) count > 5

# Complexes from histogram_complexes_time.dat, selected by composition
times, counts = sim.get_time_series({"A": 3})  # copies of the A3 complex over time
times, largest = sim.get_largest_size_time_series(include=["A"], exclude=["B"])
times, mean_size = sim.get_average_size_time_series(include=["A"])
```


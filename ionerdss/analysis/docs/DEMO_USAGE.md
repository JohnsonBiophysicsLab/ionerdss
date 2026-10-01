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
# Compute Free Energy profile for the first simulation
# Returns a Pandas DataFrame
df_fe = analyzer.simulations[0].compute_free_energy()

print(df_fe.head())
#    cluster_size  free_energy  probability
# 0             1     0.000000     0.450000
# 1             2     1.204123     0.250000
```

## 3. Plotting

### 3.1. Modern API (Recommended)

```python
import matplotlib.pyplot as plt

# Plot Free Energy
analyzer.plot.free_energy(
    simulations=[0, 1],  # Index of simulations to compare
    time_range=(100.0, 200.0)
)
plt.show()

# Plot Cluster Size Distribution
analyzer.plot.size_distribution(
    simulations="all",
    normalize=True
)
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

# Get the transition matrix for Simulation 0
# Shape: (N_sizes, N_sizes)
T_matrix = analyzer.simulations[0].data.transition_matrix

# Calculate custom metric: e.g., Eigenvalues
eigenvals = np.linalg.eigvals(T_matrix)
```

### 4.2. Filtering Data

```python
# Select data only where 'A' count > 5
filtered_sims = analyzer.filter_simulations(condition="species_A > 5")
```


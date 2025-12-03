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

The old `plot_figure` method still works but issues a DeprecationWarning.

```python
# Old style
analyzer.plot_figure(
    figure_type="line",
    x="size",
    y="free_energy",
    legend=["Sim 1", "Sim 2"]
)
```

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


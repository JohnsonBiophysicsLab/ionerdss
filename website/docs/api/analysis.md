# Analysis

`ionerdss.analysis` provides the current post-processing API for simulation outputs.

## `Analyzer`

```python
from ionerdss.analysis import Analyzer

analyzer = Analyzer("./simulation_root")
print(len(analyzer.simulations))
analyzer.plot.free_energy()
```

### Responsibilities

- Discover simulations below a root directory.
- Load transition matrix and other supported output files.
- Retrieve simulations by index or ID.
- Compute size distributions from transition matrices.
- Compute free-energy profiles and cache results.
- Expose a plotting namespace through `analyzer.plot`.

## Plotting helpers

The `Plotter` namespace currently includes:

- `free_energy()`
- `size_distribution()`
- `transitions()`
- `heatmap()`

These methods operate on a selected simulation and return Matplotlib objects from the visualization layer.

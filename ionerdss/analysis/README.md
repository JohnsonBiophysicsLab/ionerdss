# ionerdss/analysis/README.md

# Ionerdss Analysis Module

A clean, modular, and extendable framework for analyzing and visualizing biophysical simulation data.

## Overview

This module provides a set of tools to process simulation data, compute biophysical properties (e.g., copy numbers, complex sizes, free energy, transition probabilities), and generate publication-quality plots.

The architecture follows a clear separation of concerns:
1.  **Data Layer (`analysis.data`)**: Handles reading, parsing, caching, and processing raw simulation data.
2.  **Analysis Layer**: Built into the data processors, extracting meaningful statistics (e.g., `TransitionProcessor`, `CopyNumberProcessor`).
3.  **Visualization Layer (`analysis.plotting`)**: Generates plots from processed data.
4.  **Core Interface (`analysis.core.Analyzer`)**: The main entry point for users.

## Key Components

### 1. Analyzer (Main Interface)

The `Analyzer` class orchestrates the workflow. It discovers simulation directories, manages data loading, and configures plotting.

**Usage:**

```python
from ionerdss.analysis import Analyzer

# Initialize analyzer in the current directory
analyzer = Analyzer(save_dir=".")

# 1. Quick Plotting (Simplest way)
# Automatically loads data and plots
analyzer.quick_plot(
    "line_speciescopy_vs_time",
    legend=["A", "B"],
    figure_size=(12, 8)
)

# 2. Advanced Workflow (More control)
# Configure data selection
data = analyzer.load_simulations(
    simulations=[0, 1, 2],          # Select specific simulation runs
    time_frame=(0, 100)             # Select time range
)

# Configure plotting style
plotter = analyzer.set_plot(
    style="seaborn-whitegrid",
    font_size=14
)

# Generate specific plots
plotter.line_speciescopy_vs_time(data, legend=["A", "B"])
plotter.histogram_complex_size(data, legend=["A", "B"], bins=20)
```

### 2. Data Handling (`analysis.data`)

The `Data` class manages access to simulation data. It uses specialized processors:
-   `CopyNumberProcessor`: For time-series species counts.
-   `HistogramProcessor`: For complex size distributions.
-   `TransitionProcessor`: For transition matrices and lifetimes.

Data is automatically cached to improve performance on repeated analysis.

### 3. Plotting (`analysis.plotting`)

The plotting module is organized by chart type:
-   `line_plots`: Time-series and continuous variable plots.
-   `histogram_plots`: Distribution plots.
-   `heatmap_plots`: 2D density/frequency visualizations.
-   `three_d_plots`: 3D visualizations.
-   `probability_plots`: Biophysical probability analyses (free energy, association/dissociation).

## Extending the Module

### Adding a New Plot Type

1.  Create a new function in the appropriate module (e.g., `analysis/plotting/custom_plots.py`).
2.  Ensure it accepts a `Data` object as the first argument.
3.  Register it in `analysis/plotting/core.py` within the `PlotConfigure` class.

### Adding a New Data Processor

1.  Create a new processor class in `analysis/data/processors/`.
2.  Implement `read` and `process` methods.
3.  Integrate it into `analysis/data/core.py`.

## Directory Structure

```
analysis/
├── __init__.py            # Exports Analyzer
├── core.py                # Main Analyzer class
├── data/                  # Data layer
│   ├── core.py            # Data management class
│   └── processors/        # Specialized data processors
├── plotting/              # Visualization layer
│   ├── core.py            # Plot configuration
│   ├── line_plots.py
│   ├── histogram_plots.py
│   └── ...
└── legacy/                # Backward compatibility
```

## Dependencies

-   `numpy`
-   `pandas`
-   `matplotlib`
-   `seaborn`



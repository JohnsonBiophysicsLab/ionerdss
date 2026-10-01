# Design Document: `ionerdss.analysis` Refactoring

## 1. Architecture Overview

The refactored `ionerdss.analysis` module adopts a **Layered Architecture** to separate concerns, ensuring maintainability, testability, and performance.

### Layers
1.  **API Layer (`api.py`)**: The high-level entry point for users. It orchestrates the lower layers.
2.  **Processing Layer (`processing/`)**: Contains pure functions for scientific computation. It relies on NumPy/Pandas and is agnostic of file I/O.
3.  **Core Layer (`core/`)**: Defines the fundamental data structures (`Simulation`, `SimulationData`) and types used across the system.
4.  **I/O Layer (`io/`)**: Handles the messy details of parsing legacy file formats and interacting with the file system.
5.  **Visualization Layer (`visualization/`)**: Pure plotting logic using Matplotlib/Seaborn, accepting standard data structures (DataFrames). The exceptions are the movie renderers, which read a run's frames from disk: `trajectory_movie.py` (`render_trajectory_movie`, reading `parms.inp`, the `.mol` files and the PDB or XYZ frames) and the deprecated `pymol_movie.py`.
6.  **Legacy Layer (`legacy/`)**: `LegacyPlotInterface`, which maps the old `plot_figure` calls onto the API layer.

---

## 2. Class Design

### 2.1. API Layer
*   **`Analyzer`**: The main facade class.
    *   **Attributes**:
        *   `simulations`: List of `Simulation` objects, discovered under `root_dir` by `io.loader.DataLoader` (every directory up to three levels down that holds a `DATA/` folder).
        *   `plot`: A `Plotter` bound to the analyzer.
    *   **Methods**:
        *   `__init__(root_dir: str | Path)`
        *   `get_simulation(index_or_id)`, `load_simulations(simulations=None, time_frame=None)`: Look simulations up by index or ID (`time_frame` is accepted but ignored).
        *   `compute_size_distribution(sim)`, `compute_free_energy(sim, temperature=1.0)`: Delegate to `processing`.
*   **`Plotter`** (`analyzer.plot`): `free_energy(...)`, `size_distribution(...)`, `transitions(...)` and `heatmap(...)`, each taking `simulation_index` and an optional `ax`; they delegate to `processing` then `visualization`.

### 2.2. Core Layer
*   **`Simulation`**: Represents a single simulation run.
    *   **Attributes**:
        *   `path`: Path to the simulation directory.
        *   `id`: Identifier, the directory name by default.
        *   `data`: The run's `SimulationData` (transition matrices, lifetimes, copy numbers, complex histogram), read from `DATA/` on first access.
    *   **Methods**:
        *   `get_transition_matrix(time_range=None)`: Sums the per-time-point transition matrices, optionally within a time window.
        *   `get_lifetimes(cluster_size)`: All recorded lifetimes for one cluster size.
        *   `get_time_series(...)`, `get_largest_size_time_series(...)`, `get_average_size_time_series(...)`: Time series from the complex histogram.

### 2.3. Processing Layer
*   **`processing/transitions.py`** (module-level functions):
    *   `compute_size_distribution_transition_matrix(transition_matrix: np.ndarray) -> pd.DataFrame`
    *   `compute_free_energy(size_dist: pd.DataFrame, temperature: float = 1.0) -> pd.DataFrame`
    *   `compute_transition_probabilities(transition_matrix: np.ndarray, symmetric: bool = True) -> pd.DataFrame`
    *   **Algorithm**:
        *   Size distribution and free energy use **Vectorization** (NumPy): axis-wise sums over the transition matrix $T_{ij}$.
        *   Growth/shrinkage probabilities loop over sizes, because the symmetric-collision correction depends on the index.

---

## 3. Algorithms & Optimization

### 3.1. Transition Matrix Processing
*   **Old Approach**: Iterating through lines of text files, parsing integers manually, and summing in loops.
*   **New Approach**:
    1.  Parse the file once into a list of `{"time", "matrix"}` records, one dense 2D NumPy array per time point; `Simulation.get_transition_matrix` sums them, zero-padding any that differ in shape.
    2.  **Free Energy**: $G(n) = -k_B T \ln(P(n))$. Calculated via `P(n) = matrix.sum(axis=1) / total`.
    3.  **Growth/Shrinkage**: For each size, transitions to larger sizes are weighed against transitions to smaller ones; with `symmetric=True`, a transition that doubles or halves the size (e.g. monomer + monomer) counts as half.

### 3.2. Data Loading
*   **Lazy Loading**: Data files (which can be GBs) are only read when accessed, not on initialization.
*   **Caching**: `Analyzer.compute_free_energy` keeps its result in memory on `SimulationData.df_free_energy` and returns it on later calls (the cache is not keyed by `temperature`). Nothing is cached on disk.

### 3.3. Parsing
*   **Regex**: Compiled Regex patterns are used for robustly identifying data blocks in the legacy text files, handling edge cases like inconsistent spacing or typos ("transion matrix").


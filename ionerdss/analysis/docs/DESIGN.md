# Design Document: `ionerdss.analysis` Refactoring

## 1. Architecture Overview

The refactored `ionerdss.analysis` module adopts a **Layered Architecture** to separate concerns, ensuring maintainability, testability, and performance.

### Layers
1.  **API Layer (`api.py`)**: The high-level entry point for users. It orchestrates the lower layers.
2.  **Processing Layer (`processing/`)**: Contains pure functions and classes for scientific computation. It relies on NumPy/Pandas and is agnostic of file I/O.
3.  **Core Layer (`core/`)**: Defines the fundamental data structures (`Simulation`, `Trajectory`) and types used across the system.
4.  **I/O Layer (`io/`)**: Handles the messy details of parsing legacy file formats and interacting with the file system.
5.  **Visualization Layer (`visualization/`)**: Pure plotting logic using Matplotlib/Seaborn, accepting standard data structures (DataFrames).

---

## 2. Class Design

### 2.1. API Layer
*   **`Analyzer`**: The main facade class.
    *   **Attributes**:
        *   `simulations`: List of `Simulation` objects.
        *   `config`: Configuration settings.
    *   **Methods**:
        *   `__init__(root_dir: str)`
        *   `load_data(...)`: Delegates to `io`.
        *   `plot_free_energy(...)`: Delegates to `processing` then `visualization`.

### 2.2. Core Layer
*   **`Simulation`**: Represents a single simulation run.
    *   **Attributes**:
        *   `path`: Path to the simulation directory.
        *   `metadata`: Dictionary of parameters (e.g., volume, box size).
    *   **Methods**:
        *   `get_transition_matrix()`: Lazy-loads matrix data.
        *   `get_trajectory()`: Lazy-loads XYZ coordinates.

### 2.3. Processing Layer
*   **`TransitionAnalyzer`**:
    *   **Methods**:
        *   `compute_free_energy(matrix: np.ndarray) -> pd.DataFrame`
        *   `compute_flux(matrix: np.ndarray) -> pd.DataFrame`
    *   **Algorithm**:
        *   Uses **Vectorization** (NumPy) instead of Python loops.
        *   Transition Matrix $T_{ij}$ operations are matrix multiplications or axis-wise sums.

---

## 3. Algorithms & Optimization

### 3.1. Transition Matrix Processing
*   **Old Approach**: Iterating through lines of text files, parsing integers manually, and summing in loops.
*   **New Approach**:
    1.  Parse the file once into a 3D NumPy array `(Time, Size_From, Size_To)` or a list of sparse matrices.
    2.  **Free Energy**: $G(n) = -k_B T \ln(P(n))$. Calculated via `P(n) = matrix.sum(axis=0) / total`.
    3.  **Flux**: Vectorized subtraction of upper/lower triangles of the transition matrix.

### 3.2. Data Loading
*   **Lazy Loading**: Data files (which can be GBs) are only read when accessed, not on initialization.
*   **Caching**: Processed DataFrames are cached in memory (LRU Cache) or on disk (Parquet format) to speed up repeated plotting.

### 3.3. Parsing
*   **Regex**: Compiled Regex patterns are used for robustly identifying data blocks in the legacy text files, handling edge cases like inconsistent spacing or typos ("transion matrix").


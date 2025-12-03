# Developer Guide

## 1. Environment Setup

1.  **Install Dependencies**:
    ```bash
    pip install numpy pandas matplotlib seaborn scipy pytest
    ```
2.  **Install Pre-commit Hooks** (Optional but recommended):
    ```bash
    pre-commit install
    ```

## 2. Project Structure

*   `core/`: Defines data types (`types.py`) and the `Simulation` class. **Do not import high-level modules (like `api`) here** to avoid circular imports.
*   `io/`: Handles file reading and parsing. **No plotting code here.**
*   `processing/`: Pure math and scientific calculations. **No file I/O or plotting here.**
*   `visualization/`: Plotting code. **Input is always DataFrames/Arrays, never file paths.**

## 3. Extending Data Support (How to Add a New Parser)

If you have a new data file (e.g., `energy_time.dat`) you need to parse:

### Step 1: Write the Parser
Create a parsing function in `io/parser.py`.
*   **Input**: `pathlib.Path` to the file.
*   **Output**: `pd.DataFrame` or a structured Dict/List.
*   **Error Handling**: Return an empty structure if the file is missing or corrupt, and log a warning.

```python
# io/parser.py
def parse_energy_data(file_path: Path) -> pd.DataFrame:
    if not file_path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(file_path, delimiter="\t")
    except Exception as e:
        logger.error(f"Failed to parse {file_path}: {e}")
        return pd.DataFrame()
```

### Step 2: Update Core Types
Add a field for this data in `core/types.py` inside the `SimulationData` dataclass.

```python
# core/types.py
@dataclass
class SimulationData:
    # ... existing fields ...
    energy_data: Optional[pd.DataFrame] = None  # <--- Add this
```

### Step 3: Integrate into Loading Logic
Update `core/simulation.py` to identify the file and call your parser.

```python
# core/simulation.py
class Simulation:
    def load(self):
        # ... existing file definitions ...
        energy_file = self.path / "DATA" / "energy_time.dat"

        # ... Call your parser ...
        energy_df = parser.parse_energy_data(energy_file)

        # ... Store in SimulationData ...
        self._data = SimulationData(
            # ... existing fields ...
            energy_data=energy_df
        )
```

## 4. Adding a New Plot Type

1.  **Implement logic in `processing/`** (if calculation is needed):
    ```python
    # processing/new_metric.py
    def calculate_my_metric(data: np.ndarray) -> pd.DataFrame:
        # Perform calculation using vectorization
        ...
    ```
2.  **Add Plot Function in `visualization/plots.py`**:
    ```python
    def plot_my_metric(df: pd.DataFrame, ax=None, **kwargs):
        # Pure plotting logic
        ...
    ```
3.  **Expose in `api.py`**:
    Add a wrapper method in `Analyzer.plot` or `Analyzer`.

## 5. Testing

Run tests using `pytest`:

```bash
pytest tests/analysis/
```

### Writing Tests
*   **Unit Tests**: Test individual functions in `processing/` using mock NumPy arrays.
*   **Integration Tests**: Test `io/` using a small sample data file (stored in `tests/data/`).

## 6. Code Style

*   **Type Hints**: Required for all function arguments and return values.
*   **Docstrings**: Google Style.
    ```python
    def my_func(x: int) -> int:
        """
        Squares a number.

        Args:
            x (int): The input number.

        Returns:
            int: The squared number.
        """
        return x * x
    ```
*   **Formatter**: Black.
*   **Linter**: Flake8 / Ruff.


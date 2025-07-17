# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Testing
```bash
# Run all tests
pytest

# Run tests with coverage
pytest --cov=ionerdss

# Run specific test file
pytest tests/test_specific_file.py

# Run tests for a specific module
pytest tests/unit/math/

# Run tests with verbose output
pytest -v
```

### Installation for Development
```bash
# Install package in editable mode with dependencies
pip install -e .

# Install with optional dependencies
pip install -e ".[ovito_rendering,tests]"

# Install from requirements file
pip install -r env/requirements.txt
```

### Documentation
```bash
# Build documentation locally
cd website
make html
# Open website/build/html/index.html in browser

# Generate API documentation
sphinx-apidoc -o docs/source ionerdss
```

## Code Architecture

### High-Level Structure
ionerdss is a Python package for preparing NERDSS inputs and analyzing simulation outputs. It consists of several main modules:

- **model/**: Model building tools for setting up NERDSS simulations
  - `pdb_model.py`: Convert PDB structures to NERDSS format
  - `design_model.py`: Design custom molecular models
  - `complex.py`: Handle complex molecular systems
  - `platonic_solids/`: Generate platonic solid geometries for viral capsids
  - `graph_based/`: Graph-based analysis of molecular complexes

- **analysis/**: Data analysis and visualization tools
  - `core.py`: Main Analysis class with modular API
  - `data/`: Data processing pipeline with specialized processors
  - `plotting/`: Plotting functions organized by plot type
  - `legacy/`: Backward compatibility interface

- **nerdss_simulation/**: Simulation execution and management
  - `simulation.py`: Main Simulation class for running NERDSS

- **math/**: Mathematical utilities
  - `angles.py`, `coords.py`, `rotations.py`, `inertia_tensors.py`: Geometric calculations

- **nerdss_guis/**: PyQt6-based GUI applications
  - `gui.py`: Main GUI launcher
  - `nerdss.py`: PDB parsing and visualization GUI

- **ode_solver/**: ODE solving capabilities
  - `reaction_ode_solver.py`: Solve reaction kinetics ODEs
  - `reaction_string_parser.py`: Parse reaction strings

- **gillespie_simulation/**: Gillespie algorithm implementation
  - `simple_gillespie.py`: Stochastic simulation algorithm

### Key Design Patterns
- **Lazy Loading**: Main module uses LazyLoader to avoid importing heavy dependencies until needed
- **Modular Analysis API**: New analysis system with separate data processing and plotting components
- **Legacy Compatibility**: Maintains backward compatibility with older API versions
- **Optional Dependencies**: Core functionality works without heavy dependencies like ovito

### Important Files
- `pyproject.toml`: Project configuration and dependencies
- `env/requirements.txt`: Core dependencies list
- `examples/`: Jupyter notebooks demonstrating usage
- `data/`: Test data and example structures
- `tests/`: Unit tests organized by module

### Entry Points
- Import with `import ionerdss as ion`
- Main classes: `Model`, `PDBModel`, `Simulation`, `Analysis` (via `Analyzer`)
- GUI access: `ion.gui()` and `ion.pdb_gui()`
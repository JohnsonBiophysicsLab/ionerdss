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

## PDB Model Refactoring Plan

### Overview
The `pdb_model.py` file was originally a 2,470-line monolithic "God Object" that violated software design principles. We are systematically refactoring it into smaller, focused, testable components.

### Refactoring Phases

#### ✅ Phase 1: Extract Pure Functions (COMPLETED)
**Goal**: Remove utility functions and data classes from the main PDBModel class

**Completed Work:**
- ✅ Created `geometry_utils.py` - geometric calculations and signature comparisons
- ✅ Created `energy_tables.py` - residue interaction energy data  
- ✅ Created `data_classes.py` - core data structures and transformation functions
- ✅ Created `angle_utils.py` - binding angle calculations
- ✅ Updated `pdb_model.py` imports to use extracted utilities
- ✅ Removed 717 lines of duplicate code from `pdb_model.py` (29% reduction)
- ✅ All tests pass - functionality preserved

**Results**: `pdb_model.py` reduced from 2,470 lines → 1,753 lines

#### 🔄 Phase 2: Extract Service Classes (NEXT)
**Goal**: Break down the remaining large methods into focused service classes

**Planned Service Classes:**
1. **`StructureProvider`** - handles PDB/CIF file download and parsing
   - `download_pdb()` → `StructureProvider.download()`
   - `pdb_parser()` → `StructureProvider.parse()`
   - Encapsulates file handling and Biopython structure creation

2. **`CoarseGrainer`** - handles interface detection and coarse-graining
   - `coarse_grain()` method (~130 lines) → `CoarseGrainer.detect_interfaces()`
   - Processes chains, calculates COM, detects binding interfaces
   - Uses KDTree for efficient distance calculations

3. **`HomologyDetector`** - handles chain grouping and homology detection
   - `identify_homologous_chains()` → `HomologyDetector.group_chains()`
   - `_parse_pdb_header()`, `_parse_cif_header()` → `HomologyDetector.parse_headers()`
   - `_find_homologous_chains_by_alignment()` → `HomologyDetector.align_sequences()`

4. **`GeometryRegularizer`** - handles alignment and template generation
   - Core of `regularize_homologous_chains()` (~300 lines) → `GeometryRegularizer.regularize()`
   - Signature calculation, template creation, alignment logic
   - Most complex refactoring due to deep nesting and state management

5. **`ReactionBuilder`** - handles reaction generation
   - `_build_reactions()` (~100 lines) → `ReactionBuilder.build_reactions()`
   - Angle calculations, rate calculations, template creation

**Expected API after Phase 2:**
```python
class PDBModel(Model):
    def generate_model(self):
        # Clean orchestration instead of monolithic methods
        structure = StructureProvider(self.pdb_id).get_structure()
        interfaces = CoarseGrainer(structure).detect_interfaces()
        groups = HomologyDetector(interfaces).group_chains()
        templates = GeometryRegularizer(groups).create_templates()
        reactions = ReactionBuilder(templates).build_reactions()
        return ModelWriter().save(templates, reactions)
```

#### 🔄 Phase 3: Finalize and Test (FUTURE)
**Goal**: Clean up remaining code and ensure robustness

**Planned Work:**
- Extract visualization methods (`plot_*`, `save_*`) to separate utility
- Extract remaining helper methods (`_is_existing_*`, `_validate_*`)
- Add comprehensive unit tests for each service class
- Performance testing and optimization
- Documentation updates

### Current State
- **File size**: 1,753 lines (down from 2,470)
- **Functionality**: 100% preserved, all tests pass
- **Maintainability**: Significantly improved with extracted utilities
- **Next focus**: `regularize_homologous_chains()` method (300 lines, most complex)

### Key Principles
- **Preserve all functionality** - no breaking changes to existing API
- **Test-driven refactoring** - run tests after each change
- **Incremental approach** - small, safe steps with frequent commits
- **Single responsibility** - each class should have one clear purpose
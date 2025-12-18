# PDB Model Hyperparameters API Reference

Quick reference guide for configuring PDB model hyperparameters in ionerdss.

## Table of Contents

- [Quick Start](#quick-start)
- [Model Methods](#model-methods)
- [Common Configurations](#common-configurations)
- [Complete Examples](#complete-examples)
- [Parameter Reference](#parameter-reference)

## Quick Start

```python
from ionerdss.model import pdb

# Create model
model = pdb.PDBModelBuilder("1ABC")

# Set hyperparameters
model.set_hyperparameters(
    interface_detect_distance_cutoff=0.8,
    ode_enabled=True
)

# Build system (hyperparameters automatically used)
system = model.build_system(workspace_path="./workspace")
```

## Model Methods

### `model.set_hyperparameters(**kwargs)`

Set or update hyperparameters on the model. Creates new if none exist, or updates existing ones.

```python
model = pdb.PDBModelBuilder("1ABC")

# Set with defaults
model.set_hyperparameters()

# Set custom values
model.set_hyperparameters(
    interface_detect_distance_cutoff=0.8,
    interface_detect_n_residue_cutoff=5,
    chain_grouping_matching_mode="sequence"
)

# Update existing (preserves other values)
model.set_hyperparameters(ode_enabled=True)
```

### `model.export_hyperparameters(filepath)`

Export model's hyperparameters to a JSON file.

```python
model.set_hyperparameters(
    interface_detect_distance_cutoff=0.8,
    ode_enabled=True
)
model.export_hyperparameters("config.json")
```

### `model.import_hyperparameters(filepath)`

Load hyperparameters from a JSON file into the model.

```python
model = pdb.PDBModelBuilder("1ABC")
model.import_hyperparameters("config.json")
# Hyperparameters now loaded and ready to use
```

### `model.print_hyperparameters()`

Display model's current hyperparameters in a human-readable format.

```python
model.set_hyperparameters(interface_detect_distance_cutoff=0.8)
model.print_hyperparameters()
```

## Common Configurations

### High-Resolution Structures (<2.5 Å)

Tighter thresholds for well-resolved structures:

```python
model = pdb.PDBModelBuilder("1ABC")
model.set_hyperparameters(
    interface_detect_distance_cutoff=0.5,
    interface_detect_n_residue_cutoff=5,
    chain_grouping_rmsd_threshold=1.0,
    chain_grouping_seq_threshold=0.9
)
```

### Low-Resolution Structures (>3.5 Å)

More permissive thresholds for poorly resolved structures:

```python
model = pdb.PDBModelBuilder("2XYZ")
model.set_hyperparameters(
    interface_detect_distance_cutoff=1.2,
    interface_detect_n_residue_cutoff=3,
    chain_grouping_rmsd_threshold=5.0,
    chain_grouping_seq_threshold=0.3
)
```

### Ring/Cyclic Structures

Enable ring regularization:

```python
model = pdb.PDBModelBuilder("3ABC")
model.set_hyperparameters(
    ring_regularization_mode="separate",
    ring_geometry="sphere",
    min_ring_size=4
)
```

### ODE Pipeline Enabled

Enable kinetic modeling:

```python
model = pdb.PDBModelBuilder("4XYZ")
model.set_hyperparameters(
    ode_enabled=True,
    ode_time_span=(0.0, 100.0),
    ode_solver_method="BDF",
    ode_plot=True
)
```

### ProAffinity Binding Energy Prediction

Enable GNN-based affinity prediction:

```python
model = pdb.PDBModelBuilder("5ABC")
model.set_hyperparameters(
    predict_affinity=True,
    adfr_path="/path/to/prepare_receptor"  # Optional
)
```

### Steric Clash Detection

Enable automatic clash detection:

```python
model = pdb.PDBModelBuilder("6XYZ")
model.set_hyperparameters(
    steric_clash_mode="auto"
)
```

## Complete Examples

### Basic Workflow

```python
from ionerdss.model import pdb

# Create model
model = pdb.PDBModelBuilder("1ABC")

# Configure
model.set_hyperparameters(
    interface_detect_distance_cutoff=0.5,
    interface_detect_n_residue_cutoff=5,
    steric_clash_mode="auto",
    generate_visualizations=True
)

# Validate
errors = model.hyperparams.validate()
if errors:
    raise ValueError(f"Invalid configuration: {errors}")

# Save for reproducibility
model.export_hyperparameters("my_config.json")

# Build system (hyperparameters automatically used)
system = model.build_system(workspace_path="./1ABC_workspace")
```

### Multiple Features Enabled

```python
from ionerdss.model import pdb

# Create model and configure all features
model = pdb.PDBModelBuilder("8Y7S")
model.set_hyperparameters(
    # Core detection
    interface_detect_distance_cutoff=0.7,
    interface_detect_n_residue_cutoff=4,
    
    # Chain grouping
    chain_grouping_matching_mode="default",
    chain_grouping_rmsd_threshold=2.0,
    
    # Features
    steric_clash_mode="auto",
    ring_regularization_mode="uniform",
    homotypic_detection="auto",
    
    # ODE pipeline
    ode_enabled=True,
    ode_time_span=(0.0, 50.0),
    ode_solver_method="BDF",
    
    # Affinity prediction
    predict_affinity=True,
    
    # Transition matrix
    count_transition=True,
    transition_matrix_size=500
)

# Build system
system = model.build_system(
    workspace_path="./8Y7S_workspace",
    molecule_counts={"typeA": 20, "typeB": 20},
    box_nm=(150.0, 150.0, 150.0)
)
```

### Configuration Reuse

```python
from ionerdss.model import pdb

# Save configuration from first model
model1 = pdb.PDBModelBuilder("1ABC")
model1.set_hyperparameters(interface_detect_distance_cutoff=0.7)
model1.export_hyperparameters("config.json")

# Reuse in second model
model2 = pdb.PDBModelBuilder("2XYZ")
model2.import_hyperparameters("config.json")

# Adjust specific values
model2.set_hyperparameters(
    interface_detect_distance_cutoff=1.0,  # Override for lower resolution
    chain_grouping_rmsd_threshold=4.0
)

# Build with updated config
system = model2.build_system(workspace_path="./2XYZ_workspace")
```

### Multiple Independent Models

```python
from ionerdss.model import pdb

# Each model has its own configuration
model1 = pdb.PDBModelBuilder("1ABC")
model1.set_hyperparameters(
    interface_detect_distance_cutoff=0.5,  # High resolution
    ode_enabled=True
)

model2 = pdb.PDBModelBuilder("2XYZ")
model2.set_hyperparameters(
    interface_detect_distance_cutoff=1.2,  # Low resolution
    predict_affinity=True
)

# Build independently with different settings
system1 = model1.build_system(workspace_path="./workspace1")
system2 = model2.build_system(workspace_path="./workspace2")
```

## Parameter Reference

### Core Detection Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `interface_detect_distance_cutoff` | float | 0.6 nm | Contact search radius per atom pair |
| `interface_detect_n_residue_cutoff` | int | 3 | Minimum contacting residues per chain |

### Chain Grouping Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `chain_grouping_rmsd_threshold` | float | 2.0 Å | RMSD threshold for structure superposition |
| `chain_grouping_seq_threshold` | float | 0.5 | Sequence identity threshold (50%) |
| `chain_grouping_matching_mode` | str | "default" | Mode: "default", "sequence", "structure" |

### Steric Clash Detection

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `steric_clash_mode` | str | "off" | Mode: "off", "auto", "custom" |

### Template Building Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `signature_precision` | int | 6 | Decimal places for geometric signatures |
| `homodimer_distance_threshold` | float | 0.5 nm | Distance threshold for homodimer detection |
| `homodimer_angle_threshold` | float | 0.5 rad | Angle threshold for homodimer detection |

### Homotypic Detection Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `homotypic_detection` | str | "auto" | Mode: "auto", "signature", "off" |
| `homotypic_detection_residue_similarity_threshold` | float | 0.7 | Residue similarity threshold (70%) |
| `homotypic_detection_interface_radius` | float | 8.0 Å | Interface detection radius |

### Ring Regularization Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ring_regularization_mode` | str | "uniform" | Mode: "off", "separate", "uniform" |
| `ring_geometry` | str | "cylinder" | Geometry: "cylinder", "sphere" |
| `min_ring_size` | int | 3 | Minimum subunits to form a ring |

### Output Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `generate_visualizations` | bool | True | Generate visualization outputs |
| `generate_nerdss_files` | bool | True | Generate NERDSS simulation files |

### ProAffinity Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `predict_affinity` | bool | False | Enable binding energy prediction |
| `adfr_path` | str | None | Path to ADFR prepare_receptor tool |

### ODE Pipeline Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ode_enabled` | bool | False | Enable ODE kinetic modeling |
| `ode_time_span` | tuple | (0.0, 10.0) | Time span (start, end) in seconds |
| `ode_solver_method` | str | "BDF" | Solver method for stiff systems |
| `ode_atol` | float | 1e-4 | Absolute tolerance |
| `ode_plot` | bool | True | Generate plots |
| `ode_save_csv` | bool | True | Save results to CSV |
| `ode_initial_concentrations` | dict | None | Custom initial concentrations |

### Transition Matrix Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `count_transition` | bool | False | Enable transition matrix tracking |
| `transition_matrix_size` | int | 500 | Size of transition matrix |
| `transition_write` | int | None | Write interval (defaults to nItr/10) |

## Tips and Best Practices

1. **Set once, use automatically**: Hyperparameters set on the model are automatically used in `build_system()`
2. **Update incrementally**: Call `set_hyperparameters()` multiple times to update values - previous values are preserved
3. **Save configurations**: Export configurations for reproducibility and sharing
4. **Independent models**: Each model manages its own hyperparameters - no global state
5. **Validate before building**: Check `model.hyperparams.validate()` before long computations
6. **Use descriptive configs**: Save configurations with meaningful names like `high_res_config.json`

## Common Pitfalls

- ❌ Don't pass `hyperparams=` to `build_system()` unless overriding
- ❌ Don't manually instantiate `PDBModelHyperparameters` - use `model.set_hyperparameters()`
- ✅ Do use `model.set_hyperparameters()` for all configuration
- ✅ Do export/import configurations for reproducibility

## Getting Help

View complete parameter documentation:
```python
from ionerdss.model import pdb
model = pdb.PDBModelBuilder("1ABC")
help(model.set_hyperparameters)
```

See all parameters and their current values:
```python
model.print_hyperparameters()
```

## See Also

- [Migration Guide](../MIGRATION_GUIDE.md) - Updating from old API
- [Dynamic Docstring Documentation](../DYNAMIC_DOCSTRING_README.md) - How documentation is generated
- [Example Scripts](../examples/) - Complete working examples

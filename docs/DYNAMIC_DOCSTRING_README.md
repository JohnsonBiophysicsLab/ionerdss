# Dynamic Docstring Generation

## Overview

The `set_hyperparameters()` function now automatically generates its docstring from the field metadata in the `PDBModelHyperparameters` dataclass. This eliminates the need for manual copy-paste and ensures documentation stays in sync with the actual parameters.

## How It Works

### 1. Field Metadata in hyperparameters.py

Each field in the `PDBModelHyperparameters` dataclass now includes metadata with description and optional unit:

```python
@dataclass
class PDBModelHyperparameters:
    interface_detect_distance_cutoff: float = field(
        default=0.6,
        metadata={
            "description": "Contact search radius per atom pair for interface detection",
            "unit": "nm"
        }
    )
    
    ode_enabled: bool = field(
        default=False,
        metadata={
            "description": "Enable ODE pipeline for kinetic modeling"
        }
    )
```

### 2. Docstring Generation in api.py

The `_generate_hyperparameters_docstring()` function extracts this metadata and builds a complete docstring:

```python
def _generate_hyperparameters_docstring() -> str:
    """Generate docstring from PDBModelHyperparameters field metadata."""
    # Extract field metadata
    field_metadata = {}
    for field_info in fields(PDBModelHyperparameters):
        field_metadata[field_info.name] = {
            'type': field_info.type,
            'default': field_info.default,
            'metadata': field_info.metadata
        }
    
    # Build docstring with parameter descriptions, types, defaults, and units
    # ...
```

### 3. Dynamic Assignment

The docstring is assigned to the function after its definition:

```python
def set_hyperparameters(builder: 'PDBModelBuilder', **kwargs) -> PDBModelHyperparameters:
    # function body
    ...

# Dynamically set docstring from field metadata
set_hyperparameters.__doc__ = _generate_hyperparameters_docstring()
```

## Benefits

✅ **Single Source of Truth**: Parameter descriptions live only in the dataclass  
✅ **No Manual Copy-Paste**: Documentation updates automatically  
✅ **Type Safety**: Types and defaults are guaranteed to match  
✅ **Reduced Maintenance**: Changes to parameters automatically update docs  
✅ **Consistency**: Same format for all parameters  

## Adding New Parameters

To add a new hyperparameter:

1. Add it to `PDBModelHyperparameters` with metadata:
   ```python
   new_parameter: int = field(
       default=42,
       metadata={"description": "Description of what this does", "unit": "optional unit"}
   )
   ```

2. Add the field name to the appropriate category in `_generate_hyperparameters_docstring()`:
   ```python
   categories = {
       "Your Category": [
           "new_parameter",
       ],
       ...
   }
   ```

3. The docstring will automatically include it!

## Example Output

The generated docstring includes:

```
**Core Detection Parameters:**
- interface_detect_distance_cutoff (float, default=0.6): Contact search radius per atom pair for interface detection [nm]
- interface_detect_n_residue_cutoff (int, default=3): Minimum number of contacting residues (on each chain) to accept an interface [residues]

**ODE Pipeline Options:**
- ode_enabled (bool, default=False): Enable ODE pipeline for kinetic modeling
- ode_time_span (tuple, default=(0.0, 10.0)): Time span for ODE solving (start, end) [seconds]
...
```

## Viewing the Documentation

Users can access the complete documentation with:

```python
from ionerdss.model import pdb

builder = pdb.PDBModelBuilder('1ABC')
help(builder.set_hyperparameters)
```

Or:

```python
from ionerdss.model.pdb import api
help(api.set_hyperparameters)
```

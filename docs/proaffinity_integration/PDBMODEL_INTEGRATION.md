# PDBModel Integration with ProAffinity-GNN

## Overview

The `PDBModel` class has been updated to automatically use ProAffinity-GNN for binding energy predictions, with the original additive method as a backup.

## Changes Made

### 1. Import Added

```python
from .proaffinity_predictor import predict_proaffinity_binding_energy
```

### 2. New Method: `_predict_binding_energy_proaffinity()`

**Purpose**: Wrapper method to call ProAffinity-GNN for a chain pair

**Parameters**:
- `chain_id_1` (str): First chain ID
- `chain_id_2` (str): Second chain ID  
- `verbose` (bool): Enable detailed logging

**Returns**:
- `float`: Binding energy in RT units, or `np.nan` if prediction fails

**Features**:
- Automatically handles PDB ID validation
- Formats chain specification for ProAffinity
- Converts kJ/mol to RT units for consistency
- Graceful error handling

### 3. Modified Method: `_build_reactions()`

**Before**:
```python
# calculate the rates
energy = interface_1.energy
```

**After**:
```python
# Predict binding energy: Try ProAffinity first, fall back to additive method
proaffinity_energy = self._predict_binding_energy_proaffinity(
    binding_pair[0], binding_pair[1], verbose=False
)

if not np.isnan(proaffinity_energy):
    energy = proaffinity_energy
    energy_method = "ProAffinity-GNN"
    proaffinity_count += 1
    print(f"Binding pair {binding_pair[0]}-{binding_pair[1]}: Using ProAffinity-GNN (ΔG = {proaffinity_energy:.2f} RT)")
else:
    energy = interface_1.energy
    energy_method = "Additive"
    additive_count += 1
    print(f"Binding pair {binding_pair[0]}-{binding_pair[1]}: Using additive method (ΔG = {energy:.2f} RT)")

# Store which method was used
reaction.energy_method = energy_method
```

**Added Features**:
- Tracking counters for each prediction method
- Per-binding-pair logging showing which method was used
- Summary statistics at the end
- Energy method stored in reaction object

### 4. Summary Output

At the end of `_build_reactions()`, a summary is printed:

```
=== Binding Energy Prediction Summary ===
ProAffinity-GNN predictions: 2/3
Additive method (backup):    1/3
=========================================
```

## Workflow

```
For each binding pair:
  ┌─────────────────────────────┐
  │ Try ProAffinity-GNN first   │
  └──────────┬──────────────────┘
             │
             ├─ Success? ──> Use ProAffinity energy
             │
             └─ Failed? ──> Use additive method (backup)
```

## Unit Conversion

ProAffinity-GNN returns energy in **kJ/mol**, which is converted to **RT units** for consistency with the energy table:

```latex
E_{RT} = \frac{\Delta G_{kJ/mol}}{RT}
```

Where RT = 2.479 kJ/mol at 298.15 K.

## Example Output

```
Binding pair A-B: Using ProAffinity-GNN (ΔG = -14.41 RT)
Binding pair C-D: Using additive method (ΔG = -5.23 RT)
Binding pair E-F: Using ProAffinity-GNN (ΔG = -12.87 RT)

=== Binding Energy Prediction Summary ===
ProAffinity-GNN predictions: 2/3
Additive method (backup):    1/3
=========================================
```

## Advantages

1. **Automatic fallback**: If ProAffinity fails (no PDB ID, prediction error, etc.), additive method is used
2. **Transparent**: Clear logging shows which method was used for each binding pair
3. **Backward compatible**: If ProAffinity is unavailable, code still works with additive method
4. **Trackable**: `reaction.energy_method` attribute allows downstream analysis

## Testing

To test the integration:

```python
from ionerdss.model.pdb_model import PDBModel

# Initialize with PDB ID (enables ProAffinity)
model = PDBModel(pdb_id="1PPE", save_dir="./test_output")

# Run coarse graining and regularization
model.coarse_grain()
model.regularize_repeated_chains()

# Check which method was used
for reaction in model.reaction_list:
    print(f"{reaction.expression}: {reaction.energy_method}")
```

## Notes

- ProAffinity predictions are attempted for all binding pairs
- If no PDB ID is available, additive method is always used
- Energy units are automatically converted to maintain consistency
- The integration is completely non-breaking - existing workflows continue to work

---

**Version**: 1.2.0  
**Date**: 2025-11-04  
**Status**: ✅ Integrated and Tested

# Quick Reference: Affinity Prediction in coarse_grain

## TL;DR

The `coarse_grain` method now supports optional binding affinity prediction using ProAffinity-GNN.

## Basic Usage

### Option 1: Fixed Energy (Default, Fast)
```python
model = PDBModel(pdb_id="8erq", save_dir="./output")
# to use existing pdb files
# model = PDBModel(pdb_file="yourPDB.pdb", save_dir="./output")
model.coarse_grain()
# Uses fixed energy: -16 RT ≈ -39.6 kJ/mol
```

### Option 2: Predicted Energy (Accurate, Slower)
```python
model = PDBModel(pdb_id="8erq", save_dir="./output")
model.coarse_grain(
    predict_affinity=True,
    adfr_path='/path/to/ADFRsuite/bin/prepare_receptor'
)
# Predicts actual binding energies using ProAffinity-GNN
```

## New Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `predict_affinity` | bool | `False` | Enable ProAffinity prediction |
| `adfr_path` | str | `None` | Path to ADFR prepare_receptor |

## Requirements for Prediction

1. **ADFR Suite**: https://ccsb.scripps.edu/adfr/downloads/
2. **ProAffinity Model**: `ionerdss/model/model.pkl`
3. **Dependencies**: Listed in `proaffinity_requirements.txt`

## Quick Comparison

| Feature | Fixed Energy | Predicted Energy |
|---------|--------------|------------------|
| Speed | Fast | Slower (~30s per pair + ~30s for PDBQT conversion) |
| Accuracy | Approximation | ML-predicted |
| Setup | None | Requires ADFR |
| Use Case | Quick testing | Production runs |

## Example Output

With `predict_affinity=True` and `standard_output=True`:

```
Predicted binding energy for chains A-H: -45.23 kJ/mol
Predicted binding energy for chains A-L: -52.17 kJ/mol
...
```

## Access Results

```python
# After coarse_grain()
for i, chain in enumerate(model.all_chains):
    print(f"Chain {chain.id}:")
    print(f"  Partners: {model.all_interfaces[i]}")
    print(f"  Energies: {model.all_interface_energies[i]}")
```

## Error Handling

If prediction fails:
- Automatically falls back to fixed energy (-16 RT)
- Prints warning if `standard_output=True`
- Continues processing other chain pairs

## See Also

- Full documentation: `docs/AFFINITY_PREDICTION.md`
- Example script: `examples/affinity_prediction_example.py`
- Test reference: `tests/test_affinity_prediction.py`

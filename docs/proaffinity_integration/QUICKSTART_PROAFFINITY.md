# ProAffinity-GNN Quick Start Guide

## One-Command Prediction

```python
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

# Predict from PDB ID - that's it!
dG = predict_proaffinity_binding_energy("1PPE", "E,I", verbose=True)
```

## What Happens Automatically

✓ Downloads PDB file from RCSB  
✓ Filters to conventional amino acids  
✓ Converts to PDBQT format (adds hydrogens)  
✓ Runs ProAffinity-GNN model  
✓ Converts result to kJ/mol  
✓ Returns binding energy  

## Function Signature

```python
predict_proaffinity_binding_energy(
    pdb_id,                    # Required: PDB identifier (str)
    chains,                    # Required: Chain spec "A,B" or "AB,CD" (str)
    proaffinity_module_path=None,  # Optional: Path to ProAffinity module
    adfr_path=None,           # Optional: Path to prepare_receptor
    download_dir="pdbfiles",  # Optional: Where to save PDBs
    verbose=False             # Optional: Print progress
)
```

## Chain Specification Format

- **Simple dimer**: `"A,B"` (chain A binding to chain B)
- **Multi-chain**: `"AB,CD"` (chains A+B binding to chains C+D)
- **Symmetric**: `"A,A"` (homodimer)

## Return Value

- **Success**: `float` - Binding energy in kJ/mol (typically -100 to 0)
- **Error**: `np.nan` - Check error messages if verbose=True

## Complete Example

```python
import numpy as np
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

# List of PDB IDs to analyze
pdb_list = ["1PPE", "1ACB", "2PCC"]
chains_list = ["E,I", "E,I", "A,B"]

results = []
for pdb_id, chains in zip(pdb_list, chains_list):
    dG = predict_proaffinity_binding_energy(
        pdb_id=pdb_id,
        chains=chains,
        download_dir="./my_pdbs",
        verbose=True
    )
    results.append(dG)
    print(f"{pdb_id}: {dG:.2f} kJ/mol")

# Filter out failed predictions
valid_results = [r for r in results if not np.isnan(r)]
print(f"\nSuccessfully predicted {len(valid_results)}/{len(results)} structures")
```

## Unit Conversion Reference

ProAffinity-GNN outputs energy in **k_BT** units (thermal energy units).

**Conversion** (at 298.15 K):
- 1 k_BT ≈ 2.479 kJ/mol
- 1 k_BT ≈ 0.593 kcal/mol

The wrapper automatically converts to kJ/mol.

## Advanced: Custom ADFR Path

If ADFR is installed in a non-standard location:

```python
dG = predict_proaffinity_binding_energy(
    pdb_id="1ABC",
    chains="A,B",
    adfr_path="/custom/path/to/prepare_receptor"
)
```

## Batch Processing Template

```python
import pandas as pd
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

# Load dataset
df = pd.read_excel("structures.xlsx")  # Columns: PDB, Chains

# Add predictions
df['Predicted_dG_kJ'] = df.apply(
    lambda row: predict_proaffinity_binding_energy(
        pdb_id=row['PDB'],
        chains=row['Chains'],
        download_dir="./batch_pdbs",
        verbose=False
    ),
    axis=1
)

# Filter successful predictions
df_valid = df[~df['Predicted_dG_kJ'].isna()]

# Save results
df_valid.to_csv("predictions.csv", index=False)
```

## Files Generated

For each prediction (e.g., PDB ID = "1ABC"):

```
pdbfiles/
├── 1abc.pdb        # Downloaded PDB file
└── 1abc.pdbqt      # Converted PDBQT file (for ProAffinity)
```

## Error Handling

The function is designed to fail gracefully:

```python
dG = predict_proaffinity_binding_energy("INVALID", "A,B", verbose=True)
# Prints error message
# Returns: np.nan

# Check for errors:
import numpy as np
if np.isnan(dG):
    print("Prediction failed!")
else:
    print(f"Success: {dG} kJ/mol")
```

## Quick Validation

Validate your installation:

```bash
python test_proaffinity_integration.py
```

Should see:
```
Testing ProAffinity-GNN Integration
...
✓ Prediction successful!
Binding energy: [value] kJ/mol
Integration test PASSED ✓
```

## Need Help?

1. **Installation issues**: Check `PROAFFINITY_INTEGRATION.md`
2. **Usage examples**: See this file
3. **Technical details**: Read `PROAFFINITY_INTEGRATION.md`
4. **Change history**: See `CHANGELOG_PROAFFINITY.md`

---

**Quick Reference Card**
- Input: PDB ID + chains
- Output: ΔG in kJ/mol
- Automatic: Everything!
- Errors: Returns np.nan
- Docs: PROAFFINITY_INTEGRATION.md

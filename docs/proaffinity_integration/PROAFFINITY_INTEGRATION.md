# ProAffinity-GNN Integration Documentation

## Overview
This document tracks the integration of ProAffinity-GNN functionality into the ionerdss package for protein-protein binding energy prediction.

## Changes Made

### 1. ADFR Suite Installation
**Location**: `/home/workspace/GitHub/ionerdss/ADFRsuite/`

- Downloaded and installed ADFR Suite v1.0 for Linux x86_64
- Provides `prepare_receptor` tool for PDB to PDBQT conversion
- Installation size: ~103 MB

### 2. ProAffinity-GNN Files
**Original Location**: `GitHub/ionerdss/proaffinity-gnn/`
**New Location**: `GitHub/ionerdss/ionerdss/model/`

**Files moved**:
- `model.pkl` (27 MB) - Trained ProAffinity-GNN model weights
- `ProAffinity_GNN_inference.py` (28 KB) - Inference module

**Note**: The `proaffinity-gnn` folder will be removed in future versions.

### 3. Wrapper Module Created
**Location**: `ionerdss/model/proaffinity_predictor.py`

Created a self-contained wrapper module with the following functions:

#### Core Functions:

1. **`download_pdb_direct(pdb_id, download_dir, verbose)`**
   - Downloads PDB files directly from RCSB PDB
   - Returns path to downloaded file
   - Handles errors and checks for existing files

2. **`filter_pdb_file(input_pdb_path, output_pdb_path)`**
   - Filters PDB to keep only ATOM records
   - Removes non-conventional residues
   - Keeps TER records for chain termination

3. **`pdb_to_pdbqt(pdbfile, adfr_path, ph, verbose)`**
   - Converts PDB to PDBQT format using ADFR's prepare_receptor
   - Adds hydrogens at specified pH
   - Auto-locates ADFR if path not specified
   - Cleans up temporary files

4. **`kbt_to_kj_mol(dG, T_kelvin)`**
   - Converts energy from k_BT units to kJ/mol
   - Uses fundamental constants (Boltzmann, Avogadro)
   - Handles both scalar and array inputs

5. **`convert_pka_dG(pka, temperature)`**
   - Converts pKa to free energy (dG) in kJ/mol
   - Uses standard thermodynamic relationship

6. **`predict_proaffinity_binding_energy(pdb_id, chains, ...)`**
   - **Main wrapper function - fully self-contained pipeline**
   - Takes PDB ID and chain specification as input
   - Handles all steps automatically:
     - Downloads PDB
     - Filters and curates
     - Converts to PDBQT
     - Runs ProAffinity inference
     - Converts units
   - Returns binding energy in kJ/mol
   - Returns np.nan on error

7. **`run_proaffinity_from_pdbid(pdb_id, chains, **kwargs)`**
   - Simplified convenience wrapper
   - Matches Test.ipynb usage pattern
   - Forwards to main prediction function

## Usage Example

```python
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

# Predict binding affinity from PDB ID
pdb_id = "1abc"
chains = "A,B"  # or "AB,CD" for multi-chain

dG = predict_proaffinity_binding_energy(
    pdb_id=pdb_id,
    chains=chains,
    verbose=True
)

print(f"Binding energy: {dG:.2f} kJ/mol")
```

## Integration with pdb_model.py

The proaffinity predictor can be called from PDBModel class methods. Example integration:

```python
# In pdb_model.py
from .proaffinity_predictor import predict_proaffinity_binding_energy

class PDBModel(Model):
    # ... existing code ...
    
    def predict_binding_affinity_proaffinity(self, chain_pair_string):
        """
        Predict binding affinity using ProAffinity-GNN.
        
        Args:
            chain_pair_string (str): Chain specification (e.g., 'A,B')
            
        Returns:
            float: Binding energy in kJ/mol
        """
        return predict_proaffinity_binding_energy(
            pdb_id=self.pdb_id,
            chains=chain_pair_string,
            download_dir=self.save_dir,
            verbose=True
        )
```

## Dependencies

### Required packages (from proaffinity-gnn/requirements.txt):
- torch
- transformers
- numpy
- biopython
- openbabel (via ADFR)

### System requirements:
- ADFR Suite (installed at `/home/workspace/GitHub/ionerdss/ADFRsuite/`)
- ProAffinity-GNN model (at `/home/workspace/GitHub/ionerdss/proaffinity-gnn/`)

## File Structure

```
GitHub/ionerdss/
├── ADFRsuite/                          # ADFR tools installation
│   └── bin/prepare_receptor            # PDB→PDBQT converter
├── ionerdss/
│   └── model/
│       ├── proaffinity_predictor.py    # New wrapper module
│       ├── ProAffinity_GNN_inference.py # ProAffinity inference (moved here)
│       ├── model.pkl                   # Model weights (moved here)
│       └── pdb_model.py                # Existing PDBModel class
├── proaffinity-gnn/                    # ⚠️ DEPRECATED - will be removed
│   └── Test.ipynb                      # Original test notebook (reference only)
└── Documentation files (*.md)
```

## Notes

1. **ADFR Path**: The module automatically searches for ADFR in common locations:
   - `/home/workspace/GitHub/ionerdss/ADFRsuite/bin/prepare_receptor`
   - `~/ADFRsuite/bin/prepare_receptor`
   - `/usr/local/bin/prepare_receptor`

2. **Error Handling**: All functions include comprehensive error handling and return
   sensible defaults (None or np.nan) on failure.

3. **Temporary Files**: The PDB-to-PDBQT conversion creates a temporary filtered PDB
   file that is automatically cleaned up.

4. **Unit Conversion**: 
   - ProAffinity-GNN outputs binding energy in k_BT units
   - The wrapper automatically converts to kJ/mol
   - Conversion factor at 298.15 K: 1 k_BT ≈ 2.479 kJ/mol

## Testing

To test the integration:

```python
# Test with a known PDB structure
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

result = predict_proaffinity_binding_energy(
    pdb_id="1PPE",
    chains="E,I",
    verbose=True
)

print(f"Result: {result} kJ/mol")
```

## Future Enhancements

1. Add caching mechanism for downloaded PDB files
2. Support for custom PDBQT files (skip download/conversion steps)
3. Batch prediction for multiple structures
4. Integration of ionerdss and proaffinity predictions for ensemble methods
5. Add option to keep/delete intermediate files

## References

- ADFR Suite: https://ccsb.scripps.edu/adfr/
- ProAffinity-GNN: Based on the implementation in proaffinity-gnn directory
- RCSB PDB: https://www.rcsb.org/

# ProAffinity-GNN Integration for ionerdss

## Quick Start

This integration allows you to predict protein-protein binding affinities using ProAffinity-GNN directly from a PDB ID.

### Basic Usage

```python
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

# Predict binding energy from PDB ID
binding_energy = predict_proaffinity_binding_energy(
    pdb_id="1PPE",
    chains="E,I",
    verbose=True
)

print(f"Predicted ΔG: {binding_energy:.2f} kJ/mol")
```

## What's Been Integrated

### 1. ADFR Suite ✓
- **Location**: `./ADFRsuite/`
- **Purpose**: PDB to PDBQT conversion (adds hydrogens, prepares for docking)
- **Tool**: `prepare_receptor`

### 2. ProAffinity Wrapper Module ✓
- **File**: `ionerdss/model/proaffinity_predictor.py`
- **Functions**:
  - Download PDB files
  - Filter/curate PDB structures
  - Convert to PDBQT format
  - Run ProAffinity-GNN inference
  - Convert units (k_BT → kJ/mol)

### 3. Complete Documentation ✓
- `PROAFFINITY_INTEGRATION.md` - Detailed technical documentation
- `CHANGELOG_PROAFFINITY.md` - Version history and changes
- This file - Quick reference guide

## Installation Checklist

- [x] ADFR Suite downloaded and installed
- [x] proaffinity_predictor.py module created
- [x] Documentation written
- [x] Test script created
- [x] All helper functions implemented

## Pipeline Overview

```
PDB ID → Download → Filter → PDBQT → ProAffinity → kJ/mol
         (RCSB)    (ATOM)   (ADFR)   (Inference)  (Convert)
```

**Input**: PDB ID (e.g., "1ABC") + chain specification (e.g., "A,B")  
**Output**: Binding energy in kJ/mol

## Testing

Run the included test script:

```bash
cd /home/workspace/GitHub/ionerdss
python test_proaffinity_integration.py
```

Expected output:
- ✓ PDB download successful
- ✓ PDBQT conversion successful  
- ✓ ProAffinity prediction successful
- ✓ Energy in kJ/mol returned

## Detailed Documentation

- **Technical Details**: See `PROAFFINITY_INTEGRATION.md`
- **Change Log**: See `CHANGELOG_PROAFFINITY.md`
- **Test Script**: `test_proaffinity_integration.py`

## Requirements

### Python Packages
- numpy
- biopython
- torch
- transformers

### System Components
- ADFR Suite (installed at `./ADFRsuite/`)
- ProAffinity-GNN model (at `./proaffinity-gnn/`)

## Example: Integration with PDBModel

You can easily add ProAffinity prediction as a method in the PDBModel class:

```python
# In ionerdss/model/pdb_model.py

from .proaffinity_predictor import predict_proaffinity_binding_energy

class PDBModel(Model):
    # ... existing methods ...
    
    def predict_binding_with_proaffinity(self, chain_specification):
        """
        Predict binding affinity using ProAffinity-GNN.
        
        Args:
            chain_specification (str): e.g., 'A,B' or 'AB,CD'
            
        Returns:
            float: Binding energy in kJ/mol
        """
        return predict_proaffinity_binding_energy(
            pdb_id=self.pdb_id,
            chains=chain_specification,
            download_dir=self.save_dir,
            verbose=False
        )
```

Then use it:

```python
from ionerdss.model.pdb_model import PDBModel

model = PDBModel(pdb_id="1PPE", save_dir="./output")
dG = model.predict_binding_with_proaffinity("E,I")
print(f"Binding energy: {dG:.2f} kJ/mol")
```

## File Locations

```
ionerdss/
├── ADFRsuite/                              # ADFR installation
│   ├── bin/prepare_receptor                # PDB→PDBQT converter
│   └── ...
├── proaffinity-gnn/                        # ProAffinity model
│   ├── ProAffinity_GNN_inference.py        # Inference code
│   ├── model.pkl                           # Trained model
│   └── requirements.txt                    # Dependencies
├── ionerdss/model/
│   ├── pdb_model.py                        # Existing ionerdss code
│   └── proaffinity_predictor.py            # NEW: ProAffinity wrapper
├── PROAFFINITY_INTEGRATION.md              # Technical documentation
├── CHANGELOG_PROAFFINITY.md                # Version history
├── README_PROAFFINITY.md                   # This file
└── test_proaffinity_integration.py         # Test script
```

## Troubleshooting

### Common Issues

**Problem**: `prepare_receptor` not found  
**Solution**: Check ADFR installation at `./ADFRsuite/bin/prepare_receptor`

**Problem**: ProAffinity module import error  
**Solution**: Ensure `./proaffinity-gnn/` contains `ProAffinity_GNN_inference.py` and `model.pkl`

**Problem**: PDB download fails  
**Solution**: Check internet connection and verify PDB ID exists at RCSB

**Problem**: Unit conversion seems wrong  
**Solution**: Verify temperature is set correctly (default 298.15 K)

## Support

For issues specific to:
- **ADFR**: See ADFR Suite documentation
- **ProAffinity-GNN**: Check proaffinity-gnn directory
- **Integration**: See `PROAFFINITY_INTEGRATION.md`
- **ionerdss**: Refer to main ionerdss documentation

---

**Last Updated**: 2025-11-04  
**Version**: 1.0.0  
**Status**: Initial Release ✓

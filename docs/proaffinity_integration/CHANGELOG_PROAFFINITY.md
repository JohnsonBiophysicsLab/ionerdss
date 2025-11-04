# ProAffinity Integration Changelog

## Version 1.2.0 - 2025-11-04

### Added - PDBModel Integration

#### ProAffinity-GNN Now Used Automatically in PDBModel
**Modified File**: `ionerdss/model/pdb_model.py`

**What Changed**:
The `PDBModel` class now automatically uses ProAffinity-GNN for binding energy predictions, with the original Miyazawa-Jernigan additive method as a backup.

**New Method**: `_predict_binding_energy_proaffinity(chain_id_1, chain_id_2, verbose=False)`
- Wrapper to call ProAffinity-GNN for a chain pair
- Automatically converts kJ/mol to RT units
- Returns `np.nan` on failure (triggers fallback)
- Includes verbose logging option

**Modified Method**: `_build_reactions()`
- Now tries ProAffinity-GNN first for each binding pair
- Falls back to additive method if ProAffinity fails or returns NaN
- Tracks which method was used via `reaction.energy_method` attribute
- Prints per-pair results showing which method was used
- Displays summary statistics at the end

#### Example Output
```
Binding pair A-B: Using ProAffinity-GNN (ΔG = -14.41 RT)
Binding pair C-D: Using additive method (ΔG = -5.23 RT)

=== Binding Energy Prediction Summary ===
ProAffinity-GNN predictions: 1/2
Additive method (backup):    1/2
=========================================
```

#### Benefits
- ✅ More accurate predictions when ProAffinity-GNN works
- ✅ Automatic fallback ensures reliability
- ✅ Transparent logging shows which method was used
- ✅ Fully backward compatible
- ✅ No breaking changes to existing code

---

## Version 1.1.0 - 2025-11-04

### Changed - File Organization

#### ProAffinity-GNN Files Relocated
**Moved from**: `GitHub/ionerdss/proaffinity-gnn/`  
**Moved to**: `GitHub/ionerdss/ionerdss/model/`

**Files relocated**:
- `model.pkl` (27 MB) - ProAffinity-GNN trained model weights
- `ProAffinity_GNN_inference.py` (28 KB) - Core inference module

**Rationale**: 
- Better integration with ionerdss package structure
- Cleaner imports (relative imports from same module)
- The `proaffinity-gnn/` folder is now deprecated and will be removed in future versions

#### Updated Import Paths
**Old approach** (v1.0.0):
```python
# Had to add proaffinity-gnn to sys.path
sys.path.insert(0, '/path/to/proaffinity-gnn')
from ProAffinity_GNN_inference import run_proaffinity_inference
```

**New approach** (v1.1.0):
```python
# Clean relative import from same package
from .ProAffinity_GNN_inference import run_proaffinity_inference
```

#### Modified Files
1. **`ionerdss/model/proaffinity_predictor.py`**:
   - Changed from absolute import to relative import
   - Updated default model weights path to use `__file__` directory
   - Removed `proaffinity_module_path` parameter
   - Added `model_weights_path` parameter instead
   - Now automatically locates model.pkl in same directory

2. **`test_proaffinity_integration.py`**:
   - Removed references to proaffinity-gnn folder
   - Simplified to use direct import from ionerdss.model

3. **Documentation files**:
   - Updated all file paths in documentation
   - Added deprecation notices for proaffinity-gnn folder
   - Updated file structure diagrams

### Migration Guide

If you have existing code using v1.0.0:

**Before (v1.0.0)**:
```python
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

dG = predict_proaffinity_binding_energy(
    pdb_id="1PPE",
    chains="E,I",
    proaffinity_module_path="/path/to/proaffinity-gnn"  # OLD
)
```

**After (v1.1.0)**:
```python
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

dG = predict_proaffinity_binding_energy(
    pdb_id="1PPE",
    chains="E,I",
    model_weights_path=None  # NEW - auto-detects, or specify custom path
)
```

### Benefits
- ✅ Cleaner package structure
- ✅ No need to modify sys.path
- ✅ Easier to distribute as part of ionerdss package
- ✅ Relative imports work properly
- ✅ Model weights automatically located

---

## Version 1.0.0 - 2025-11-04

### Added

#### New Module: `ionerdss/model/proaffinity_predictor.py`
- Created self-contained wrapper for ProAffinity-GNN binding affinity prediction
- Implements complete pipeline from PDB ID to binding energy prediction

#### Functions Added:
1. `download_pdb_direct()` - Direct PDB file download from RCSB
2. `filter_pdb_file()` - Filter PDB to conventional residues only
3. `pdb_to_pdbqt()` - Convert PDB to PDBQT using ADFR's prepare_receptor
4. `kbt_to_kj_mol()` - Convert energy from k_BT to kJ/mol units
5. `convert_pka_dG()` - Convert pKa to free energy in kJ/mol
6. `predict_proaffinity_binding_energy()` - Main prediction pipeline
7. `run_proaffinity_from_pdbid()` - Simplified wrapper interface

### Dependencies Installed

#### ADFR Suite v1.0
- **Location**: `/home/workspace/GitHub/ionerdss/ADFRsuite/`
- **Key Tool**: `prepare_receptor` for PDB to PDBQT conversion
- **License**: Academic use (as per installation)
- **Downloaded from**: https://ccsb.scripps.edu/adfr/

### Documentation

#### Created Files:
1. `PROAFFINITY_INTEGRATION.md` - Complete integration documentation
   - Overview of integration
   - Function descriptions
   - Usage examples
   - File structure
   - Testing instructions
   - Future enhancements

2. `CHANGELOG_PROAFFINITY.md` - This file
   - Version tracking
   - Change documentation
   - Dependencies list

3. `test_proaffinity_integration.py` - Test script
   - Validates installation
   - Tests basic prediction
   - Checks result validity

### Implementation Details

#### Pipeline Steps:
1. **PDB Download**: Uses urllib to fetch from RCSB PDB
2. **PDB Curation**: Filters ATOM records, removes non-standard residues
3. **PDBQT Conversion**: Uses ADFR's prepare_receptor with hydrogen addition
4. **Inference**: Calls ProAffinity_GNN_inference.run_proaffinity_inference()
5. **Unit Conversion**: Converts output from k_BT to kJ/mol

#### Key Features:
- Automatic ADFR path detection
- Comprehensive error handling
- Verbose logging option
- Automatic cleanup of temporary files
- Returns np.nan on errors for easy filtering

### Usage

#### Basic Usage:
```python
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

dG = predict_proaffinity_binding_energy(
    pdb_id="1ABC",
    chains="A,B",
    verbose=True
)
```

#### Integration with PDBModel:
Users can add a method to the PDBModel class to call ProAffinity predictions:

```python
def predict_binding_affinity_proaffinity(self, chain_pair_string):
    from .proaffinity_predictor import predict_proaffinity_binding_energy
    return predict_proaffinity_binding_energy(
        pdb_id=self.pdb_id,
        chains=chain_pair_string,
        download_dir=self.save_dir
    )
```

### Testing

Run the test script:
```bash
cd /home/workspace/GitHub/ionerdss
python test_proaffinity_integration.py
```

Expected output:
- Downloads PDB file
- Converts to PDBQT
- Runs ProAffinity inference
- Returns binding energy in kJ/mol

### Requirements

#### Python Packages (from proaffinity-gnn/requirements.txt):
- torch
- transformers  
- numpy
- biopython
- Additional packages as specified in ProAffinity-GNN requirements

#### System Requirements:
- ADFR Suite (installed)
- ProAffinity-GNN model files (model.pkl)
- Sufficient disk space for PDB files

### Known Limitations

1. Requires ProAffinity-GNN model files to be present
2. ADFR uses Python 2.7 internally but is called as subprocess
3. May require adjustments for different ProAffinity-GNN versions
4. Download speed depends on RCSB PDB server

### Future Work

#### Planned Enhancements:
- [ ] Add caching for downloaded PDB files
- [ ] Support batch predictions
- [ ] Add ensemble predictions (ionerdss + ProAffinity)
- [ ] Create visualization tools for predictions
- [ ] Add confidence metrics
- [ ] Support for custom PDB files (not just PDB IDs)
- [ ] Add multi-chain complex handling
- [ ] Create example notebooks demonstrating integration

#### Integration Opportunities:
- Combine with ionerdss coarse-graining for multi-scale analysis
- Use ProAffinity for validation of ionerdss predictions
- Create hybrid scoring function
- Add to existing PDBModel workflow

### Notes

- All unit conversions verified against Test.ipynb
- Error handling ensures graceful failures
- Verbose mode provides detailed progress information
- Temporary files automatically cleaned up
- Compatible with existing ionerdss architecture

---

## Maintenance

### Contact
For issues or questions about this integration, please refer to:
- ProAffinity-GNN repository documentation
- ADFR Suite documentation
- ionerdss main repository

### Version History
- v1.0.0 (2025-11-04): Initial integration
  - Core functionality implemented
  - ADFR installed
  - Documentation created
  - Test script added

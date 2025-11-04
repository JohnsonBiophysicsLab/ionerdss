# Migration Guide: v1.0.0 → v1.1.0

## Summary of Changes

ProAffinity-GNN files have been **relocated** from `proaffinity-gnn/` to `ionerdss/model/` for better package integration.

## What Changed?

### File Locations

| File | v1.0.0 Location | v1.1.0 Location | Status |
|------|----------------|----------------|---------|
| `model.pkl` | `proaffinity-gnn/` | `ionerdss/model/` | ✅ Moved |
| `ProAffinity_GNN_inference.py` | `proaffinity-gnn/` | `ionerdss/model/` | ✅ Moved |
| `Test.ipynb` | `proaffinity-gnn/` | `proaffinity-gnn/` | ⚠️ Deprecated folder |

### Import Changes

**v1.0.0** (absolute import with sys.path manipulation):
```python
import sys
sys.path.insert(0, '/path/to/proaffinity-gnn')
from ProAffinity_GNN_inference import run_proaffinity_inference
```

**v1.1.0** (clean relative import):
```python
# Inside ionerdss/model/proaffinity_predictor.py
from .ProAffinity_GNN_inference import run_proaffinity_inference
```

### API Changes

#### Function Parameter Update

**Old** (v1.0.0):
```python
predict_proaffinity_binding_energy(
    pdb_id="1PPE",
    chains="E,I",
    proaffinity_module_path="/path/to/proaffinity-gnn"  # REMOVED
)
```

**New** (v1.1.0):
```python
predict_proaffinity_binding_energy(
    pdb_id="1PPE",
    chains="E,I",
    model_weights_path=None  # NEW - optional, auto-detects
)
```

## Do I Need to Change My Code?

### If you're using the high-level API: **NO CHANGES NEEDED** ✅

The main function signature is **backward compatible**:

```python
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

# This works in both v1.0.0 and v1.1.0
dG = predict_proaffinity_binding_energy("1PPE", "E,I", verbose=True)
```

### If you were specifying `proaffinity_module_path`: **MINOR UPDATE NEEDED** ⚠️

**Before** (v1.0.0):
```python
dG = predict_proaffinity_binding_energy(
    pdb_id="1PPE",
    chains="E,I",
    proaffinity_module_path="/custom/path/to/proaffinity-gnn"
)
```

**After** (v1.1.0):
```python
dG = predict_proaffinity_binding_energy(
    pdb_id="1PPE",
    chains="E,I",
    model_weights_path="/custom/path/to/model.pkl"  # Point to model.pkl directly
)
```

### If you were importing ProAffinity directly: **UPDATE NEEDED** ⚠️

**Before** (v1.0.0):
```python
import sys
sys.path.insert(0, '/home/workspace/GitHub/ionerdss/proaffinity-gnn')
from ProAffinity_GNN_inference import run_proaffinity_inference

dG = run_proaffinity_inference("file.pdbqt", "A,B", weights_path="./model.pkl")
```

**After** (v1.1.0):
```python
from ionerdss.model.ProAffinity_GNN_inference import run_proaffinity_inference

dG = run_proaffinity_inference(
    "file.pdbqt", 
    "A,B", 
    weights_path="/path/to/ionerdss/model/model.pkl"
)
```

## Benefits of v1.1.0

1. ✅ **Cleaner imports**: No more sys.path manipulation
2. ✅ **Better organization**: All model files in one place
3. ✅ **Easier distribution**: Works as a proper Python package
4. ✅ **Auto-detection**: Model weights automatically located
5. ✅ **Future-proof**: Prepares for proaffinity-gnn folder removal

## Deprecation Notice

⚠️ **The `proaffinity-gnn/` folder is deprecated and will be removed in a future version.**

### What to do:
1. Update your code as described above
2. Do not rely on files in `proaffinity-gnn/` folder
3. Use imports from `ionerdss.model` instead

### Timeline:
- **v1.1.0** (current): proaffinity-gnn/ folder still exists but is not used
- **v2.0.0** (future): proaffinity-gnn/ folder will be removed

## Verification

### Check Your Setup

Run this test to verify everything works:

```python
from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy

# Should work without errors
dG = predict_proaffinity_binding_energy("1PPE", "E,I", verbose=True)
print(f"Result: {dG:.2f} kJ/mol")
```

Or run the test script:
```bash
python /home/workspace/GitHub/ionerdss/test_proaffinity_integration.py
```

Expected: "Integration test PASSED ✓"

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'ProAffinity_GNN_inference'"

**Cause**: Old code trying to import from proaffinity-gnn folder

**Solution**: Update imports as shown above

---

### Issue: "Model weights not found"

**Cause**: Custom model_weights_path pointing to old location

**Solution**: Update path to point to `ionerdss/model/model.pkl`

```python
# Old
model_weights_path="/path/to/proaffinity-gnn/model.pkl"

# New
model_weights_path="/path/to/ionerdss/model/model.pkl"
```

---

### Issue: Code works but gives deprecation warning

**Cause**: Using `proaffinity_module_path` parameter (deprecated)

**Solution**: Switch to `model_weights_path` parameter

---

## Questions?

- **Technical docs**: See `PROAFFINITY_INTEGRATION.md`
- **Changelog**: See `CHANGELOG_PROAFFINITY.md`
- **Quick start**: See `QUICKSTART_PROAFFINITY.md`

---

**Migration Date**: 2025-11-04  
**From Version**: 1.0.0  
**To Version**: 1.1.0  
**Breaking Changes**: None (backward compatible with minor deprecations)

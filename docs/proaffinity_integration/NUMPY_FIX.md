# NumPy Compatibility Fix

## Issue

When running ProAffinity-GNN with the original requirements.txt, you may encounter this error:

```
A module that was compiled using NumPy 1.x cannot be run in
NumPy 2.1.2 as it may crash. To support both 1.x and 2.x
versions of NumPy, modules must be compiled with NumPy 2.0.
```

## Root Cause

- **PyTorch 2.2.2** was compiled with **NumPy 1.x**
- The original requirements.txt specified `numpy==1.24.1`
- However, during installation, pip may automatically install NumPy 2.x as a dependency of other packages
- NumPy 2.x has breaking changes that are incompatible with PyTorch 2.2.2

## Solution

**Downgrade NumPy to version < 2.0**

```bash
conda activate proaffinity
pip install "numpy<2.0"
```

This will install NumPy 1.26.4, which is compatible with PyTorch 2.2.2.

## Verification

After the fix, verify that everything works:

```bash
conda activate proaffinity

# Check NumPy version
python -c "import numpy; print(f'NumPy version: {numpy.__version__}')"
# Should output: NumPy version: 1.26.4

# Check PyTorch works
python -c "import torch; print(f'PyTorch version: {torch.__version__}')"
# Should output: PyTorch version: 2.2.2+cu121

# Run the integration test
python test_proaffinity_integration.py
# Should pass without NumPy errors
```

## Updated Requirements

The corrected requirements should pin NumPy to < 2.0:

```txt
numpy<2.0  # Instead of numpy==1.24.1
torch==2.2.2
torchvision==0.17.2
torchaudio==2.2.2
torch_geometric==2.3.0
transformers==4.38
scikit-learn
scipy
```

## Status

✅ **FIXED** - The proaffinity conda environment now has NumPy 1.26.4 installed and all tests pass.

---

**Issue Discovered**: 2025-11-04  
**Resolution**: Downgrade to numpy<2.0  
**Status**: ✅ Resolved

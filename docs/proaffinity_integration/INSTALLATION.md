# ProAffinity-GNN Installation Guide

## Prerequisites

- Miniconda or Anaconda installed
- Python 3.10
- CUDA 12.1 (optional, for GPU support)

## Installation Steps

### 1. Create Conda Environment

```bash
conda create -n proaffinity python=3.10 -y
conda activate proaffinity
```

### 2. Install PyTorch with CUDA Support

**For CUDA 12.1** (recommended):
```bash
pip install torch==2.2.2 torchvision==0.17.2 torchaudio==2.2.2 --index-url https://download.pytorch.org/whl/cu121
```

**For CPU only**:
```bash
pip install torch==2.2.2 torchvision==0.17.2 torchaudio==2.2.2 --index-url https://download.pytorch.org/whl/cpu
```

**For other CUDA versions**:
Visit https://pytorch.org/get-started/previous-versions/ and find the appropriate command for your CUDA version.

### 3. Install NumPy (Compatible Version)

**IMPORTANT**: Must be < 2.0 for PyTorch 2.2.2 compatibility

```bash
pip install "numpy<2.0"
```

### 4. Install Remaining Dependencies

```bash
pip install torch_geometric==2.3.0 transformers==4.38 scikit-learn scipy
```

### 5. Verify Installation

```bash
# Test imports
python -c "import torch; print(f'PyTorch: {torch.__version__}')"
python -c "import numpy; print(f'NumPy: {numpy.__version__}')"
python -c "import torch_geometric; print('torch_geometric OK')"
python -c "import transformers; print('transformers OK')"

# Run integration test
cd /home/workspace/GitHub/ionerdss
python test_proaffinity_integration.py
```

## Alternative: One-Command Installation

```bash
conda create -n proaffinity python=3.10 -y && \
conda activate proaffinity && \
pip install torch==2.2.2 torchvision==0.17.2 torchaudio==2.2.2 --index-url https://download.pytorch.org/whl/cu121 && \
pip install "numpy<2.0" && \
pip install torch_geometric==2.3.0 transformers==4.38 scikit-learn scipy
```

## Troubleshooting

### Issue: NumPy Version Error

If you see:
```
A module that was compiled using NumPy 1.x cannot be run in NumPy 2.1.2
```

**Fix**:
```bash
pip install "numpy<2.0"
```

### Issue: CUDA Not Available

If `torch.cuda.is_available()` returns `False`:
1. Check your CUDA installation: `nvidia-smi`
2. Ensure you installed the CUDA version of PyTorch (cu121, not cpu)
3. The model can still run on CPU (slower)

### Issue: Model Download Errors

If ESM model download fails:
```bash
# Set up Hugging Face cache
export HF_HOME=/home/workspace/.cache/huggingface
pip install --upgrade huggingface-hub
```

## Package Versions (Verified)

```
Python: 3.10
PyTorch: 2.2.2+cu121
NumPy: 1.26.4 (< 2.0)
torch_geometric: 2.3.0
transformers: 4.38
scikit-learn: 1.7.2
scipy: 1.15.3
```

## System Requirements

- **RAM**: Minimum 8GB, Recommended 16GB+
- **Storage**: ~5GB for conda environment + models
- **GPU**: Optional, but recommended for faster inference
  - Tested with CUDA 12.1
  - Model size: ~27MB
  - ESM model: ~1.3GB (downloaded on first use)

## Next Steps

After installation:
1. Read: file `'START_HERE.md'`
2. Try: file `'../QUICKSTART_PROAFFINITY.md'`
3. Test: `python test_proaffinity_integration.py`

---

**Last Updated**: 2025-11-04  
**Tested On**: Python 3.10, PyTorch 2.2.2, CUDA 12.1  
**Status**: ✅ Verified Working

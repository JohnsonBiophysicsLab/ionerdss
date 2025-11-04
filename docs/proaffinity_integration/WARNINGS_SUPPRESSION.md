# Warnings Suppression Guide

## Overview

The ProAffinity-GNN integration uses several deep learning libraries that may emit benign warnings. This document explains these warnings and how to suppress them.

## Known Benign Warnings

### 1. HuggingFace Hub Warning
```
FutureWarning: `resume_download` is deprecated and will be removed in version 1.0.0.
```

**Explanation**: This warning is from huggingface_hub and indicates a deprecated parameter. It doesn't affect functionality.

**Status**: ✅ Suppressed in code

### 2. ESM Model Weights Warning
```
Some weights of EsmModel were not initialized from the model checkpoint...
You should probably TRAIN this model on a down-stream task...
```

**Explanation**: ESM2 model has some weights that aren't used for our specific inference task. This is expected and doesn't affect predictions.

**Status**: ⚠️ Partially suppressed (printed during model loading)

## Suppression Methods

### Method 1: In-Code Suppression (Current Approach)

Both `proaffinity_predictor.py` and `ProAffinity_GNN_inference.py` include:

```python
import warnings
warnings.filterwarnings('ignore', category=FutureWarning, module='huggingface_hub.file_download')
warnings.filterwarnings('ignore', message='.*Some weights of EsmModel were not initialized.*')
warnings.filterwarnings('ignore', message='.*You should probably TRAIN this model.*')

import logging
logging.getLogger("transformers").setLevel(logging.ERROR)
```

### Method 2: Environment Variables (Recommended for Scripts)

Set these before importing:

```python
import os
os.environ['TRANSFORMERS_VERBOSITY'] = 'error'
os.environ['TRANSFORMERS_NO_ADVISORY_WARNINGS'] = '1'
```

Or in bash:

```bash
export TRANSFORMERS_VERBOSITY=error
export TRANSFORMERS_NO_ADVISORY_WARNINGS=1
python your_script.py
```

### Method 3: Redirect stderr (For Clean Output)

```bash
python your_script.py 2>/dev/null
```

Or in Python:

```python
import sys
import os

# Redirect stderr temporarily
class SuppressStderr:
    def __enter__(self):
        self.old_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')
        return self
    
    def __exit__(self, *args):
        sys.stderr.close()
        sys.stderr = self.old_stderr

with SuppressStderr():
    from ionerdss.model.proaffinity_predictor import predict_proaffinity_binding_energy
    result = predict_proaffinity_binding_energy('1PPE', 'E,I')
```

## Why Some Warnings Persist

The ESM model warning is printed by transformers **during model initialization**, before our warning filters can intercept it. This is a limitation of how transformers logs these messages.

## Recommendation

For production use, we recommend:

1. **Accept the warnings**: They're benign and don't affect results
2. **Use environment variables**: Set `TRANSFORMERS_VERBOSITY=error` 
3. **Redirect output**: Use `2>/dev/null` if you need clean output

## Impact on Functionality

✅ **None** - All warnings are purely informational and don't affect:
- Prediction accuracy
- Model performance
- Reliability
- Results

The warnings simply inform you that:
1. Some API methods will change in future versions
2. Some model weights aren't used (expected for our use case)

---

**Last Updated**: 2025-11-04  
**Status**: Warnings are benign and can be safely ignored

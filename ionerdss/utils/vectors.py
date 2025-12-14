"""
ionerdss.utils.vectors

Uitlity functions to process vectors in the form of numpy arrays
"""

import numpy as np
from numpy.typing import ArrayLike, NDArray

def get_magnitude(vector: ArrayLike) -> float:
    """
    Compute the Euclidean magnitude of a vector.

    Parameters
    ----------
    vector : ArrayLike
        Any array-like object (list, tuple, or np.ndarray) of numeric values.

    Returns
    -------
    float
        The Euclidean norm of the vector.
    """
    arr = np.asarray(vector, dtype=float)
    return float(np.sqrt(np.sum(arr ** 2)))


def convert_to_unit(vector: ArrayLike, tol: float = 1e-12) -> NDArray[np.float64]:
    """
    Normalize a vector to unit length. Returns a zero vector if norm < tol.

    Parameters
    ----------
    vector : ArrayLike
        Any array-like object (list, tuple, or np.ndarray) of numeric values.
    tol : float, optional
        Small threshold below which the vector is treated as zero.

    Returns
    -------
    np.ndarray
        Normalized vector as a NumPy float64 array.
    """
    arr = np.asarray(vector, dtype=float)
    norm = np.sqrt(np.sum(arr ** 2))
    if norm < tol:
        return np.zeros_like(arr, dtype=float)
    return arr / norm
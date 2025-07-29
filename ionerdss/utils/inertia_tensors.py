"""
inertia_tools.py

This module provides utility functions for analyzing geometric and rotational
symmetries of point-based structures (e.g., molecular assemblies) using inertia tensors
and vector algebra. These tools are particularly useful for identifying symmetry axes
and classifying point group symmetries in 3D configurations such as molecular complexes,
polyhedral assemblies, and coarse-grained models.

Note: This file includes code adapted from the 'pointgroup' package,
originally authored by Abel Carreras (https://github.com/abelcarreras/pointgroup),
and is licensed under the MIT License (which is attached to the end of this docstring)

Included Functions
------------------

1. get_inertia_tensor(coords, tol=1e-12)
    Computes the normalized moment of inertia tensor for a set of 3D points,
    assuming uniform mass. Returns a symmetric 3×3 tensor that characterizes
    how mass is distributed relative to the center of mass.

2. get_degeneracy(eigenvalues, tolerance=0.1)
    Estimates the degeneracy (repeatedness) of eigenvalues within a specified tolerance.
    Useful for inferring rotational symmetry based on isotropy of the inertia tensor.

3. get_non_degenerated(eigenvalues, tolerance=0.1)
    Identifies which eigenvalue in a set of three is non-degenerate, assuming
    two are nearly equal. Often used to detect symmetry axes.

4. get_perpendicular_vector(vector, normalize=True, tol=1e-8)
    Computes a vector orthogonal to a given input vector. Supports general n-dimensional input.
    For 3D vectors, uses the cross product; for nD, falls back to Gram-Schmidt projection.

Use Cases
---------
- Classifying molecular shapes and point group symmetries from 3D atomic coordinates
- Determining principal axes of rotation
- Automatically choosing reference frames based on geometry
- Computing orientation-aware quantities from molecular subunits

Dependencies
------------
- numpy

Examples
--------
>>> coords = np.array([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0]])
>>> I = get_inertia_tensor(coords)
>>> eigvals, eigvecs = np.linalg.eigh(I)
>>> degeneracy = get_degeneracy(eigvals)
>>> main_axis_idx = get_non_degenerated(eigvals)
>>> perp = get_perpendicular_vector(eigvecs[:, main_axis_idx])




The MIT License (MIT)

Copyright (c) 2023 Efrem Bernuz and Abel Carreras

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import numpy as np

def get_inertia_tensor(coords, tol = 1e-12):
    """
    Compute the normalized moment of inertia tensor assuming uniform mass for all sites.

    This function calculates the inertia tensor using the standard formula for point masses:
        I = Σ [‖r_i‖² * I_3 - r_i ⊗ r_i]
    where:
        - r_i is the position vector of the i-th point (assumed relative to the center of mass),
        - I_3 is the 3×3 identity matrix,
        - ⊗ denotes the outer product.

    The result is then normalized by the scalar sum:
        Σ ‖r_i‖²

    The moment of inertia tensor is important in determining
    symmetry because it encodes how mass (or structure) is
    distributed relative to an origin—typically the center of mass.
    Its eigenvalues and eigenvectors reveal rotational properties
    that are tightly linked to the object's symmetry group.
    
    Parameters
    ----------
    coords : ndarray of shape (N, 3)
        Cartesian coordinates of the points (atoms), assumed centered at origin.

    Returns
    -------
    inertia_tensor : ndarray of shape (3, 3)
        Normalized inertia tensor (unitless, assumes equal mass per atom).
    """
    coords = np.asarray(coords)
    inertia_tensor = np.zeros((3, 3))
    total_inertia = 0.0

    for c in coords:
        r2 = np.dot(c, c)
        inertia_tensor += np.identity(3) * r2 - np.outer(c, c)
        total_inertia += r2

    if abs(total_inertia) > tol:
        inertia_tensor /= total_inertia

    return inertia_tensor

def get_degeneracy(eigenvalues, tolerance=0.1):
    """
    Estimate the degeneracy of eigenvalues within a specified tolerance.

    Degeneracy refers to the number of times the same (or nearly the same) eigenvalue appears.
    This function scans through the list of eigenvalues and returns the maximum count of
    eigenvalues that are equal within the given numerical tolerance.

    Parameters
    ----------
    eigenvalues : list or array-like of float
        A sequence of eigenvalues (typically from an inertia tensor or other symmetric matrix).
    
    tolerance : float, optional
        Numerical tolerance used to determine whether two eigenvalues are considered equal.
        Default is 0.1.

    Returns
    -------
    int
        The estimated degeneracy (i.e., how many eigenvalues are approximately equal).
        Returns 1 if no two eigenvalues are close enough to be considered degenerate.
    
    Examples
    --------
    >>> get_degeneracy([1.0, 1.0, 2.0])
    2

    >>> get_degeneracy([1.0, 1.2, 1.4], tolerance=0.01)
    1
    """
    for ev1 in eigenvalues:
        single_deg = 0
        for ev2 in eigenvalues:
            if abs(ev1 - ev2) < tolerance:
                single_deg += 1
        if single_deg > 1:
            return single_deg
    return 1

def get_non_degenerated(eigenvalues, tolerance=0.1):
    """
    Identify the index of a non-degenerate eigenvalue in a set of eigenvalues.

    This function assumes that two of the eigenvalues are approximately equal (degenerate),
    and one is different (non-degenerate). It returns the index of the non-degenerate value.

    Parameters
    ----------
    eigenvalues : list or array-like of float
        A list or array of eigenvalues (typically of an inertia tensor).
        Expected to contain exactly 3 values.

    tolerance : float, optional
        Numerical tolerance for considering two eigenvalues as equal (degenerate). Default is 0.1.

    Returns
    -------
    int
        The index of the non-degenerate eigenvalue in the input list.

    Raises
    ------
    Exception
        If no non-degenerate eigenvalue is found.

    Examples
    --------
    >>> get_non_degenerated([1.0, 1.0, 2.0])
    2
    >>> get_non_degenerated([2.0, 1.0, 2.0])
    1
    """
    eigenvalues = list(eigenvalues)
    n = len(eigenvalues)

    for i in range(n):
        degenerate_count = sum(
            abs(eigenvalues[i] - eigenvalues[j]) < tolerance
            for j in range(n) if j != i
        )
        if degenerate_count == 0:
            return i

    raise RuntimeError("Could not identify a unique non-degenerate eigenvalue based on the provided tolerance.")

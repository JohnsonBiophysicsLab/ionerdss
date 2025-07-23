"""
angles.py

This module provides utility functions for angles related calculateds

Note: This file includes code adapted from the 'pointgroup' package,
originally authored by Abel Carreras (https://github.com/abelcarreras/pointgroup),
and is licensed under the MIT License (which is attached to the end of this docstring)

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

def angles_between_vector_and_vectors(reference_vec, targets, tol=1e-5):
    """
    Compute angles (in radians) between a reference vector and each row in a matrix.

    Parameters
    ----------
    reference_vec : array_like of shape (3,)
        The reference 3D vector.
    targets : ndarray of shape (N, 3)
        Array of target 3D vectors to compute angles with respect to.
    tol : float
        Threshold below which vector norms are treated as zero.

    Returns
    -------
    angles : ndarray of shape (N,)
        Array of angles (in radians) between reference vector and each target vector.
    """
    targets = np.asarray(targets)
    ref_norm = np.linalg.norm(reference_vec)
    target_norms = np.linalg.norm(targets, axis=1)

    dot_products = np.dot(targets, reference_vec)

    angles = []
    for dot, target_norm in zip(dot_products, target_norms):
        denom = target_norm * ref_norm
        if denom < tol:
            angles.append(0.0)
        else:
            cos_theta = np.clip(dot / denom, -1.0, 1.0)
            angles.append(np.arccos(cos_theta))
    return np.array(angles)

def angles_from_points(p1, p2, p3, complementary = True):
    """
    Calculates the angles at point p2 formed by vectors p1->p2 and p2->p3.

    Args:
        p1 (Coords): The first points.
        p2 (Coords): The vertex points where the angle is calculated.
        p3 (Coords): The third points.

    Returns:
        float: The angle in degrees.
    """
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)
    p3 = np.asarray(p3)

    v1 = p2 - p1
    v2 = p3 - p2
    theta = np.arccos(np.einsum("ij,ij->i", v1, v2) / (np.linalg.norm(v1, axis = 1) * np.linalg.norm(v2, axis = 1)))
    if complementary:
        return theta
    return np.pi - theta

def dihedrals_from_points(p1, p2, p3, p4, tol=1e-8):
    """
    Compute dihedral angles (in radians) from sets of 4 points in 3D.

    Vectorized version that accepts arrays of points.

    Parameters
    ----------
    p1, p2, p3, p4 : ndarray of shape (N, 3)
        Arrays of N points representing sequential atoms or coordinates.

    tol : float, optional
        Threshold below which vectors are considered degenerate.

    Returns
    -------
    dihedrals : ndarray of shape (N,)
        Array of dihedral angles in radians for each set of 4 points.
    """
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)
    p3 = np.asarray(p3)
    p4 = np.asarray(p4)

    b0 = p2 - p1
    b1 = p3 - p2
    b2 = p4 - p3

    # Normalize b1 to prevent numerical instability
    b1_norm = np.linalg.norm(b1, axis=1, keepdims=True)
    b1_unit = np.divide(b1, b1_norm, where=b1_norm > tol)

    # Orthogonal components
    v = b0 - (np.sum(b0 * b1_unit, axis=1, keepdims=True)) * b1_unit
    w = b2 - (np.sum(b2 * b1_unit, axis=1, keepdims=True)) * b1_unit

    x = np.sum(v * w, axis=1)
    y = np.sum(np.cross(b1_unit, v) * w, axis=1)

    return np.arctan2(y, x)

def absolute_error_to_angle(error, points, tol=1e-8):
    """
    Convert absolute coordinate errors to angular errors (in radians),
    assuming origin-centered radial vectors.

    Parameters
    ----------
    error : float
        Absolute positional error (e.g., in angstroms or nanometers).
    points : ndarray of shape (N, 3)
        Array of 3D coordinates representing points from origin.
    tol : float
        Minimum radius threshold to avoid divide-by-zero.

    Returns
    -------
    angle_errors : ndarray of shape (N,)
        Angular errors in radians for each point.
    """
    points = np.asarray(points)
    radii = np.linalg.norm(points, axis=1)
    clipped_radii = np.clip(radii, tol, None)
    return error / clipped_radii

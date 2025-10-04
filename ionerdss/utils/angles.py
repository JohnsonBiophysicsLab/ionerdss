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

from ionerdss.utils.coords import Coords

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

def angles_from_points(p1, p2, p3):
    """
    Compute the angle formed at p2 between the vectors p1->p2 and p3->p2.

    Supports both scalar (1×3) and batch (N×3) inputs.

    Returns:
        float or np.ndarray: Angle(s) in radians.
    """
    def to_array(x):
        if isinstance(x, Coords):
            return x.as_array()
        elif isinstance(x, list):
            return np.array([c.as_array() if isinstance(c, Coords) else np.asarray(c) for c in x], dtype=float)
        x = np.asarray(x)
        if x.ndim == 1:
            return x
        elif x.ndim == 2 and x.shape[1] == 3:
            return x
        elif x.ndim == 2 and x.shape[0] == 3:
            return x.T
        else:
            raise ValueError("Input must be shape (3,) or (N, 3)")

    v1 = to_array(p2) - to_array(p1)
    v2 = to_array(p3) - to_array(p2)

    v1 = np.atleast_2d(v1)
    v2 = np.atleast_2d(v2)

    dot = np.einsum("ij,ij->i", v1, v2)
    norm = np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1)

    cos_theta = np.clip(dot / norm, -1.0, 1.0)
    angles = np.arccos(cos_theta)

    return angles[0] if angles.shape[0] == 1 else angles

def dihedrals_from_points(p1, p2, p3, p4, tol=1e-8):
    """
    Compute signed dihedral (torsion) angles for sequences of four 3D points.

    This vectorized implementation accepts N quadruplets and returns N angles
    in radians. For each quadruplet (p1, p2, p3, p4), let:
        b0 = p2 - p1
        b1 = p3 - p2
        b2 = p4 - p3
        b1̂ = b1 / ||b1||
        v  = b0 - (b0·b1̂) b1̂       # component of b0 orthogonal to b1
        w  = b2 - (b2·b1̂) b1̂       # component of b2 orthogonal to b1

    The dihedral angle φ is computed as
        φ = atan2( (b1̂ × v) · w , v · w )
    and lies in the range (-π, π].

    Sign convention (direction):
      • The rotation axis is +b1 (from p2 toward p3).
      • φ is positive if the rotation that carries v into w is counter-clockwise
        when looking along +b1 (i.e., applying the right-hand rule with thumb
        pointing from p2 to p3). Equivalently, φ > 0 when ((b1̂ × v) · w) > 0.

    Parameters
    ----------
    p1, p2, p3, p4 : array_like, shape (N, 3)
        Arrays of N Cartesian points forming N consecutive quadruplets.
        All four inputs must have the same leading dimension N.
    tol : float, optional
        Threshold used to guard normalization of b1. If ||b1|| <= tol, the
        corresponding dihedral is treated as degenerate.

    Returns
    -------
    dihedrals : ndarray, shape (N,)
        Signed dihedral angles (radians) in (-π, π].

    Notes
    -----
    • Degenerate configurations (e.g., ||b1|| ≈ 0 or v ≈ 0 or w ≈ 0) yield
      x ≈ 0 and y ≈ 0 and thus return φ ≈ 0 with this implementation.
      Increase `tol` if you prefer to classify near-collinear cases as
      degenerate more aggressively.
    • This formulation is numerically stable and avoids explicit plane normals.

    Examples
    --------
    >>> import numpy as np
    >>> # Single quadruplet (N=1) — returns array([phi])
    >>> p1 = np.array([[0., 0., 0.]])
    >>> p2 = np.array([[1., 0., 0.]])
    >>> p3 = np.array([[1., 1., 0.]])
    >>> p4 = np.array([[1., 1., 1.]])
    >>> dihedrals_from_points(p1, p2, p3, p4)  # rotates about +x from xy-plane toward +z
    array([1.57079633])  # ≈ +π/2 (right-hand rule about +x)
    """
    p1 = np.asarray(p1)
    p2 = np.asarray(p2)
    p3 = np.asarray(p3)
    p4 = np.asarray(p4)

    # Force cast to 2D tensor
    def ensure_2d(*pts):
        out = []
        for x in pts:
            a = np.asarray(x)
            if a.ndim == 1:
                a = np.array(a[None, :])          # make it (1, 3)
            out.append(a)
        return out

    p1, p2, p3, p4 = ensure_2d(p1, p2, p3, p4)

    # Get relative vectors
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

    # Get the dihedral from angle between two norm vectors
    dihedrals = np.arctan(y, x)

    return dihedrals[0] if dihedrals.shape[0] == 1 else dihedrals

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

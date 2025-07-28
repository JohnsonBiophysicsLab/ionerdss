"""
rigid_transform.py

This module provides functions for computing and applying optimal rigid body transformations 
between two sets of 3D points using the **Procrustes alignment method**.

### Procrustes Transform:
The goal is to find the best-fit rigid transformation — consisting of a rotation matrix `R` and 
a translation vector `t` — that maps one point cloud (source) onto another (target), minimizing
the root-mean-square deviation (RMSD) between them.

### Mathematical Principle:
Given two sets of 3D points:
    - `source_points` (N×3)
    - `target_points` (N×3)

We:
1. Center both point sets to their respective centroids.
2. Compute the cross-covariance matrix.
3. Use **Singular Value Decomposition (SVD)** to extract the optimal rotation.
4. Compute the optimal translation based on the aligned centroids.
5. Ensure the rotation is proper (i.e., det(R) = +1, avoiding reflections).

This is useful in structural biology (e.g., protein fitting), robotics, computer vision, etc.
"""

import numpy as np

def rigid_transform_3d(source_points, target_points):
    """
    Compute the optimal rigid transformation (rotation R and translation t)
    that aligns two 3D point sets using Procrustes alignment.

    Parameters
    ----------
    source_points : ndarray of shape (N, 3)
        The source point cloud.
    target_points : ndarray of shape (N, 3)
        The target point cloud.

    Returns
    -------
    rotation_matrix : ndarray of shape (3, 3)
        Optimal rotation matrix that aligns source to target.
    translation_vector : ndarray of shape (3,)
        Translation vector from aligned source to target.
    """
    source_points = np.asarray(source_points)
    target_points = np.asarray(target_points)

    assert source_points.shape == target_points.shape, "Point sets must be the same shape"

    # Compute centroids
    source_centroid = source_points.mean(axis=0)
    target_centroid = target_points.mean(axis=0)

    # Center the point sets
    source_centered = source_points - source_centroid
    target_centered = target_points - target_centroid

    # Compute covariance matrix
    covariance_matrix = source_centered.T @ target_centered

    # Singular Value Decomposition
    U, _, Vt = np.linalg.svd(covariance_matrix)
    rotation_matrix = Vt.T @ U.T

    # Ensure proper rotation (det = +1)
    if np.linalg.det(rotation_matrix) < 0:
        Vt[-1, :] *= -1
        rotation_matrix = Vt.T @ U.T

    # Compute translation
    translation_vector = target_centroid - rotation_matrix @ source_centroid

    return rotation_matrix, translation_vector


def apply_rigid_transform(rotation_matrix, translation_vector, points):
    """
    Apply a rigid transformation (rotation + translation) to 3D points.

    Parameters
    ----------
    rotation_matrix : ndarray of shape (3, 3)
        Rotation matrix.
    translation_vector : ndarray of shape (3,)
        Translation vector.
    points : ndarray of shape (3,) or (N, 3)
        Single point or array of points to transform.

    Returns
    -------
    transformed_points : ndarray
        Transformed point(s), same shape as input.
    """
    points = np.asarray(points)
    return (rotation_matrix @ points.T).T + translation_vector

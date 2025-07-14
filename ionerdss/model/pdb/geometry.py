"""
geometry.py

Provides geometric transformation utilities used throughout the NERDSS coarse-graining
pipeline, including 3D rigid alignment of homologous chains, vector transformations,
steric clash detection, and angle measurements.

Key Features
------------
- Rigid-body alignment using singular value decomposition (SVD), ensuring optimal
  rotation and translation that minimizes RMSD between two point clouds.
- Transformation application for atomic or coarse-grained coordinates.
- Utility functions for angular analysis between unit vectors or planes.
- Steric clash detection between coarse-grained spheres.

Functions
---------
- rigid_transform_chains(chain1_coords, chain2_coords): Convenience wrapper for rigid_transform_3d.
- check_steric_clashes(pos1, pos2, r1, r2, buffer=0.0): Returns True if two spheres overlap.
- angle_between(v1, v2): Computes angle (in radians) between two 3D vectors.
- dihedral(p1, p2, p3, p4): Returns dihedral angle (in radians) defined by 4 points.
"""

import numpy as np
from ionerdss.math.rigid_transform import rigid_transform_3d

def rigid_transform_chains(chain1_coords, chain2_coords):
    """
    Convenience wrapper to align one chain to another using Cα or COM coordinates.

    Parameters
    ----------
    chain1_coords : ndarray (N, 3)
        Original coordinates to be aligned.
    chain2_coords : ndarray (N, 3)
        Target coordinates.

    Returns
    -------
    R : ndarray (3, 3)
        Rotation matrix.
    t : ndarray (3,)
        Translation vector.
    """
    return rigid_transform_3d(chain1_coords, chain2_coords)


def check_steric_clashes(pos1, pos2, r1, r2, buffer=0.0):
    """
    Returns True if two coarse-grained spheres overlap beyond allowed buffer.

    Parameters
    ----------
    pos1 : array-like of shape (3,)
        Center of first molecule.
    pos2 : array-like of shape (3,)
        Center of second molecule.
    r1 : float
        Radius of first molecule.
    r2 : float
        Radius of second molecule.
    buffer : float, optional
        Extra distance allowed without clash. Default is 0.0 nm.

    Returns
    -------
    bool
        True if spheres are clashing.
    """
    d = np.linalg.norm(np.asarray(pos1) - np.asarray(pos2))
    return d < (r1 + r2 - buffer)

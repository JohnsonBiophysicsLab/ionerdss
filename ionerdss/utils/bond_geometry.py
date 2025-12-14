"""
bond_geometry.py

Author: yying7@jh.edu
(method confirmed by Dr. Sam Foley)

This module provides geometry utilities for analyzing the spatial configuration
between two interacting protein chains. Specifically, it computes a set of five
geometric angles:

    - θ₁, θ₂: angles between the interaction vector and the inter-chain vector
    - φ₁, φ₂: signed angles between the local normal vector and the interaction plane
    - ω: torsion angle between interaction vectors around the inter-chain axis

These angles are useful in structural analysis, symmetry-aware coarse-graining,
and identifying rigid-body relationships in macromolecular assemblies.

Main Functions
--------------
- compute_bond_angles_and_length : Compute five orientation-defining angles with explicit normal points
- generate_normals_from_binding_sites : Generate local normal vectors orthogonal to interaction vectors
- compute_bond_angles_and_length_auto : Automatically construct normals and compute angles

All vectors are assumed to be in Cartesian coordinates, and all angles are returned in radians.
"""

import math
import numpy as np

from ionerdss.utils.vectors import get_magnitude, convert_to_unit

def compute_bond_angles_and_length(com1, com2,
                                   bind_site1, bind_site2,
                                   normal_point1, normal_point2,
                                   tol=1e-10):
    """
    Compute five orientation-defining angles between two protein chains.
    
    Inputs
    ------
    com1 : array-like
        Center of mass of chain 1
    com2 : array-like
        Center of mass of chain 2
    bind_site1 : array-like
        Interaction site on chain 1
    bind_site2 : array-like
        Interaction site on chain 2
    normal_point1 : array-like
        Point along the local normal vector from com1
    normal_point2 : array-like
        Point along the local normal vector from com2
    tol : float, default to 1e-10
        Tolerance such that values with differences lower than this value are considered the same

    Returns
    -------
    tuple
        theta1, theta2 : angle between interaction vector and inter-chain vector (in radians)
        phi1, phi2     : signed angle between normal and interaction plane (in radians)
        omega          : signed torsion angle between v1 and v2 around inter-chain vector (in radians)
        sigma_magnitude: distance between interaction sites
        normal_point1, normal_point2 :  passed-through inputs for traceability
    """

    # Interaction vectors from center of mass to interaction site
    v1 = [bind_site1[i] - com1[i] for i in range(3)]
    v2 = [bind_site2[i] - com2[i] for i in range(3)]
    sigma1 = [bind_site1[i] - bind_site2[i] for i in range(3)]
    sigma2 = [-x for x in sigma1]
    sigma_magnitude = get_magnitude(sigma1)

    # Theta: angle between interaction vector and inter-chain vector
    theta1 = math.acos(np.dot(v1, sigma1) / (get_magnitude(v1) * get_magnitude(sigma1)))
    theta2 = math.acos(np.dot(v2, sigma2) / (get_magnitude(v2) * get_magnitude(sigma2)))

    # Local normal vectors
    normal1 = convert_to_unit([normal_point1[i] - com1[i] for i in range(3)])
    normal2 = convert_to_unit([normal_point2[i] - com2[i] for i in range(3)])

    # Torsion-like vectors for phi angles
    torsion1_plane = convert_to_unit(np.cross(v1, sigma1))
    torsion1_normal = convert_to_unit(np.cross(v1, normal1))
    torsion2_plane = convert_to_unit(np.cross(v2, sigma2))
    torsion2_normal = convert_to_unit(np.cross(v2, normal2))

    phi1 = math.acos(np.dot(torsion1_plane, torsion1_normal))
    phi2 = math.acos(np.dot(torsion2_plane, torsion2_normal))

    # Sign correction for phi1
    v1_unit = convert_to_unit(v1)
    n1_proj = [normal1[i] - v1_unit[i] * np.dot(v1_unit, normal1) for i in range(3)]
    sigma1_proj = [sigma1[i] - v1_unit[i] * np.dot(v1_unit, sigma1) for i in range(3)]
    phi1_dir = convert_to_unit(np.cross(sigma1_proj, n1_proj))
    if abs(v1_unit[0] - phi1_dir[0]) < tol:
        phi1 = -phi1
    elif abs(v1_unit[0] + phi1_dir[0]) >= tol:
        print("Warning: Unable to determine phi1 sign.")

    # Sign correction for phi2
    v2_unit = convert_to_unit(v2)
    n2_proj = [normal2[i] - v2_unit[i] * np.dot(v2_unit, normal2) for i in range(3)]
    sigma2_proj = [sigma2[i] - v2_unit[i] * np.dot(v2_unit, sigma2) for i in range(3)]
    phi2_dir = convert_to_unit(np.cross(sigma2_proj, n2_proj))
    if abs(v2_unit[0] - phi2_dir[0]) < tol:
        phi2 = -phi2
    elif abs(v2_unit[0] + phi2_dir[0]) >= tol:
        print("Warning: Unable to determine phi2 sign.")

    # Omega: torsion angle between v1 and v2 around sigma1
    a1 = np.cross(sigma1, v1) / get_magnitude(np.cross(sigma1, v1))
    a2 = np.cross(sigma1, v2) / get_magnitude(np.cross(sigma1, v2))
    omega = math.acos(np.dot(a1, a2))

    # Omega: torsion angle between v1 and v2 around sigma1
    cross1 = np.cross(sigma1, v1)
    cross2 = np.cross(sigma1, v2)

    if get_magnitude(cross1) < tol or get_magnitude(cross2) < tol:
        omega = 0.0
    else:
        a1 = cross1 / get_magnitude(cross1)
        a2 = cross2 / get_magnitude(cross2)
        omega = math.acos(np.clip(np.dot(a1, a2), -1.0, 1.0))

        # Sign correction for omega
        sigma1_unit = convert_to_unit(sigma1)
        v1_proj = [v1[i] - sigma1_unit[i] * np.dot(sigma1_unit, v1) for i in range(3)]
        v2_proj = [v2[i] - sigma1_unit[i] * np.dot(sigma1_unit, v2) for i in range(3)]
        omega_dir = convert_to_unit(np.cross(v1_proj, v2_proj))

        if abs(sigma1_unit[0] - omega_dir[0]) < tol:
            omega = -omega
        elif abs(sigma1_unit[0] + omega_dir[0]) >= tol:
            print("Warning: Unable to determine omega sign.")

    return theta1, theta2, phi1, phi2, omega, sigma_magnitude, normal_point1, normal_point2

def generate_normals_from_binding_sites(com1, com2, bind_site1, bind_site2, reference_axis=None, tol=1e-6):
    """
    Generate normal vectors orthogonal to binding vectors to avoid gimbal lock.
    
    Parameters
    ----------
    com1, com2 : array-like
        Centers of mass for chains 1 and 2.
    bind_site1, bind_site2 : array-like
        Binding sites for chains 1 and 2.
    reference_axis : array-like, optional
        Axis to use for normal generation (default: [0, 0, 1]).
    tol : float
        Threshold to switch reference axis if nearly parallel.

    Returns
    -------
    normal_point1, normal_point2 : list of float
        Normal endpoints located at unit distance from COM1 and COM2.
    """
    def compute_normal(com, bind_site):
        v = np.array(bind_site) - np.array(com)
        v_unit = v / np.linalg.norm(v)

        # Default to Z axis
        if reference_axis is None:
            ref = np.array([0.0, 0.0, 1.0])
        else:
            ref = np.array(reference_axis)

        # Switch if nearly colinear
        if abs(np.dot(v_unit, ref)) > 1.0 - tol:
            ref = np.array([0.0, 1.0, 0.0])

        # Generate orthogonal normal vector
        normal = np.cross(v_unit, ref)
        normal = normal / np.linalg.norm(normal)

        # Return a point along the normal direction, unit distance from COM
        return list(np.array(com) + normal)

    normal_point1 = compute_normal(com1, bind_site1)
    normal_point2 = compute_normal(com2, bind_site2)

    return normal_point1, normal_point2

def compute_bond_angles_and_length_auto(com1, com2,
                                        bind_site1, bind_site2,
                                        tol=1e-10,
                                        reference_axis=None,
                                        normal_distance=1.0):
    """
    Automatically generate normal vectors and compute bond angles and distance.

    This function combines:
    - `generate_normals_from_binding_sites`: to construct consistent normal vectors
      orthogonal to interaction vectors (avoiding gimbal lock)
    - `compute_bond_angles_and_length`: to compute the five orientation-defining angles

    Parameters
    ----------
    com1 : array-like
        Center of mass of chain 1
    com2 : array-like
        Center of mass of chain 2
    bind_site1 : array-like
        Interaction site on chain 1
    bind_site2 : array-like
        Interaction site on chain 2
    tol : float, optional
        Tolerance used for determining sign ambiguity in angles (default = 1e-10)
    reference_axis : array-like, optional
        Axis to use when computing normal vectors (default = [0, 0, 1])
    normal_distance : float, optional
        Distance from COM to generate the normal point (default = 1.0)

    Returns
    -------
    tuple
        theta1, theta2 : angle between interaction vector and inter-chain vector (in radians)
        phi1, phi2     : signed angle between normal and interaction plane (in radians)
        omega          : signed torsion angle between v1 and v2 around inter-chain vector (in radians)
        sigma_magnitude: distance between interaction sites
        normal_point1, normal_point2: automatically assigned normals
    """
    # Generate normal points orthogonal to binding vectors
    normal_point1, normal_point2 = generate_normals_from_binding_sites(
        com1, com2, bind_site1, bind_site2,
        reference_axis=reference_axis,
        tol=tol
    )

    # Optionally scale normal vector to desired length
    def scale_normal(com, normal_point):
        direction = np.array(normal_point) - np.array(com)
        direction = direction / np.linalg.norm(direction)
        return list(np.array(com) + normal_distance * direction)

    normal_point1 = scale_normal(com1, normal_point1)
    normal_point2 = scale_normal(com2, normal_point2)

    # Compute angles and length
    return compute_bond_angles_and_length(
        com1=com1,
        com2=com2,
        bind_site1=bind_site1,
        bind_site2=bind_site2,
        normal_point1=normal_point1,
        normal_point2=normal_point2,
        tol=tol
    )

if __name__ == "__main__":
    # confirmed by Dr. Foley!
    print("slimer!")
    com2 = np.array([3.3587652772, -3.358754019, 5.00000079])
    bind_site2 = np.array([0.707106781, -0.707106781, 1.9999999732051])
    com1 = np.array([0.0, 0.0, 0.0])
    bind_site1 = np.array([0.0, 0.0, 2.0])
    normal_point2 = com2 + np.array([0,0,1])
    normal_point1 = com1 + np.array([1,1,0])
    
    print(compute_bond_angles_and_length(
        com1=com1,
        com2=com2,
        bind_site1=bind_site1,
        bind_site2=bind_site2,
        normal_point1=normal_point1,
        normal_point2=normal_point2
    ))
    
    
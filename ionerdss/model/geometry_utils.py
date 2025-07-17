"""
Geometry utility functions for PDB model processing.

This module contains pure functions for geometric calculations used in
coarse-graining and interface analysis.
"""

import math
import numpy as np


def calc_angle(P, Q, R):
    """
    Calculate the angle at point Q formed by vectors P->Q and R->Q.

    Args:
        P (Coords): The first point.
        Q (Coords): The vertex point where the angle is calculated.
        R (Coords): The third point.

    Returns:
        float: The angle in degrees.
    """
    v1 = [(Q - P).x, (Q - P).y, (Q - P).z]
    v2 = [(R - Q).x, (R - Q).y, (R - Q).z]
    theta = np.degrees(math.acos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
    return theta


def sig_are_similar(sig1, sig2, dist_threshold_intra, dist_threshold_inter, angle_threshold):
    """
    Compare two groups of interface interaction geometry signatures.

    Args:
        sig1 (dict): The first interface signature.
        sig2 (dict): The second interface signature.
        dist_threshold_intra (float): Distance threshold for intra-molecular comparisons.
        dist_threshold_inter (float): Distance threshold for inter-molecular comparisons.
        angle_threshold (float): Angle threshold for comparisons.

    Returns:
        bool: True if the signatures are similar within the given thresholds, False otherwise.
    """
    for key in ("dA", "dB"):
        if abs(sig1[key] - sig2[key]) > dist_threshold_intra:
            return False
    for key in ("dAB",):
        if abs(sig1[key] - sig2[key]) > dist_threshold_inter:
            return False
    for key in ("thetaA", "thetaB"):
        if abs(sig1[key] - sig2[key]) > angle_threshold:
            return False
    return True


def sig_difference(sig1, sig2):
    """
    Compute the sum of relative differences between two signatures.
    
    Args:
        sig1 (dict): First signature with keys: dA, dB, dAB, thetaA, thetaB
        sig2 (dict): Second signature with keys: dA, dB, dAB, thetaA, thetaB
    
    Returns:
        float: Total relative difference between signatures
    """
    total_diff = 0.0
    for key in ("dA", "dB", "dAB", "thetaA", "thetaB"):
        denom = abs(sig1[key]) if abs(sig1[key]) > 1e-6 else 1.0  # avoid divide-by-zero
        total_diff += abs(sig1[key] - sig2[key]) / denom
    return total_diff


def compute_diffusion_constants_nm_us(R_nm, T=298.0, eta=1e-3):
    """
    Compute translational and rotational diffusion constants.
    
    Args:
        R_nm (float): radius in nanometers
        T (float): temperature in Kelvin
        eta (float): viscosity in Pa·s (default: water)
    
    Returns:
        tuple: (D_t in nm^2/μs, D_r in rad^2/μs)
    """
    kB = 1.380649e-23  # J/K
    R_m = R_nm * 1e-9  # convert nm to meters

    # Diffusion constants
    D_t_m2_per_s = kB * T / (6 * np.pi * eta * R_m)
    D_r_rad2_per_s = kB * T / (8 * np.pi * eta * R_m**3)

    # Convert units
    D_t_nm2_per_us = D_t_m2_per_s * 1e12  # m²/s → nm²/μs
    D_r_rad2_per_us = D_r_rad2_per_s * 1e-6  # rad²/s → rad²/μs

    return D_t_nm2_per_us, D_r_rad2_per_us
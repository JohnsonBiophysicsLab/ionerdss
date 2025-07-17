"""
Angle calculation utilities for binding geometry.

This module contains functions for calculating binding angles between molecules,
which are critical for determining the geometry of molecular interactions.
"""

import math
import numpy as np


def unit(x: np.ndarray, eps=10**-6) -> np.ndarray:
    """
    Normalize a vector to unit length, handling numerical precision errors.

    Args:
        x (np.ndarray): Input vector.
        eps (float, optional): Small threshold for numerical stability.

    Returns:
        np.ndarray: Unit vector.
    """
    x_norm = np.linalg.norm(x)
    if abs(x_norm-1) < eps:
        return x
    elif x_norm < eps:
        return np.zeros(3)
    else:
        return x/x_norm


def _clip_cosine_value(x: float, eps=10**-6) -> float:
    """
    Ensure cosine values remain in the range [-1, 1] for numerical stability.

    Args:
        x (float): Input value.
        eps (float, optional): Small numerical threshold.

    Returns:
        float: Corrected value within [-1, 1].
    """
    if x < -1 and abs(x+1) < eps:
        return -1
    elif x > 1 and abs(x-1) < eps:
        return 1
    elif -1 <= x <= 1:
        return x
    else:
        raise ValueError(f'{x} is out of the range of sin/cos')


def calculate_phi(v: np.ndarray, n: np.ndarray, sigma: np.ndarray, eps=10**-6) -> float:
    """
    Compute the phi angle given three vectors.

    Args:
        v (np.ndarray): Direction vector.
        n (np.ndarray): Normal vector.
        sigma (np.ndarray): Sigma direction.
        eps (float, optional): Small numerical threshold.

    Returns:
        float: Computed phi angle.
    """
    # calculate phi
    t1 = unit(np.cross(v, sigma))
    t2 = unit(np.cross(v, n))
    phi = math.acos(_clip_cosine_value(np.dot(t1, t2)))

    # determine the sign of phi (+/-)
    v_uni = unit(v)
    n_proj = n - v_uni * np.dot(v_uni, n)
    sigma_proj = sigma - v_uni * np.dot(v_uni, sigma)
    phi_dir = unit(np.cross(sigma_proj, n_proj))

    if np.dot(v_uni, phi_dir) > 0:
        phi = -phi
    else:
        phi = phi
    
    return phi


def angles(com1, com2, int_site1, int_site2, normal_point1, normal_point2, eps=10**-6):
    """
    Compute binding angles for two molecules based on their center-of-mass (COM),
    interface sites, and normal vectors.

    This function determines five angles (theta1, theta2, phi1, phi2, omega) that describe
    the relative orientation between two interacting molecules.

    Args:
        com1 (array-like): Coordinates of the center-of-mass of molecule 1.
        com2 (array-like): Coordinates of the center-of-mass of molecule 2.
        int_site1 (array-like): Coordinates of the interface site on molecule 1.
        int_site2 (array-like): Coordinates of the interface site on molecule 2.
        normal_point1 (array-like): A point defining the orientation of molecule 1.
        normal_point2 (array-like): A point defining the orientation of molecule 2.
        eps (float, optional): Small numerical threshold to prevent division errors. Default is 1e-6.

    Returns:
        tuple:
            - str: Theta1 (binding angle in radians, formatted as string).
            - str: Theta2 (binding angle in radians, formatted as string).
            - str: Phi1 (torsion angle in radians, formatted as string, or 'nan' if undefined).
            - str: Phi2 (torsion angle in radians, formatted as string, or 'nan' if undefined).
            - str: Omega (twist angle in radians, formatted as string).
            - float: Distance between interface sites (sigma magnitude).
    """
    # Convert sequences into arrays for convenience
    com1 = np.array(com1)
    com2 = np.array(com2)
    int_site1 = np.array(int_site1)
    int_site2 = np.array(int_site2)
    normal_point1 = np.array(normal_point1)
    normal_point2 = np.array(normal_point2)

    # Get Vectors
    v1 = int_site1 - com1  # from COM to interface (particle 1)
    v2 = int_site2 - com2  # from COM to interface (particle 2)
    sigma1 = int_site1 - int_site2  # sigma, from p2 to p1
    sigma2 = int_site2 - int_site1  # sigma, from p1 to p2
    n1 = unit(normal_point1 - com1)  # normal vector for p1
    n2 = unit(normal_point2 - com2)  # normal vector for p2

    # Calculate the magnitude of sigma
    sigma_mag = np.linalg.norm(sigma1)

    # Calculate theta1 and theta2
    theta1 = math.acos(_clip_cosine_value(np.dot(unit(v1), unit(sigma1))))
    theta2 = math.acos(_clip_cosine_value(np.dot(unit(v2), unit(sigma2))))

    # Calculate phi1 and phi2
    phi1 = calculate_phi(v1, n1, sigma1, eps)
    phi2 = calculate_phi(v2, n2, sigma2, eps)

    # Calculate omega
    omega = math.acos(_clip_cosine_value(np.dot(unit(n1), unit(n2))))

    # Determine the sign of omega
    omega_dir = unit(np.cross(n1, n2))
    if np.dot(unit(sigma1), omega_dir) > 0:
        omega = -omega
    else:
        omega = omega

    # Format outputs
    theta1_str = f"{theta1:.6f}"
    theta2_str = f"{theta2:.6f}"
    phi1_str = f"{phi1:.6f}" if not math.isnan(phi1) else "nan"
    phi2_str = f"{phi2:.6f}" if not math.isnan(phi2) else "nan"
    omega_str = f"{omega:.6f}"

    return theta1_str, theta2_str, phi1_str, phi2_str, omega_str, sigma_mag
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
    v1 = int_site1 - com1 # from COM to interface (particle 1)
    v2 = int_site2 - com2  # from COM to interface (particle 2)
    sigma1 = int_site1 - int_site2 # sigma, from p2 to p1
    sigma2 = int_site2 - int_site1  # sigma, from p1 to p2
    n1 = unit(normal_point1 - com1) # normal vector for p1
    n2 = unit(normal_point2 - com2) # normal vector for p2

    # Calculate the magnitude of sigma
    sigma_magnitude = np.linalg.norm(sigma1)

    # Calculate theta1 and theta2
    costheta1 = np.dot(v1, sigma1) / np.linalg.norm(v1) / np.linalg.norm(sigma1)
    costheta2 = np.dot(v2, sigma2) / np.linalg.norm(v2) / np.linalg.norm(sigma2)
    theta1 = math.acos(_clip_cosine_value(costheta1))
    theta2 = math.acos(_clip_cosine_value(costheta2))

    # check geometry
    errormsg = ''
    iferror = False # determine if v // n
    if np.linalg.norm(np.cross(n1, v1)) < eps:
        iferror = True
        errormsg += f'\n\tn1 ({n1}) and v1 ({v1}) parallel, phi1 not available'
    if np.linalg.norm(np.cross(n2, v2)) < eps:
        iferror = True
        errormsg += f'\n\tn2 ({n2}) and v2 ({v2}) parallel, phi2 not available'
    if iferror:
        raise ValueError(errormsg)

    # determine if phi1 exists (v1 // sigma1 ?)
    if np.linalg.norm(np.cross(sigma1, v1)) < eps:
        phi1 = float('nan')
        # omega_parallel = True
        omega_t1 = unit(np.cross(sigma1, n1))
    else:
        phi1 = calculate_phi(v1, n1, sigma1)
        omega_t1 = unit(np.cross(sigma1, v1))

    # determine if phi2 exists (v2 // sigma2 ?)
    if np.linalg.norm(np.cross(sigma2, v2)) < eps:
        phi2 = float('nan')
        # omega_parallel = True
        omega_t2 = unit(np.cross(sigma1, n2))
    else:
        phi2 = calculate_phi(v2, n2, sigma2)
        omega_t2 = unit(np.cross(sigma1, v2))

    # calculate omega (both cases are same)
    omega = math.acos(_clip_cosine_value(np.dot(omega_t1, omega_t2)))
    # determine the sign of omega (+/-)
    sigma1_uni = unit(sigma1)
    sigma1xomega_t1 = np.cross(sigma1, omega_t1)
    sigma1xomega_t2 = np.cross(sigma1, omega_t2)
    omega_dir = unit(np.cross(sigma1xomega_t1, sigma1xomega_t2))
    if np.dot(sigma1_uni, omega_dir) > 0:
        omega = -omega
    else:
        omega = omega

    if abs(theta1 - np.pi) < eps:
        theta1 = 'M_PI'
    else:
        theta1 = "%.6f" % theta1
    if abs(theta2 - np.pi) < eps:
        theta2 = 'M_PI'
    else:
        theta2 = "%.6f" % theta2
    if abs(phi1 - np.pi) < eps:
        phi1 = 'M_PI'
    else:
        phi1 = "%.6f" % phi1
    if abs(phi2 - np.pi) < eps:
        phi2 = 'M_PI'
    else:
        phi2 = "%.6f" % phi2
    if abs(omega - np.pi) < eps:
        omega = 'M_PI'
    else:
        omega = "%.6f" % omega

    return theta1, theta2, phi1, phi2, omega, sigma_magnitude
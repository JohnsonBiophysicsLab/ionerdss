"""
geometry.py

Geometric utilities for spherical capsid construction and analysis of HIV Gag subunits.
Includes fitting routines, spherical coordinate transforms, and angular measurements.
"""

import numpy as np

def calculateAngles(c1, c2, p1, p2, n1=None, n2=None):
    """
    Compute orientation angles between two molecules and their interaction sites.

    Parameters
    ----------
    c1, c2 : np.ndarray
        Centers of mass of molecules 1 and 2.
    p1, p2 : np.ndarray
        Interaction site positions on molecules 1 and 2.
    n1, n2 : np.ndarray or None
        Optional normal vectors (defaults to [0, 0, 1]).

    Returns
    -------
    tuple of float
        (theta1, theta2, phi1, phi2, omega) in radians.
    """
    if n1 is None: n1 = np.array([0, 0, 1])
    if n2 is None: n2 = np.array([0, 0, 1])

    v1 = p1 - c1
    v2 = p2 - c2
    sigma = p1 - p2

    theta1 = np.arccos(np.dot(v1, sigma) / (np.linalg.norm(v1) * np.linalg.norm(sigma)))
    theta2 = np.arccos(np.dot(v2, -sigma) / (np.linalg.norm(v2) * np.linalg.norm(sigma)))

    # phi1
    t1 = np.cross(v1, sigma)
    t2 = np.cross(v1, n1)
    phi1 = np.arccos(np.dot(t1/np.linalg.norm(t1), t2/np.linalg.norm(t2)))
    if np.dot(np.cross(v1, t1), t2) > 0:
        phi1 = -phi1

    # phi2
    t1 = np.cross(v2, -sigma)
    t2 = np.cross(v2, n2)
    phi2 = np.arccos(np.dot(t1/np.linalg.norm(t1), t2/np.linalg.norm(t2)))
    if np.dot(np.cross(v2, t1), t2) > 0:
        phi2 = -phi2

    # omega
    if not np.isclose(np.linalg.norm(np.cross(v1, sigma)), 0) and not np.isclose(np.linalg.norm(np.cross(v2, -sigma)), 0):
        t1 = np.cross(sigma, v1)
        t2 = np.cross(sigma, v2)
    else:
        t1 = np.cross(sigma, n1)
        t2 = np.cross(sigma, n2)

    omega = np.arccos(np.dot(t1, t2) / (np.linalg.norm(t1) * np.linalg.norm(t2)))
    if np.dot(np.cross(sigma, t1), t2) > 0:
        omega = -omega

    return theta1, theta2, phi1, phi2, omega

def calculate_rmsd(centers, xyzR):
    """RMSD from sphere surface."""
    x0, y0, z0, r0 = xyzR
    diffs = np.linalg.norm(centers - np.array([x0, y0, z0]), axis=1) - r0
    return np.sum(diffs ** 2)

def calculate_gradient(centers, xyzR):
    """Gradient of RMSD sphere fitting objective."""
    x0, y0, z0, r0 = xyzR
    grad = np.zeros(4)
    for xi, yi, zi in centers:
        ri = np.linalg.norm([xi - x0, yi - y0, zi - z0])
        if ri == 0: continue
        diff = ri - r0
        grad[0] += (-2.0/ri) * diff * (xi - x0)
        grad[1] += (-2.0/ri) * diff * (yi - y0)
        grad[2] += (-2.0/ri) * diff * (zi - z0)
        grad[3] += -2.0 * diff
    return grad

def determine_gagTemplate_structure(numGag, positionsVec):
    """Compute template subunit internal frame and average internal coordinates."""
    coeffs = np.zeros((5, 3, numGag))
    basis = np.zeros((3, 3, numGag))

    for i in range(numGag):
        center = positionsVec[6*i]
        interfaces = positionsVec[6*i+1:6*i+6]
        vec1 = center / np.linalg.norm(center)
        vec2 = interfaces[0] - center
        vec3 = np.cross(vec1, vec2)
        vec3 /= np.linalg.norm(vec3)
        vec2 = np.cross(vec3, vec1)
        vec2 /= np.linalg.norm(vec2)

        basis[0,:,i] = vec1
        basis[1,:,i] = vec2
        basis[2,:,i] = vec3

        A = np.stack([vec1, vec2, vec3])
        for j in range(5):
            rel = interfaces[j] - center
            coeffs[j,:,i] = np.dot(rel, np.linalg.inv(A))

    mean_coeffs = np.mean(coeffs, axis=2)
    chosen_idx = 0
    gag_center = positionsVec[6*chosen_idx]
    b1, b2, b3 = basis[:, :, chosen_idx]

    template = np.zeros((6, 3))
    template[0] = gag_center
    for i in range(5):
        template[i+1] = (mean_coeffs[i,0]*b1 + mean_coeffs[i,1]*b2 + mean_coeffs[i,2]*b3 + gag_center)
    return template

def xyz_to_sphere_coordinates(pos):
    """Convert Cartesian coords to spherical (theta, phi, r)."""
    x, y, z = pos
    r = np.linalg.norm(pos)
    theta = np.arccos(z / r)
    phi = np.arccos(x / (r * np.sin(theta)))
    if y < 0:
        phi = 2*np.pi - phi
    return [theta, phi, r]

def translate_gags_on_sphere(hexmer, from_center, to_center):
    """Reproject a hexamer from one spherical location to another."""
    vec1 = from_center / np.linalg.norm(from_center)
    vec2 = to_center - from_center
    vec3 = np.cross(vec1, vec2)
    vec3 /= np.linalg.norm(vec3)
    vec2 = np.cross(vec3, vec1)
    vec2 /= np.linalg.norm(vec2)

    coeffs = np.array([np.dot(hexmer[i] - from_center, np.linalg.inv(np.stack([vec1, vec2, vec3])))
                       if np.linalg.norm(hexmer[i] - from_center) > 1e-10 else [0, 0, 0]
                       for i in range(hexmer.shape[0])])

    vec1 = to_center / np.linalg.norm(to_center)
    vec2 = to_center - from_center
    vec3 = np.cross(vec1, vec2)
    vec3 /= np.linalg.norm(vec3)
    vec2 = np.cross(vec3, vec1)
    vec2 /= np.linalg.norm(vec2)

    newhexmer = np.array([coeffs[i,0]*vec1 + coeffs[i,1]*vec2 + coeffs[i,2]*vec3 + to_center
                          for i in range(hexmer.shape[0])])
    return newhexmer

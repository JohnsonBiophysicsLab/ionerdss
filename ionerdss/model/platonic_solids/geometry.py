"""Geometry utilities for Platonic solids generation."""

import math
import numpy as np
from typing import Tuple

def angle_cal(COM1: np.ndarray, leg1: np.ndarray, COM2: np.ndarray, leg2: np.ndarray) -> Tuple[float, float, float, float, float]:
    """Calculates angles between vectors based on given inputs.

    Args:
        COM1: Center of Mass (COM) for the first leg.
        leg1: Endpoint of the first leg.
        COM2: Center of Mass (COM) for the second leg.
        leg2: Endpoint of the second leg.

    Returns:
        tuple: (theta1, theta2, phi1, phi2, omega) in radians rounded to 8 decimal places.
    """
    n = 8
    # Ensure inputs are numpy arrays
    c1 = np.array(COM1)
    p1 = np.array(leg1)
    c2 = np.array(COM2)
    p2 = np.array(leg2)
    
    v1 = p1 - c1
    v2 = p2 - c2
    sig1 = p1 - p2
    sig2 = -sig1
    
    # Calculate angles
    # Note: Added error handling for floating point precision issues in arccos
    
    def safe_acos(x):
        return math.acos(max(-1.0, min(1.0, x)))

    dot_theta1 = np.dot(v1, sig1) / (np.linalg.norm(v1) * np.linalg.norm(sig1))
    theta1 = round(safe_acos(dot_theta1), n)
    
    dot_theta2 = np.dot(v2, sig2) / (np.linalg.norm(v2) * np.linalg.norm(sig2))
    theta2 = round(safe_acos(dot_theta2), n)
    
    t1 = np.cross(v1, sig1)
    t2 = np.cross(v1, c1) 
    t1_hat = t1 / np.linalg.norm(t1)
    t2_hat = t2 / np.linalg.norm(t2)
    phi1 = round(safe_acos(np.around(np.dot(t1_hat, t2_hat), n)), n)
    
    t3 = np.cross(v2, sig2)
    t4 = np.cross(v2, c2)
    t3_hat = t3 / np.linalg.norm(t3)
    t4_hat = t4 / np.linalg.norm(t4)
    phi2 = round(safe_acos(np.around(np.dot(t3_hat, t4_hat), n)), n)
    
    t1_ = np.cross(sig1, v1)
    t2_ = np.cross(sig1, v2)
    t1__hat = t1_ / np.linalg.norm(t1_)
    t2__hat = t2_ / np.linalg.norm(t2_)
    omega = round(safe_acos(np.around(np.dot(t1__hat, t2__hat), n)), n)
    
    return theta1, theta2, phi1, phi2, omega

def distance(a: np.ndarray, b: np.ndarray) -> float:
    """Compute Euclidean distance between two points."""
    return float(np.linalg.norm(np.array(a) - np.array(b)))

def mid_pt(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Compute mid-point between two points."""
    return (np.array(a) + np.array(b)) / 2.0

"""
steric_clash.py

This file contains functions for steric clash checking
"""

from scipy.spatial import KDTree

def check_clashes_between_two_sets(points_1, points_2,
                                   cutoff: float = 3.5,
                                   number_threshold: int = 2):
    """
    Detects steric clashes between two sets of molecular points.

    Args:
        points_1 (np.ndarray): N x 3 coordinates for the first molecule.
        points_2 (np.ndarray): M x 3 coordinates for the second molecule.
        cutoff (float, optional): Distance threshold (default: 3.5 Å).
        threshold (int, optional): Minimum number of close contacts to flag a clash (default: 2).

    Returns:
        bool: True if a steric clash is detected, False otherwise.
    """
    tree = KDTree(points_2)
    clashes = tree.query_ball_point(points_1, r=cutoff)
    return any(len(clash) >= number_threshold for clash in clashes)

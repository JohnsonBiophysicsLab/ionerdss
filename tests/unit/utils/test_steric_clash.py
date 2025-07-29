"""
Unit tests for `check_clashes_between_two_sets` in `steric_clash.py`.

This function determines whether two sets of 3D molecular coordinates are in steric clash
based on a spatial distance cutoff and a minimum number of nearby neighbors.
"""

import unittest
import numpy as np
from scipy.spatial import KDTree
# Replace this import with the actual path if needed
# from your_module.steric_clash import check_clashes_between_two_sets

from ionerdss.utils.steric_clash import check_clashes_between_two_sets


class TestStericClash(unittest.TestCase):

    def test_no_clash(self):
        """Points far apart should not clash."""
        points_1 = np.array([[0.0, 0.0, 0.0],
                             [5.0, 5.0, 5.0]])
        points_2 = np.array([[10.0, 10.0, 10.0],
                             [15.0, 15.0, 15.0]])
        result = check_clashes_between_two_sets(points_1, points_2)
        self.assertFalse(result)

    def test_clash_detected(self):
        """Close points should result in a clash when threshold is met."""
        points_1 = np.array([[0.0, 0.0, 0.0],
                             [1.0, 1.0, 1.0]])
        points_2 = np.array([[0.5, 0.0, 0.0],
                             [1.1, 1.0, 1.0]])
        result = check_clashes_between_two_sets(points_1, points_2, cutoff=2.0, number_threshold=2)
        self.assertTrue(result)

    def test_below_threshold(self):
        """Should not clash if number of neighbors is below threshold."""
        points_1 = np.array([[0.0, 0.0, 0.0]])
        points_2 = np.array([[0.5, 0.0, 0.0],
                             [5.0, 5.0, 5.0]])
        result = check_clashes_between_two_sets(points_1, points_2, cutoff=2.0, number_threshold=2)
        self.assertFalse(result)

    def test_exact_threshold(self):
        """Should clash when the number of close neighbors equals the threshold."""
        points_1 = np.array([[0.0, 0.0, 0.0]])
        points_2 = np.array([[0.5, 0.0, 0.0],
                             [0.0, 0.5, 0.0]])
        result = check_clashes_between_two_sets(points_1, points_2, cutoff=1.0, number_threshold=2)
        self.assertTrue(result)

    def test_zero_points(self):
        """Should not crash on empty input."""
        points_1 = np.zeros((0, 3))
        points_2 = np.zeros((0, 3))
        result = check_clashes_between_two_sets(points_1, points_2)
        self.assertFalse(result)

if __name__ == "__main__":
    unittest.main()

import unittest
import numpy as np
from ionerdss.math.angles import (
    absolute_error_to_angle,
    angles_between_vector_and_vectors,
    angles_from_points,
    
    )

class TestAngles(unittest.TestCase):

    def test_absolute_error_to_angle_basic(self):
        points = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 4]])
        error = 0.1
        expected = np.array([0.1 / 1, 0.1 / 2, 0.1 / 4])
        result = absolute_error_to_angle(error, points)
        np.testing.assert_allclose(result, expected, rtol=1e-6)

    def test_absolute_error_to_angle_with_small_radius(self):
        points = np.array([[0, 0, 0], [1e-12, 0, 0]])
        error = 0.1
        result = absolute_error_to_angle(error, points, tol=1e-8)
        self.assertTrue(np.all(result >= 0))
        self.assertAlmostEqual(result[0], 0.1 / 1e-8)
        self.assertAlmostEqual(result[1], 0.1 / 1e-8)

    def test_angles_between_vector_and_vectors_basic(self):
        ref = np.array([1, 0, 0])
        angle_deg = 47.3
        angle_rad = np.deg2rad(angle_deg)
        target_47_3 = np.array([np.cos(angle_rad), np.sin(angle_rad), 0])

        targets = np.array([
            [1, 0, 0],       # 0 degrees
            [0, 1, 0],       # 90 degrees
            [-1, 0, 0],      # 180 degrees
            [1, 1, 0],       # ~45 degrees
            target_47_3,     # 47.3 degrees
        ])

        result = angles_between_vector_and_vectors(ref, targets)
        expected = np.array([
            0.0,
            np.pi / 2,
            np.pi,
            np.pi / 4,
            angle_rad,
        ])
        np.testing.assert_allclose(result, expected, rtol=1e-6)

    def test_angles_between_vector_and_zero_vector(self):
        ref = np.array([1, 0, 0])
        targets = np.array([
            [0, 0, 0],  # zero norm
            [0, 0, 0],  # zero norm
        ])
        result = angles_between_vector_and_vectors(ref, targets)
        expected = np.array([0.0, 0.0])
        np.testing.assert_array_equal(result, expected)

    def test_angles_between_vector_and_vectors_with_tolerance(self):
        ref = np.array([1e-12, 0, 0])
        targets = np.array([[1, 0, 0]])
        result = angles_between_vector_and_vectors(ref, targets, tol=1e-8)
        expected = np.array([0.0])  # Denominator too small, treated as 0
        np.testing.assert_array_equal(result, expected)

    def test_angles_from_points(self):

        # initialize list
        p1 = []
        p2 = []
        p3 = []

        # 90 degrees at vertex (1, 0, 0]
        p1.append([1, 1, 0])
        p2.append([1, 0, 0])
        p3.append([2, 0, 0])

        # 180 degrees at center
        p1.append([0, 0, 0])
        p2.append([1, 0, 0])
        p3.append([2, 0, 0])

        # 60 degrees at any vertex of equilateral triangle
        p1.append([0, 0, 0])
        p2.append([1, 0, 0])
        p3.append([0.5, np.sqrt(3)/2, 0])

        # colinear inward (0 deg)
        p1.append([2, 0, 0])
        p2.append([0, 0, 0])
        p3.append([1, 0, 0])

        # tetrahedron vertices
        # Approx 70.53° angle between two bonds in tetrahedral geometry
        p1.append([1, 1, 1])
        p2.append([0, 0, 0])
        p3.append([1, -1, -1])

        # test numerical stability
        p1.append([-1, 0, 0])
        p2.append([0, 0, 0])
        p3.append([1, 1e-7, 0])

        # test
        angle = angles_from_points(p1, p2, p3)
        np.testing.assert_allclose(angle, np.array([np.pi/2.0,
                                                0.0,
                                                2 * np.pi/3.0,
                                                np.pi,
                                                np.pi-1.910633236249,
                                                9.884312e-08]))

if __name__ == "__main__":
    unittest.main()

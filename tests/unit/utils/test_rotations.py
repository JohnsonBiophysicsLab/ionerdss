"""
test_rotations.py

Unit tests for `math/rotations.py`, which defines tools for working with 3D
rotational symmetry using axis-angle representations.

Tests cover:
- Construction of proper Cₙ rotation matrices
- Construction of improper Sₙ operations (rotation + reflection)
- Orthogonality and determinant of generated matrices
- Error handling for invalid input axes
"""

import unittest
import numpy as np
from math import isclose
from ionerdss.utils.rotations import (
    rotation_matrix,
    Rotation,
    ImproperRotation
)

class TestRotationMatrix(unittest.TestCase):

    def test_rotation_identity(self):
        axis = [1, 0, 0]
        mat = rotation_matrix(axis, 0.0)
        np.testing.assert_allclose(mat, np.eye(3), atol=1e-8)

    def test_rotation_known_angle(self):
        axis = [0, 0, 1]
        angle = np.pi / 2  # 90 degrees
        vec = np.array([1, 0, 0])
        mat = rotation_matrix(axis, angle)
        rotated = mat @ vec
        expected = np.array([0, 1, 0])
        np.testing.assert_allclose(rotated, expected, atol=1e-8)

    def test_invalid_axis(self):
        with self.assertRaises(ValueError):
            rotation_matrix([0, 0, 0], np.pi)

class TestCnRotation(unittest.TestCase):

    def test_order_4_rotation(self):
        rot = Rotation([0, 0, 1], order=4)
        mat = rot.get_matrix()

        # 90 degree rotation about Z axis
        vec = np.array([1, 0, 0])
        rotated = mat @ vec
        expected = np.array([0, 1, 0])
        np.testing.assert_allclose(rotated, expected, atol=1e-8)

    def test_rotation_matrix_properties(self):
        mat = Rotation([1, 1, 1], order=3).get_matrix()
        should_be_identity = mat @ mat @ mat
        np.testing.assert_allclose(should_be_identity, np.eye(3), atol=1e-8)

        # Orthogonality
        np.testing.assert_allclose(mat.T @ mat, np.eye(3), atol=1e-8)
        # Determinant = 1 (proper rotation)
        self.assertTrue(isclose(np.linalg.det(mat), 1.0, abs_tol=1e-8))

class TestSnImproperRotation(unittest.TestCase):

    def test_improper_rotation_reflects(self):
        imp = ImproperRotation([0, 0, 1], order=2)
        mat = imp.get_matrix()

        # 180-degree rotation + reflection across XY plane (Z -> -Z)
        vec = np.array([1, 0, 1])
        result = mat @ vec
        expected = np.array([-1, 0, -1])
        np.testing.assert_allclose(result, expected, atol=1e-8)

    def test_improper_matrix_properties(self):
        mat = ImproperRotation([1, 0, 0], order=2).get_matrix()
        # Should still be orthogonal
        np.testing.assert_allclose(mat.T @ mat, np.eye(3), atol=1e-8)
        # Improper => det = -1
        self.assertTrue(isclose(np.linalg.det(mat), -1.0, abs_tol=1e-8))


if __name__ == "__main__":
    unittest.main()

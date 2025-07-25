"""
Unit tests for rigid transformation utilities: `rigid_transform_3d` and `apply_transform`.

These tests verify:
- That the optimal rotation and translation returned by `rigid_transform_3d`
  align two point clouds correctly using Procrustes alignment.
- That `apply_transform` correctly applies a rigid body transform to a single point or array of points.
- That the rotation matrix is proper (i.e., no reflection), with determinant +1.
"""

import unittest
import numpy as np
from numpy.testing import assert_allclose


from ionerdss.math.rigid_transform import (
    rigid_transform_3d,
    apply_transform,
    )

class TestRigidTransform3D(unittest.TestCase):

    def test_identity_transform(self):
        """Rigid transform between identical point clouds should return identity R and zero t."""
        P = np.random.rand(10, 3)
        Q = P.copy()
        R, t = rigid_transform_3d(P, Q)

        assert_allclose(R, np.eye(3), atol=1e-6)
        assert_allclose(t, np.zeros(3), atol=1e-6)

    def test_translation_only(self):
        """Should correctly recover a pure translation."""
        P = np.random.rand(5, 3)
        t_true = np.array([1.0, -2.0, 3.0])
        Q = P + t_true

        R, t = rigid_transform_3d(P, Q)

        assert_allclose(R, np.eye(3), atol=1e-6)
        assert_allclose(t, t_true, atol=1e-6)

    def test_rotation_and_translation(self):
        """Should recover known rotation and translation."""
        angle = np.pi / 4  # 45 degrees around z-axis
        R_true = np.array([
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle),  np.cos(angle), 0],
            [0,              0,             1]
        ])
        t_true = np.array([1.0, 2.0, -1.0])

        P = np.random.rand(6, 3)
        Q = apply_transform(R_true, t_true, P)

        R, t = rigid_transform_3d(P, Q)

        assert_allclose(R, R_true, atol=1e-6)
        assert_allclose(t, t_true, atol=1e-6)

    def test_apply_transform_single_point(self):
        """apply_transform should work for single 3D point."""
        R = np.eye(3)
        t = np.array([1.0, 2.0, 3.0])
        x = np.array([0.5, 0.5, 0.5])
        y = apply_transform(R, t, x)
        expected = x + t
        assert_allclose(y, expected, atol=1e-6)

    def test_apply_transform_multiple_points(self):
        """apply_transform should work for multiple 3D points."""
        R = np.eye(3)
        t = np.array([1.0, -1.0, 0.0])
        X = np.array([[0.0, 0.0, 0.0],
                      [1.0, 2.0, 3.0]])
        Y = apply_transform(R, t, X)
        expected = X + t
        assert_allclose(Y, expected, atol=1e-6)

    def test_no_reflection_in_transform(self):
        """The computed rotation matrix should be proper (no reflection)."""
        P = np.random.rand(7, 3)
        Q = P + np.array([2.0, -1.0, 0.5])
        R, _ = rigid_transform_3d(P, Q)
        det = np.linalg.det(R)
        self.assertAlmostEqual(det, 1.0, places=6)

if __name__ == "__main__":
    unittest.main()

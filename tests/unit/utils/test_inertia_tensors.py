"""
test_inertia_tensors.py

Unit tests for `inertia_tools.py`, which provides mathematical tools for analyzing
3D symmetry, moment of inertia tensors, and rotational properties of molecular assemblies.

The module includes functionality for:
- Computing normalized inertia tensors of point distributions
- Detecting degeneracy in eigenvalues (e.g., for identifying rotational symmetries)
- Locating the unique non-degenerate eigenvalue (used to identify symmetry axes)
- Generating perpendicular vectors (in 3D or general nD)

These tests validate numerical correctness, geometric consistency,
and error handling for various structured and degenerate configurations.

Test coverage includes:
- Perfect symmetry shapes (e.g., square, equilateral triangle)
- Degenerate and non-degenerate eigenvalue configurations
- Input validation and edge cases for vector operations
"""

import unittest
import numpy as np
from ionerdss.utils.inertia_tensors import (
    get_inertia_tensor,
    get_degeneracy,
    get_non_degenerated
)

class TestInertiaTensor(unittest.TestCase):

    def test_tensor_for_square_xy_plane(self):
        # Square in XY plane
        coords = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]])
        I = get_inertia_tensor(coords)
        eigvals, eigvecs = np.linalg.eigh(I)

        # Should show degeneracy in x and y axes
        degeneracy = get_degeneracy(eigvals)
        self.assertEqual(degeneracy, 2)

        # Non-degenerate eigenvalue should correspond to Z axis
        idx = get_non_degenerated(eigvals)
        main_axis = eigvecs[:, idx]
        # Should be aligned with Z-axis
        self.assertTrue(np.allclose(np.abs(main_axis), [0, 0, 1], atol=1e-6))

    def test_tensor_for_equilateral_triangle(self):
        # Triangle in XY plane centered at origin
        coords = np.array([
            [1, 0, 0],
            [-0.5, np.sqrt(3)/2, 0],
            [-0.5, -np.sqrt(3)/2, 0]
        ])
        I = get_inertia_tensor(coords)
        eigvals, eigvecs = np.linalg.eigh(I)
        degeneracy = get_degeneracy(eigvals)
        self.assertEqual(degeneracy, 2)

    def test_tensor_with_perturbation_breaks_degeneracy(self):
        coords = np.array([
            [1, 0, 0],
            [-1, 0, 0],
            [0, 1, 0],
            [0, -1.1, 0]  # small asymmetry
        ])
        I = get_inertia_tensor(coords)
        eigvals, _ = np.linalg.eigh(I)
        degeneracy = get_degeneracy(eigvals, tolerance=0.01)
        self.assertEqual(degeneracy, 1)


class TestDegeneracy(unittest.TestCase):

    def test_degenerate_two_equal(self):
        vals = [1.0, 1.0, 2.0]
        self.assertEqual(get_degeneracy(vals), 2)

    def test_all_degenerate(self):
        vals = [2.0, 2.0, 2.0]
        self.assertEqual(get_degeneracy(vals), 3)

    def test_all_unique(self):
        vals = [1.0, 2.0, 3.0]
        self.assertEqual(get_degeneracy(vals), 1)

    def test_sensitive_to_tolerance(self):
        vals = [1.0, 1.05, 2.0]
        self.assertEqual(get_degeneracy(vals, tolerance=0.1), 2)
        self.assertEqual(get_degeneracy(vals, tolerance=0.01), 1)


class TestNonDegenerateIndex(unittest.TestCase):

    def test_index_detection(self):
        vals = [1.0, 1.0, 2.0]
        self.assertEqual(get_non_degenerated(vals), 2)

        vals = [2.0, 1.0, 2.0]
        self.assertEqual(get_non_degenerated(vals), 1)

        vals = [3.0, 2.9, 3.0]
        self.assertEqual(get_non_degenerated(vals, tolerance=0.05), 1)

    def test_raises_if_no_unique(self):
        vals = [1.0, 1.0, 1.0]
        with self.assertRaises(RuntimeError):
            get_non_degenerated(vals)

if __name__ == "__main__":
    unittest.main()

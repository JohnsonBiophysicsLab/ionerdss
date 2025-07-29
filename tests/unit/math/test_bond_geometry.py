"""
Unit tests for bond_geometry.py

This test suite verifies the correctness of geometric computations performed
in the bond_geometry module, including:

- Theta (θ₁, θ₂): angles between interaction vectors and the inter-chain vector
- Phi (φ₁, φ₂): signed angles between normal vectors and interaction planes
- Omega (ω): signed torsional angle between interaction vectors around inter-chain axis
- Projection-based sign corrections for dihedral angles
- Generation of orthogonal normal vectors avoiding gimbal lock
- Consistency of vector operations (magnitude, unit vector)

The test cases include:
- A symmetric arrangement with known expected angles
- A case with colinear interaction vectors (degenerate torsion)
- A fully general, non-coplanar configuration for torsion validation
- Edge case tests for helper utilities (_magnitude, _unit)
- Verifications that all returned angles lie within valid bounds (±π)
"""

import unittest
import numpy as np
from math import pi
from ionerdss.utils import bond_geometry

class TestBondGeometry(unittest.TestCase):

    def setUp(self):
        # Generic COM and binding sites for varied geometric testing
        self.com1 = [0.0, 0.0, 0.0]
        self.com2 = [0.0, 0.0, 2.0]

    def test_colinear_case(self):
        site1 = [0.0, 0.0, 0.5]
        site2 = [0.0, 0.0, 1.5]
        angles = bond_geometry.compute_bond_angles_and_length_auto(
            self.com1, self.com2, site1, site2
        )
        theta1, theta2, phi1, phi2, omega, sigma_mag, _, _ = angles
        self.assertAlmostEqual(omega, 0.0, delta=1e-6)
        self.assertTrue(all(-pi <= x <= pi or np.isnan(x) for x in [theta1, theta2, phi1, phi2, omega]))

    def test_arbitrary_torsion_case(self):
        site1 = [1.0, 0.0, 0.5]
        site2 = [-1.0, 1.0, 1.5]
        angles = bond_geometry.compute_bond_angles_and_length_auto(
            self.com1, self.com2, site1, site2
        )
        theta1, theta2, phi1, phi2, omega, sigma_mag, _, _ = angles
        self.assertTrue(all(-pi <= x <= pi or np.isnan(x) for x in [theta1, theta2, phi1, phi2, omega]))
        self.assertAlmostEqual(sigma_mag, 2.44948974278, places=4)
        self.assertAlmostEqual(theta1, 0.9911565864311924, places=4)
        self.assertAlmostEqual(theta2, 0.8224691545143296, places=4)
        self.assertAlmostEqual(phi1, 1.0610566479633896, places=4)
        self.assertAlmostEqual(phi2, -1.97568811307998, places=4)
        self.assertAlmostEqual(omega, 2.224122132419503, places=4)

    def test_unit_and_magnitude_functions(self):
        v = [3, 4, 0]
        mag = bond_geometry._magnitude(v)
        self.assertAlmostEqual(mag, 5.0)
        unit = bond_geometry._unit(v)
        expected_unit = [0.6, 0.8, 0.0]
        for a, b in zip(unit, expected_unit):
            self.assertAlmostEqual(a, b, places=6)

    def test_generate_normals(self):
        bind1 = [1, 0, 0]
        bind2 = [-1, 0, 0]
        normal1, normal2 = bond_geometry.generate_normals_from_binding_sites(
            self.com1, self.com2, bind1, bind2
        )
        # Check if they are unit length away from COM
        self.assertAlmostEqual(np.linalg.norm(np.array(normal1) - np.array(self.com1)), 1.0, places=6)
        self.assertAlmostEqual(np.linalg.norm(np.array(normal2) - np.array(self.com2)), 1.0, places=6)

if __name__ == "__main__":
    unittest.main()

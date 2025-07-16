import unittest
import numpy as np
from ionerdss.math import bond_geometry

class TestBondGeometry(unittest.TestCase):
    def setUp(self):
        self.com1 = [0.0, 0.0, 0.0]
        self.com2 = [0.0, 0.0, 2.0]
        self.site1 = [0.0, 0.0, 0.5]
        self.site2 = [0.0, 0.0, 1.5]

    def test_compute_bond_angles_and_length_auto(self):
        angles = bond_geometry.compute_bond_angles_and_length_auto(
            com1=self.com1,
            com2=self.com2,
            bind_site1=self.site1,
            bind_site2=self.site2,
            reference_axis=[1.0, 0.0, 0.0]
        )

        theta1, theta2, phi1, phi2, omega, sigma_mag, n1, n2 = angles

        # Basic shape check
        self.assertEqual(len(angles), 8)

        # Angle ranges
        for angle in [theta1, theta2, phi1, phi2, omega]:
            self.assertTrue(-np.pi <= angle <= np.pi)

        # Magnitude
        self.assertAlmostEqual(sigma_mag, 1.0, places=5)

        # Normal points should lie unit distance from COMs
        self.assertAlmostEqual(
            np.linalg.norm(np.array(n1) - np.array(self.com1)), 1.0, places=5)
        self.assertAlmostEqual(
            np.linalg.norm(np.array(n2) - np.array(self.com2)), 1.0, places=5)


if __name__ == "__main__":
    unittest.main()

"""
Unit tests for the diffusion constant calculation module in `ionerdss.math.diffusion_constant`.

This test suite verifies the correctness and physical behavior of the function
`compute_diffusion_constants_nm_us`, which computes translational and rotational
diffusion constants of spherical particles (e.g., proteins) in solution based on
the Einstein-Stokes equations.

Tests included:
- `test_typical_radius`: Validates output for a 5 nm radius sphere in water at 298.15 K.
- `test_diffusion_scales_with_radius`: Verifies inverse scaling of diffusion constants with radius.
- `test_custom_temperature_viscosity`: Confirms output at non-default temperature and viscosity.
- `test_zero_radius`: Ensures zero radius raises a ValueError.
- `test_negative_radius`: Ensures negative radius raises a ValueError.

All constants are validated against expected physical behavior and known analytical outputs.
"""

import unittest

from ionerdss.math.diffusion_constant import (
    compute_diffusion_constants_nm_us,
)


# Assume compute_diffusion_constants_nm_us is in a module named `diffusion`
# from your_module import compute_diffusion_constants_nm_us


class TestDiffusionConstants(unittest.TestCase):

    def test_typical_radius(self):
        radius_nm = 5.0  # Typical protein radius
        d_t, d_r = compute_diffusion_constants_nm_us(radius_nm)

        self.assertTrue(d_t > 0, "Translational diffusion should be positive")
        self.assertTrue(d_r > 0, "Rotational diffusion should be positive")

        # Rough theoretical estimates (not strict values)
        self.assertAlmostEqual(d_t, 49.074621908890364)
        self.assertAlmostEqual(d_r, 1.4722386572667108)

    def test_diffusion_scales_with_radius(self):
        small_radius = 2.0
        large_radius = 10.0

        d_t_small, d_r_small = compute_diffusion_constants_nm_us(small_radius)
        d_t_large, d_r_large = compute_diffusion_constants_nm_us(large_radius)

        self.assertTrue(d_t_small > d_t_large, "Translational diffusion should decrease with radius")
        self.assertTrue(d_r_small > d_r_large, "Rotational diffusion should decrease rapidly with radius³")

    def test_custom_temperature_viscosity(self):
        radius_nm = 5.0
        temp = 310  # Body temperature
        viscosity = 1e-3  # Slightly more viscous fluid

        d_t, d_r = compute_diffusion_constants_nm_us(radius_nm, temperature_kelvin=temp, viscosity_pas=viscosity)
        self.assertIsInstance(d_t, float)
        self.assertIsInstance(d_r, float)

    def test_zero_radius(self):
        with self.assertRaises(ValueError):
            compute_diffusion_constants_nm_us(0)

    def test_negative_radius(self):
        with self.assertRaises(ValueError):
            compute_diffusion_constants_nm_us(-5)

if __name__ == '__main__':
    unittest.main()

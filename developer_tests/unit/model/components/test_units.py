"""
Unit tests for ionerdss.model.components.units module.

Tests the Units class including initialization, serialization, and
unit consistency across the simulation framework.
"""

import unittest

# Import the module under test
from ionerdss.model.components.units import Units


class TestUnits(unittest.TestCase):
    """Test the Units class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create default units instance
        self.default_units = Units()

        # Create custom units instance
        self.custom_units = Units(
            coords="Å",
            binding_radius="nm",
            D_trans="nm^2/s",
            D_rot="deg^2/s"
        )

    def test_init_default(self):
        """Test initialization with default values."""
        units = Units()

        self.assertEqual(units.coords, "nm")
        self.assertEqual(units.binding_radius, "nm")
        self.assertEqual(units.D_trans, "nm^2/µs")
        self.assertEqual(units.D_rot, "rad^2/µs")

    def test_init_custom(self):
        """Test initialization with custom values."""
        units = Units(
            coords="Å",
            binding_radius="pm",
            D_trans="m^2/s",
            D_rot="deg^2/ms"
        )

        self.assertEqual(units.coords, "Å")
        self.assertEqual(units.binding_radius, "pm")
        self.assertEqual(units.D_trans, "m^2/s")
        self.assertEqual(units.D_rot, "deg^2/ms")

    def test_init_partial_custom(self):
        """Test initialization with some custom values."""
        units = Units(coords="Å", D_trans="nm^2/s")

        self.assertEqual(units.coords, "Å")
        self.assertEqual(units.binding_radius, "nm")  # default
        self.assertEqual(units.D_trans, "nm^2/s")
        self.assertEqual(units.D_rot, "rad^2/µs")  # default

    def test_to_dict_default(self):
        """Test dictionary serialization with default values."""
        result = self.default_units.to_dict()

        expected = {
            "coords": "nm",
            "binding_radius": "nm",
            "diffusion_translation": "nm^2/µs",
            "diffusion_rotation": "rad^2/µs"
        }

        self.assertEqual(result, expected)

    def test_to_dict_custom(self):
        """Test dictionary serialization with custom values."""
        result = self.custom_units.to_dict()

        expected = {
            "coords": "Å",
            "binding_radius": "nm",
            "diffusion_translation": "nm^2/s",
            "diffusion_rotation": "deg^2/s"
        }

        self.assertEqual(result, expected)

    def test_to_dict_key_mapping(self):
        """Test that to_dict uses correct key names for diffusion units."""
        units = Units(D_trans="test_trans", D_rot="test_rot")
        result = units.to_dict()

        # Verify the key mapping
        self.assertEqual(result["diffusion_translation"], "test_trans")
        self.assertEqual(result["diffusion_rotation"], "test_rot")
        self.assertNotIn("D_trans", result)
        self.assertNotIn("D_rot", result)

    def test_from_dict_empty(self):
        """Test creating instance from empty dictionary."""
        units = Units.from_dict({})

        # Should return default values
        self.assertEqual(units.coords, "nm")
        self.assertEqual(units.binding_radius, "nm")
        self.assertEqual(units.D_trans, "nm^2/µs")
        self.assertEqual(units.D_rot, "rad^2/µs")

    def test_from_dict_none(self):
        """Test creating instance from None."""
        units = Units.from_dict(None)

        # Should return default values
        self.assertEqual(units.coords, "nm")
        self.assertEqual(units.binding_radius, "nm")
        self.assertEqual(units.D_trans, "nm^2/µs")
        self.assertEqual(units.D_rot, "rad^2/µs")

    def test_from_dict_complete(self):
        """Test creating instance from complete dictionary."""
        data = {
            "coords": "Å",
            "binding_radius": "pm",
            "diffusion_translation": "m^2/s",
            "diffusion_rotation": "deg^2/ms"
        }

        units = Units.from_dict(data)

        self.assertEqual(units.coords, "Å")
        self.assertEqual(units.binding_radius, "pm")
        self.assertEqual(units.D_trans, "m^2/s")
        self.assertEqual(units.D_rot, "deg^2/ms")

    def test_from_dict_partial(self):
        """Test creating instance from partial dictionary."""
        data = {
            "coords": "Å",
            "diffusion_translation": "nm^2/s"
        }

        units = Units.from_dict(data)

        self.assertEqual(units.coords, "Å")
        self.assertEqual(units.binding_radius, "nm")  # default
        self.assertEqual(units.D_trans, "nm^2/s")
        self.assertEqual(units.D_rot, "rad^2/µs")  # default

    def test_from_dict_extra_keys(self):
        """Test creating instance from dictionary with extra keys."""
        data = {
            "coords": "nm",
            "binding_radius": "nm",
            "diffusion_translation": "nm^2/µs",
            "diffusion_rotation": "rad^2/µs",
            "extra_key": "ignored_value",
            "another_extra": 123
        }

        units = Units.from_dict(data)

        # Should ignore extra keys and work normally
        self.assertEqual(units.coords, "nm")
        self.assertEqual(units.binding_radius, "nm")
        self.assertEqual(units.D_trans, "nm^2/µs")
        self.assertEqual(units.D_rot, "rad^2/µs")

    def test_serialization_roundtrip(self):
        """Test that serialization and deserialization preserve data."""
        # Test with default units
        original = Units()
        dict_repr = original.to_dict()
        restored = Units.from_dict(dict_repr)

        self.assertEqual(original.coords, restored.coords)
        self.assertEqual(original.binding_radius, restored.binding_radius)
        self.assertEqual(original.D_trans, restored.D_trans)
        self.assertEqual(original.D_rot, restored.D_rot)

        # Test with custom units
        original_custom = Units(
            coords="Å",
            binding_radius="pm",
            D_trans="m^2/s",
            D_rot="deg^2/ms"
        )
        dict_repr_custom = original_custom.to_dict()
        restored_custom = Units.from_dict(dict_repr_custom)

        self.assertEqual(original_custom.coords, restored_custom.coords)
        self.assertEqual(original_custom.binding_radius,
                         restored_custom.binding_radius)
        self.assertEqual(original_custom.D_trans, restored_custom.D_trans)
        self.assertEqual(original_custom.D_rot, restored_custom.D_rot)


class TestUnitsStandards(unittest.TestCase):
    """Test units standards and conventions."""

    def test_default_units_are_ionerdss_standard(self):
        """Test that default units match ionerdss standards."""
        units = Units()

        # Test ionerdss standard units
        self.assertEqual(units.coords, "nm",
                         "Default coordinate unit should be nm")
        self.assertEqual(units.binding_radius, "nm",
                         "Default binding radius unit should be nm")
        self.assertEqual(units.D_trans, "nm^2/µs",
                         "Default translational diffusion unit should be nm^2/µs")
        self.assertEqual(units.D_rot, "rad^2/µs",
                         "Default rotational diffusion unit should be rad^2/µs")

    def test_pdb_compatibility_units(self):
        """Test units that would be compatible with PDB files."""
        pdb_units = Units(coords="Å")

        self.assertEqual(pdb_units.coords, "Å")
        # Other units should remain at ionerdss defaults
        self.assertEqual(pdb_units.binding_radius, "nm")
        self.assertEqual(pdb_units.D_trans, "nm^2/µs")
        self.assertEqual(pdb_units.D_rot, "rad^2/µs")

    def test_common_alternative_units(self):
        """Test common alternative unit combinations."""
        # SI base units
        si_units = Units(
            coords="m",
            binding_radius="m",
            D_trans="m^2/s",
            D_rot="rad^2/s"
        )

        self.assertEqual(si_units.coords, "m")
        self.assertEqual(si_units.D_trans, "m^2/s")

        # Angstrom-based units (common in structural biology)
        angstrom_units = Units(
            coords="Å",
            binding_radius="Å",
            D_trans="Å^2/ns",
            D_rot="rad^2/ns"
        )

        self.assertEqual(angstrom_units.coords, "Å")
        self.assertEqual(angstrom_units.D_trans, "Å^2/ns")


class TestUnitsEdgeCases(unittest.TestCase):
    """Test edge cases and special scenarios."""

    def test_empty_string_units(self):
        """Test behavior with empty string units."""
        units = Units(
            coords="",
            binding_radius="",
            D_trans="",
            D_rot=""
        )

        self.assertEqual(units.coords, "")
        self.assertEqual(units.binding_radius, "")
        self.assertEqual(units.D_trans, "")
        self.assertEqual(units.D_rot, "")

    def test_unicode_units(self):
        """Test units with Unicode characters."""
        unicode_units = Units(
            coords="nm",
            binding_radius="nm",
            D_trans="nm²/μs",  # Using Unicode superscript and mu
            D_rot="rad²/μs"
        )

        self.assertEqual(unicode_units.D_trans, "nm²/μs")
        self.assertEqual(unicode_units.D_rot, "rad²/μs")

    def test_long_unit_strings(self):
        """Test with very long unit strings."""
        long_unit = "very_long_unit_name_that_might_cause_issues_in_some_systems"
        units = Units(coords=long_unit)

        self.assertEqual(units.coords, long_unit)

    def test_special_characters_in_units(self):
        """Test units with special characters."""
        special_units = Units(
            coords="nm",
            binding_radius="nm",
            D_trans="nm^2/µs·K^-1",  # Complex unit with multiple symbols
            D_rot="rad^2/(µs·mol)"
        )

        self.assertEqual(special_units.D_trans, "nm^2/µs·K^-1")
        self.assertEqual(special_units.D_rot, "rad^2/(µs·mol)")


class TestUnitsDocumentation(unittest.TestCase):
    """Test that units behave as documented."""

    def test_pdb_conversion_note(self):
        """Test the documented PDB conversion requirement."""
        # Default ionerdss units
        ionerdss_units = Units()
        self.assertEqual(ionerdss_units.coords, "nm")

        # PDB-compatible units
        pdb_units = Units(coords="Å")
        self.assertEqual(pdb_units.coords, "Å")

        # Document the conversion factor (1 nm = 10 Å)
        # This test serves as documentation of the conversion requirement
        nm_value = 1.0
        angstrom_value = nm_value * 10.0

        self.assertEqual(angstrom_value, 10.0)
        self.assertEqual(nm_value, angstrom_value / 10.0)

    def test_diffusion_time_units(self):
        """Test different time units for diffusion constants."""
        # Default: microseconds
        default_units = Units()
        self.assertIn("µs", default_units.D_trans)
        self.assertIn("µs", default_units.D_rot)

        # Alternative: seconds
        second_units = Units(
            D_trans="nm^2/s",
            D_rot="rad^2/s"
        )
        self.assertIn("/s", second_units.D_trans)
        self.assertIn("/s", second_units.D_rot)

    def test_angle_units(self):
        """Test different angle units for rotational diffusion."""
        # Default: radians
        default_units = Units()
        self.assertIn("rad", default_units.D_rot)

        # Alternative: degrees
        degree_units = Units(D_rot="deg^2/µs")
        self.assertIn("deg", degree_units.D_rot)


class TestUnitsIntegration(unittest.TestCase):
    """Integration tests for Units with other components."""

    def test_configuration_file_format(self):
        """Test units in configuration file format."""
        # Simulate a configuration dictionary
        config = {
            "units": {
                "coords": "nm",
                "binding_radius": "nm",
                "diffusion_translation": "nm^2/µs",
                "diffusion_rotation": "rad^2/µs"
            },
            "other_config": "value"
        }

        units = Units.from_dict(config["units"])

        self.assertEqual(units.coords, "nm")
        self.assertEqual(units.D_trans, "nm^2/µs")

        # Test serialization back to config format
        serialized = {"units": units.to_dict()}
        self.assertEqual(serialized["units"]["coords"], "nm")
        self.assertEqual(serialized["units"]
                         ["diffusion_translation"], "nm^2/µs")

    def test_units_consistency_check(self):
        """Test that units are consistent across related quantities."""
        units = Units()

        # Coordinate-related units should be consistent
        self.assertEqual(units.coords, "nm")
        self.assertEqual(units.binding_radius, "nm")

        # Diffusion units should have consistent length scale
        self.assertIn("nm", units.D_trans)
        # Time scale should be consistent between translational and rotational
        if "µs" in units.D_trans:
            self.assertIn("µs", units.D_rot)


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)

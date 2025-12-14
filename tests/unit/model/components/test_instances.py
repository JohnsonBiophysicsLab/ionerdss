"""
Unit tests for ionerdss.model.components.instances module.

Tests the InterfaceInstance and MoleculeInstance classes including their
serialization, naming conventions, and relationships.
"""

import unittest
from unittest.mock import Mock
import numpy as np

# Import the module under test
from ionerdss.model.components.instances import (
    InterfaceInstance,
    MoleculeInstance,
    _replace_underscore_with_dash
)
from ionerdss.model.components.types import MoleculeType, InterfaceType


class TestHelperFunctions(unittest.TestCase):
    """Test helper functions in the instances module."""

    def test_replace_underscore_with_dash(self):
        """Test underscore replacement function."""
        # Basic replacement
        self.assertEqual(_replace_underscore_with_dash("A_12"), "A-12")

        # Multiple underscores
        self.assertEqual(_replace_underscore_with_dash("A_B_C"), "A-B-C")

        # No underscores
        self.assertEqual(_replace_underscore_with_dash("ABC"), "ABC")

        # Empty string
        self.assertEqual(_replace_underscore_with_dash(""), "")

        # Only underscores
        self.assertEqual(_replace_underscore_with_dash("___"), "---")


class TestInterfaceInstance(unittest.TestCase):
    """Test the InterfaceInstance class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock objects
        self.mock_interface_type = Mock(spec=InterfaceType)
        self.mock_interface_type.get_name.return_value = "A_B_1"

        self.mock_molecule_instance = Mock(spec=MoleculeInstance)
        self.mock_partner_interface = Mock(spec=InterfaceInstance)

        # Create test coordinates
        self.test_coord = np.array([1.0, 2.0, 3.0])

        # Create a basic interface instance (only absolute_coord is required)
        self.interface = InterfaceInstance(
            absolute_coord=self.test_coord,
            partner_interface=self.mock_partner_interface,
            interface_type=self.mock_interface_type,
            this_mol=self.mock_molecule_instance,
            this_mol_name="A",
            partner_mol_name="B",
            interface_index=1,
            residues=[1, 2, 3],
            energy=-5.0
        )

    def test_init_minimal(self):
        """Test initialization with only required field."""
        interface = InterfaceInstance(absolute_coord=self.test_coord)

        np.testing.assert_array_equal(
            interface.absolute_coord, self.test_coord)
        self.assertIsNone(interface.partner_interface)
        self.assertIsNone(interface.interface_type)
        self.assertIsNone(interface.this_mol)
        self.assertEqual(interface.this_mol_name, "unnamed")
        self.assertEqual(interface.partner_mol_name, "unnamed")
        self.assertEqual(interface.interface_index, 0)
        self.assertEqual(interface.residues, [])
        self.assertEqual(interface.energy, -1.0)
        self.assertEqual(interface.signature, {})

    def test_get_name(self):
        """Test interface name generation."""
        name = self.interface.get_name()
        self.assertEqual(name, "A_B_1")

    def test_get_name_with_different_values(self):
        """Test name generation with different values."""
        interface = InterfaceInstance(
            absolute_coord=self.test_coord,
            this_mol_name="ProteinX",
            partner_mol_name="ProteinY",
            interface_index=42
        )

        name = interface.get_name()
        self.assertEqual(name, "ProteinX_ProteinY_42")

    def test_set_name_valid(self):
        """Test setting name from valid string."""
        self.interface.set_name("C_D_5")

        self.assertEqual(self.interface.this_mol_name, "C")
        self.assertEqual(self.interface.partner_mol_name, "D")
        self.assertEqual(self.interface.interface_index, 5)

    def test_set_name_invalid_format(self):
        """Test setting name with invalid format."""
        with self.assertRaises(IndexError):
            self.interface.set_name("A_B")  # Missing index

        with self.assertRaises(IndexError):
            self.interface.set_name("A")  # Missing partner and index

    def test_set_name_invalid_index(self):
        """Test setting name with non-numeric index."""
        with self.assertRaises(ValueError):
            self.interface.set_name("A_B_invalid")

    def test_to_dict(self):
        """Test dictionary serialization."""
        result = self.interface.to_dict()

        # Check that result contains expected keys
        self.assertEqual(result["name"], "A_B_1")
        self.assertEqual(result["type"], "A_B_1")
        self.assertEqual(result["coord"], [1.0, 2.0, 3.0])
        self.assertEqual(result["residues"], [1, 2, 3])
        self.assertEqual(result["energy"], -5.0)
        self.mock_interface_type.get_name.assert_called_once()

    def test_to_dict_no_interface_type(self):
        """Test dictionary serialization when interface_type is None."""
        interface = InterfaceInstance(absolute_coord=self.test_coord)
        result = interface.to_dict()

        self.assertEqual(result["type"], "unknown")

    def test_from_dict_complete(self):
        """Test creating instance from complete dictionary."""
        data = {
            "this_mol_name": "X",
            "partner_mol_name": "Y",
            "interface_index": 3,
            "this_mol": self.mock_molecule_instance,
            "partner_interface": self.mock_partner_interface,
            "interface_type": self.mock_interface_type,
            "absolute_coords": [1.0, 2.0, 3.0],
            "residues": [4, 5, 6],
            "energy": -10.0,
            "signature": {"key": "value"}
        }

        interface = InterfaceInstance.from_dict(data)

        self.assertEqual(interface.this_mol_name, "X")
        self.assertEqual(interface.partner_mol_name, "Y")
        self.assertEqual(interface.interface_index, 3)
        self.assertEqual(interface.this_mol, self.mock_molecule_instance)
        self.assertEqual(interface.partner_interface,
                         self.mock_partner_interface)
        self.assertEqual(interface.interface_type, self.mock_interface_type)
        np.testing.assert_array_equal(
            interface.absolute_coord, np.array([1.0, 2.0, 3.0]))
        self.assertEqual(interface.residues, [4, 5, 6])
        self.assertEqual(interface.energy, -10.0)
        self.assertEqual(interface.signature, {"key": "value"})

    def test_from_dict_minimal(self):
        """Test creating instance from minimal dictionary."""
        data = {
            "absolute_coords": [1.0, 2.0, 3.0]
        }

        interface = InterfaceInstance.from_dict(data)

        np.testing.assert_array_equal(
            interface.absolute_coord, np.array([1.0, 2.0, 3.0]))
        self.assertEqual(interface.this_mol_name, "unnamed")
        self.assertEqual(interface.partner_mol_name, "unnamed")
        self.assertEqual(interface.interface_index, 0)
        self.assertIsNone(interface.this_mol)
        self.assertIsNone(interface.partner_interface)
        self.assertIsNone(interface.interface_type)
        self.assertEqual(interface.energy, -1.0)
        self.assertEqual(interface.signature, {})

    def test_from_dict_missing_required_key(self):
        """Test from_dict with missing required key."""
        data = {"this_mol_name": "A"}

        with self.assertRaises(KeyError):
            InterfaceInstance.from_dict(data)


class TestMoleculeInstance(unittest.TestCase):
    """Test the MoleculeInstance class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock objects
        self.mock_molecule_type = Mock(spec=MoleculeType)
        self.mock_molecule_type.name = "TestMolecule"

        # Create test coordinates
        self.test_com = np.array([0.0, 1.0, 2.0])
        self.test_norm = np.array([0.0, 0.0, 1.0])
        self.test_ref1 = np.array([1.0, 0.0, 0.0])
        self.test_ref2 = np.array([0.0, 0.0, 1.0])

        # Create mock interfaces
        self.mock_interface1 = Mock(spec=InterfaceInstance)
        self.mock_interface1.get_name.return_value = "A_B_1"

        self.mock_interface2 = Mock(spec=InterfaceInstance)
        self.mock_interface2.get_name.return_value = "A_C_1"

        # Create mock partner molecules with name attribute
        self.mock_partner1 = Mock(spec=MoleculeInstance)
        self.mock_partner1.name = "PartnerMol1"
        
        self.mock_partner2 = Mock(spec=MoleculeInstance)
        self.mock_partner2.name = "PartnerMol2"

        # Create molecule instance
        self.molecule = MoleculeInstance(
            name="TestMol_001",
            norm=self.test_norm,
            ref1=self.test_ref1,
            ref2=self.test_ref2,
            com=self.test_com,
            molecule_type=self.mock_molecule_type,
            interfaces_neighbors_map={
                self.mock_interface1: self.mock_partner1,
                self.mock_interface2: self.mock_partner2
            }
        )

    def test_init_minimal(self):
        """Test initialization with only required fields."""
        molecule = MoleculeInstance(
            name="Test",
            norm=self.test_norm,
            ref1=self.test_ref1,
            ref2=self.test_ref2,
            com=self.test_com
        )

        self.assertEqual(molecule.name, "Test")
        np.testing.assert_array_equal(molecule.norm, self.test_norm)
        np.testing.assert_array_equal(molecule.ref1, self.test_ref1)
        np.testing.assert_array_equal(molecule.ref2, self.test_ref2)
        np.testing.assert_array_equal(molecule.com, self.test_com)
        self.assertIsNone(molecule.molecule_type)
        self.assertEqual(molecule.interfaces_neighbors_map, {})

    def test_init_complete(self):
        """Test initialization with all fields."""
        self.assertEqual(self.molecule.name, "TestMol_001")
        self.assertEqual(self.molecule.molecule_type, self.mock_molecule_type)
        np.testing.assert_array_equal(self.molecule.norm, self.test_norm)
        np.testing.assert_array_equal(self.molecule.com, self.test_com)
        self.assertEqual(len(self.molecule.interfaces_neighbors_map), 2)

    def test_to_dict(self):
        """Test dictionary serialization."""
        result = self.molecule.to_dict()

        # Just check the basic structure is correct
        self.assertEqual(result["name"], "TestMol_001")
        self.assertEqual(result["type"], "TestMolecule")
        self.assertEqual(result["com"], [0.0, 1.0, 2.0])
        self.assertEqual(result["norm"], [0.0, 0.0, 1.0])
        self.assertEqual(result["ref1"], [1.0, 0.0, 0.0])
        self.assertEqual(result["ref2"], [0.0, 0.0, 1.0])
        self.assertIn("interfaces", result)
        self.assertEqual(len(result["interfaces"]), 2)

        self.mock_interface1.get_name.assert_called_once()
        self.mock_interface2.get_name.assert_called_once()

    def test_to_dict_no_molecule_type(self):
        """Test dictionary serialization when molecule_type is None."""
        molecule = MoleculeInstance(
            name="Test",
            norm=self.test_norm,
            ref1=self.test_ref1,
            ref2=self.test_ref2,
            com=self.test_com
        )
        result = molecule.to_dict()

        self.assertEqual(result["type"], "unknown")

    def test_from_dict_complete(self):
        """Test creating instance from complete dictionary."""
        data = {
            "name": "NewMol_002",
            # interfaces_neighbors_map is intentionally omitted - it's rebuilt by system
            "molecule_type": self.mock_molecule_type,
            "norm": [1.0, 0.0, 0.0],
            "ref1": [1.0, 0.0, 0.0],
            "ref2": [0.0, 0.0, 1.0],
            "coord": [5.0, 6.0, 7.0]
        }

        molecule = MoleculeInstance.from_dict(data)

        self.assertEqual(molecule.name, "NewMol-002")  # Underscore replaced
        # interfaces_neighbors_map is always empty from from_dict - rebuilt later by system
        self.assertEqual(molecule.interfaces_neighbors_map, {})
        self.assertEqual(molecule.molecule_type, self.mock_molecule_type)
        np.testing.assert_array_equal(molecule.norm, np.array([1.0, 0.0, 0.0]))
        np.testing.assert_array_equal(molecule.com, np.array([5.0, 6.0, 7.0]))

    def test_from_dict_minimal(self):
        """Test creating instance from minimal dictionary."""
        # This should work now since all required fields have defaults in from_dict
        data = {}

        molecule = MoleculeInstance.from_dict(data)

        self.assertEqual(molecule.name, "unnamed")
        self.assertIsNone(molecule.molecule_type)
        np.testing.assert_array_equal(molecule.norm, np.array([0.0, 0.0, 1.0]))
        np.testing.assert_array_equal(molecule.com, np.array([0.0, 0.0, 0.0]))
        self.assertEqual(molecule.interfaces_neighbors_map, {})

    def test_from_dict_underscore_replacement(self):
        """Test that underscores in names are replaced with dashes."""
        data = {"name": "Protein_A_123"}

        molecule = MoleculeInstance.from_dict(data)

        self.assertEqual(molecule.name, "Protein-A-123")

    def test_from_dict_with_numpy_arrays(self):
        """Test from_dict properly handles numpy array creation."""
        data = {
            "norm": [0.5, 0.5, 0.707],
            "coord": [10.0, 20.0, 30.0]
        }

        molecule = MoleculeInstance.from_dict(data)

        self.assertIsInstance(molecule.norm, np.ndarray)
        self.assertIsInstance(molecule.com, np.ndarray)
        np.testing.assert_array_almost_equal(molecule.norm, [0.5, 0.5, 0.707])
        np.testing.assert_array_equal(molecule.com, [10.0, 20.0, 30.0])


class TestIntegration(unittest.TestCase):
    """Integration tests for instances working together."""

    def test_interface_molecule_relationship(self):
        """Test that interfaces and molecules reference each other correctly."""
        # Create molecule type
        mol_type = Mock(spec=MoleculeType)
        mol_type.name = "TestType"

        # Create molecule instance
        molecule = MoleculeInstance(
            name="TestMol",
            norm=np.array([0.0, 0.0, 1.0]),
            ref1=np.array([1.0, 0.0, 0.0]),
            ref2=np.array([0.0, 0.0, 1.0]),
            com=np.array([0.0, 0.0, 0.0]),
            molecule_type=mol_type
        )

        # Create interface type
        interface_type = Mock(spec=InterfaceType)
        interface_type.get_name.return_value = "A_B_1"

        # Create interface instance
        interface = InterfaceInstance(
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            this_mol_name="A",
            partner_mol_name="B",
            interface_index=1,
            partner_interface=None,  # Will set later
            interface_type=interface_type,
            this_mol=molecule
        )

        # Test relationships
        self.assertEqual(interface.this_mol, molecule)
        self.assertEqual(interface.interface_type, interface_type)
        self.assertEqual(interface.get_name(), "A_B_1")

    def test_bidirectional_interface_binding(self):
        """Test that two interfaces can reference each other."""
        # Create two interface instances
        interface1 = InterfaceInstance(
            absolute_coord=np.array([0.0, 0.0, 0.0]),
            this_mol_name="A",
            partner_mol_name="B",
            interface_index=1,
            partner_interface=None,  # Set later
            interface_type=Mock(spec=InterfaceType),
            this_mol=Mock(spec=MoleculeInstance)
        )

        interface2 = InterfaceInstance(
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            this_mol_name="B",
            partner_mol_name="A",
            interface_index=1,
            partner_interface=interface1,
            interface_type=Mock(spec=InterfaceType),
            this_mol=Mock(spec=MoleculeInstance)
        )

        # Set up bidirectional relationship
        interface1.partner_interface = interface2

        # Test bidirectional references
        self.assertEqual(interface1.partner_interface, interface2)
        self.assertEqual(interface2.partner_interface, interface1)


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)

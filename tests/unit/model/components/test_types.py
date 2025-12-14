"""
Unit tests for ionerdss.model.components.types module.

Tests the InterfaceType and MoleculeType classes including their
serialization, naming conventions, and physical property calculations.
"""

import unittest
from unittest.mock import Mock, patch
import numpy as np

# Import the module under test
from ionerdss.model.components.types import InterfaceType, MoleculeType


class TestInterfaceType(unittest.TestCase):
    """Test the InterfaceType class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create test coordinates
        self.test_absolute_coord = np.array([1.0, 2.0, 3.0])
        self.test_local_coord = np.array([0.5, 1.0, 1.5])

        # Create basic interface type
        self.interface_type = InterfaceType(
            this_mol_type_name="ProteinA",
            partner_mol_type_name="ProteinB",
            interface_index=1,
            absolute_coord=self.test_absolute_coord,
            local_coord=self.test_local_coord,
            energy=-5.0,
            required_free=["A_C_1", "A_D_1"]
        )

    def test_init_required_fields(self):
        """Test initialization with only required fields."""
        interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=self.test_absolute_coord,
            local_coord=self.test_local_coord
        )

        self.assertEqual(interface.this_mol_type_name, "A")
        self.assertEqual(interface.partner_mol_type_name, "B")
        self.assertEqual(interface.interface_index, 1)
        np.testing.assert_array_equal(
            interface.absolute_coord, self.test_absolute_coord)
        np.testing.assert_array_equal(
            interface.local_coord, self.test_local_coord)
        self.assertIsNone(interface.partner_interface_type)
        self.assertIsNone(interface.this_mol_type)
        self.assertIsNone(interface.partner_mol_type)
        self.assertEqual(interface.required_free, [])
        self.assertEqual(interface.energy, -1.0)
        self.assertEqual(interface.signature, {})

    def test_init_all_fields(self):
        """Test initialization with all fields."""
        mock_partner_interface = Mock(spec=InterfaceType)
        mock_this_mol_type = Mock(spec=MoleculeType)
        mock_partner_mol_type = Mock(spec=MoleculeType)

        interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=2,
            absolute_coord=self.test_absolute_coord,
            local_coord=self.test_local_coord,
            partner_interface_type=mock_partner_interface,
            this_mol_type=mock_this_mol_type,
            partner_mol_type=mock_partner_mol_type,
            required_free=["test"],
            energy=-10.0,
            signature={"key": "value"}
        )

        self.assertEqual(interface.partner_interface_type,
                         mock_partner_interface)
        self.assertEqual(interface.this_mol_type, mock_this_mol_type)
        self.assertEqual(interface.partner_mol_type, mock_partner_mol_type)
        self.assertEqual(interface.required_free, ["test"])
        self.assertEqual(interface.energy, -10.0)
        self.assertEqual(interface.signature, {"key": "value"})

    def test_get_name(self):
        """Test interface name generation."""
        name = self.interface_type.get_name()
        self.assertEqual(name, "ProteinA_ProteinB_1")

    def test_get_name_different_values(self):
        """Test name generation with different values."""
        interface = InterfaceType(
            this_mol_type_name="X",
            partner_mol_type_name="Y",
            interface_index=42,
            absolute_coord=self.test_absolute_coord,
            local_coord=self.test_local_coord
        )

        name = interface.get_name()
        self.assertEqual(name, "X_Y_42")

    def test_set_name_valid(self):
        """Test setting name from valid string."""
        self.interface_type.set_name("NewMol_NewPartner_5")

        self.assertEqual(self.interface_type.this_mol_type_name, "NewMol")
        self.assertEqual(
            self.interface_type.partner_mol_type_name, "NewPartner")
        self.assertEqual(self.interface_type.interface_index, 5)

    def test_set_name_invalid_format(self):
        """Test setting name with invalid format."""
        with self.assertRaises(IndexError):
            self.interface_type.set_name("A_B")  # Missing index

        with self.assertRaises(IndexError):
            self.interface_type.set_name("A")  # Missing partner and index

    def test_set_name_invalid_index(self):
        """Test setting name with non-numeric index."""
        with self.assertRaises(ValueError):
            self.interface_type.set_name("A_B_invalid")

    def test_to_dict_minimal(self):
        """Test dictionary serialization with minimal data."""
        interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=self.test_absolute_coord,
            local_coord=self.test_local_coord
        )

        result = interface.to_dict()

        expected = {
            "name": "A_B_1",
            "partner_interface_type": None,
            "this_mol_type": None,
            "absolute_coord": [1.0, 2.0, 3.0],
            "local_coord": [0.5, 1.0, 1.5],
            "required_free": [],
            "energy": -1.0,
            "signature": {}
        }

        self.assertEqual(result, expected)

    def test_to_dict_complete(self):
        """Test dictionary serialization with complete data."""
        # Set up mocks
        mock_partner_interface = Mock(spec=InterfaceType)
        mock_partner_interface.get_name.return_value = "B_A_1"

        mock_this_mol_type = Mock(spec=MoleculeType)
        mock_this_mol_type.name = "MolTypeA"

        self.interface_type.partner_interface_type = mock_partner_interface
        self.interface_type.this_mol_type = mock_this_mol_type

        result = self.interface_type.to_dict()

        self.assertEqual(result["name"], "ProteinA_ProteinB_1")
        self.assertEqual(result["partner_interface_type"], "B_A_1")
        self.assertEqual(result["this_mol_type"], "MolTypeA")
        self.assertEqual(result["absolute_coord"], [1.0, 2.0, 3.0])
        self.assertEqual(result["local_coord"], [0.5, 1.0, 1.5])
        self.assertEqual(result["required_free"], ["A_C_1", "A_D_1"])
        self.assertEqual(result["energy"], -5.0)
        self.assertEqual(result["signature"], {})

    def test_from_dict_minimal(self):
        """Test creating instance from minimal dictionary."""
        data = {
            "this_mol_type_name": "X",
            "partner_mol_type_name": "Y",
            "interface_index": 3,
            "absolute_coord": [2.0, 3.0, 4.0],
            "local_coord": [1.0, 1.5, 2.0]
        }

        interface = InterfaceType.from_dict(data)

        self.assertEqual(interface.this_mol_type_name, "X")
        self.assertEqual(interface.partner_mol_type_name, "Y")
        self.assertEqual(interface.interface_index, 3)
        np.testing.assert_array_equal(
            interface.absolute_coord, np.array([2.0, 3.0, 4.0]))
        np.testing.assert_array_equal(
            interface.local_coord, np.array([1.0, 1.5, 2.0]))
        self.assertIsNone(interface.partner_interface_type)
        self.assertEqual(interface.required_free, [])
        self.assertEqual(interface.energy, -1.0)
        self.assertEqual(interface.signature, {})

    def test_from_dict_complete(self):
        """Test creating instance from complete dictionary."""
        mock_partner_interface = Mock(spec=InterfaceType)
        mock_this_mol_type = Mock(spec=MoleculeType)
        mock_partner_mol_type = Mock(spec=MoleculeType)

        data = {
            "this_mol_type_name": "A",
            "partner_mol_type_name": "B",
            "interface_index": 2,
            "absolute_coord": [1.0, 2.0, 3.0],
            "local_coord": [0.5, 1.0, 1.5],
            "partner_interface_type": mock_partner_interface,
            "this_mol_type": mock_this_mol_type,
            "partner_mol_type": mock_partner_mol_type,
            "required_free": ["test1", "test2"],
            "energy": -15.0,
            "signature": {"test": "data"}
        }

        interface = InterfaceType.from_dict(data)

        self.assertEqual(interface.partner_interface_type,
                         mock_partner_interface)
        self.assertEqual(interface.this_mol_type, mock_this_mol_type)
        self.assertEqual(interface.partner_mol_type, mock_partner_mol_type)
        self.assertEqual(interface.required_free, ["test1", "test2"])
        self.assertEqual(interface.energy, -15.0)
        self.assertEqual(interface.signature, {"test": "data"})

    def test_from_dict_missing_required_keys(self):
        """Test from_dict with missing required keys."""
        # Missing this_mol_type_name
        data = {
            "partner_mol_type_name": "B",
            "interface_index": 1,
            "absolute_coord": [1.0, 2.0, 3.0],
            "local_coord": [0.5, 1.0, 1.5]
        }

        with self.assertRaises(KeyError):
            InterfaceType.from_dict(data)

        # Missing absolute_coord
        data = {
            "this_mol_type_name": "A",
            "partner_mol_type_name": "B",
            "interface_index": 1,
            "local_coord": [0.5, 1.0, 1.5]
        }

        with self.assertRaises(KeyError):
            InterfaceType.from_dict(data)


class TestMoleculeType(unittest.TestCase):
    """Test the MoleculeType class."""

    def setUp(self):
        """Set up test fixtures."""
        self.molecule_type = MoleculeType(
            name="TestProtein",
            radius_nm=2.5,
            D_t_nm2_us=1.0,
            D_r_rad2_us=0.5,
            interfaces_neighbors_map={
                "interface1": "neighbor1", "interface2": "neighbor2"}
        )

    def test_init_minimal(self):
        """Test initialization with only required fields."""
        mol_type = MoleculeType(name="MinimalProtein")

        self.assertEqual(mol_type.name, "MinimalProtein")
        self.assertEqual(mol_type.interfaces_neighbors_map, {})
        self.assertEqual(mol_type.radius_nm, 0.0)
        self.assertEqual(mol_type.D_t_nm2_us, 0.0)
        self.assertEqual(mol_type.D_r_rad2_us, 0.0)

    def test_init_complete(self):
        """Test initialization with all fields."""
        self.assertEqual(self.molecule_type.name, "TestProtein")
        self.assertEqual(self.molecule_type.radius_nm, 2.5)
        self.assertEqual(self.molecule_type.D_t_nm2_us, 1.0)
        self.assertEqual(self.molecule_type.D_r_rad2_us, 0.5)
        self.assertEqual(len(self.molecule_type.interfaces_neighbors_map), 2)

    @patch('ionerdss.model.components.types.compute_diffusion_constants_nm_us')
    def test_set_diffusion_constants_from_radius(self, mock_compute):
        """Test diffusion constant calculation from radius."""
        # Set up mock return values
        mock_compute.return_value = (1.5, 0.75)

        mol_type = MoleculeType(name="TestMol", radius_nm=3.0)
        mol_type.set_diffusion_constants_from_radius()

        # Verify the function was called with correct parameters
        mock_compute.assert_called_once_with(radius_nm=3.0)

        # Verify the values were set
        self.assertEqual(mol_type.D_t_nm2_us, 1.5)
        self.assertEqual(mol_type.D_r_rad2_us, 0.75)

    def test_to_dict(self):
        """Test dictionary serialization."""
        result = self.molecule_type.to_dict()

        expected = {
            "name": "TestProtein",
            "radius": 2.5,
            "diffusion_translation": 1.0,
            "diffusion_rotation": 0.5,
            "interfaces_neighbors_map": {"interface1": "neighbor1", "interface2": "neighbor2"},
            "ref1_local": [1.0, 0.0, 0.0],
            "ref2_local": [0.0, 0.0, 1.0]
        }

        self.assertEqual(result, expected)

    def test_to_dict_minimal(self):
        """Test dictionary serialization with minimal data."""
        mol_type = MoleculeType(name="MinimalMol")
        result = mol_type.to_dict()

        expected = {
            "name": "MinimalMol",
            "radius": 0.0,
            "diffusion_translation": 0.0,
            "diffusion_rotation": 0.0,
            "interfaces_neighbors_map": {},
            "ref1_local": [1.0, 0.0, 0.0],
            "ref2_local": [0.0, 0.0, 1.0]
        }

        self.assertEqual(result, expected)

    def test_from_dict_minimal(self):
        """Test creating instance from minimal dictionary."""
        data = {"name": "FromDictMol"}

        mol_type = MoleculeType.from_dict(data)

        self.assertEqual(mol_type.name, "FromDictMol")
        self.assertEqual(mol_type.radius_nm, 0.0)
        self.assertEqual(mol_type.D_t_nm2_us, 0.0)
        self.assertEqual(mol_type.D_r_rad2_us, 0.0)
        self.assertEqual(mol_type.interfaces_neighbors_map, {})

    def test_from_dict_complete(self):
        """Test creating instance from complete dictionary."""
        data = {
            "name": "CompleteFromDict",
            "radius": 3.5,
            "diffusion_translation": 2.0,
            "diffusion_rotation": 1.0,
            "interfaces_neighbors_map": {"int1": "neigh1"}
        }

        mol_type = MoleculeType.from_dict(data)

        self.assertEqual(mol_type.name, "CompleteFromDict")
        self.assertEqual(mol_type.radius_nm, 3.5)
        self.assertEqual(mol_type.D_t_nm2_us, 2.0)
        self.assertEqual(mol_type.D_r_rad2_us, 1.0)
        self.assertEqual(mol_type.interfaces_neighbors_map, {"int1": "neigh1"})

    def test_from_dict_missing_name(self):
        """Test from_dict with missing required name field."""
        data = {"radius": 2.0}

        with self.assertRaises(KeyError):
            MoleculeType.from_dict(data)

    def test_from_dict_type_conversion(self):
        """Test that from_dict properly converts string numbers to floats."""
        data = {
            "name": "TypeConversion",
            "radius": "2.5",
            "diffusion_translation": "1.5",
            "diffusion_rotation": "0.75"
        }

        mol_type = MoleculeType.from_dict(data)

        self.assertIsInstance(mol_type.radius_nm, float)
        self.assertIsInstance(mol_type.D_t_nm2_us, float)
        self.assertIsInstance(mol_type.D_r_rad2_us, float)
        self.assertEqual(mol_type.radius_nm, 2.5)
        self.assertEqual(mol_type.D_t_nm2_us, 1.5)
        self.assertEqual(mol_type.D_r_rad2_us, 0.75)


class TestIntegration(unittest.TestCase):
    """Integration tests for types working together."""

    def test_interface_molecule_relationship(self):
        """Test that interfaces and molecule types can reference each other."""
        # Create molecule types
        mol_type_A = MoleculeType(name="MolA", radius_nm=2.0)
        mol_type_B = MoleculeType(name="MolB", radius_nm=1.5)

        # Create interface types
        interface_A_B = InterfaceType(
            this_mol_type_name="MolA",
            partner_mol_type_name="MolB",
            interface_index=1,
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            local_coord=np.array([0.5, 0.0, 0.0]),
            this_mol_type=mol_type_A,
            partner_mol_type=mol_type_B,
            energy=-10.0
        )

        interface_B_A = InterfaceType(
            this_mol_type_name="MolB",
            partner_mol_type_name="MolA",
            interface_index=1,
            absolute_coord=np.array([-1.0, 0.0, 0.0]),
            local_coord=np.array([-0.5, 0.0, 0.0]),
            this_mol_type=mol_type_B,
            partner_mol_type=mol_type_A,
            partner_interface_type=interface_A_B,
            energy=-10.0
        )

        # Set up bidirectional relationship
        interface_A_B.partner_interface_type = interface_B_A

        # Test relationships
        self.assertEqual(interface_A_B.this_mol_type, mol_type_A)
        self.assertEqual(interface_A_B.partner_mol_type, mol_type_B)
        self.assertEqual(interface_A_B.partner_interface_type, interface_B_A)
        self.assertEqual(interface_B_A.partner_interface_type, interface_A_B)

        # Test naming
        self.assertEqual(interface_A_B.get_name(), "MolA_MolB_1")
        self.assertEqual(interface_B_A.get_name(), "MolB_MolA_1")

    def test_serialization_roundtrip_molecule(self):
        """Test that molecule serialization and deserialization preserve data."""
        # Create original molecule type
        original_mol = MoleculeType(
            name="RoundtripTest",
            radius_nm=2.5,
            D_t_nm2_us=1.0,
            D_r_rad2_us=0.5,
            interfaces_neighbors_map={"test": "neighbor"}
        )

        # Serialize to dict
        mol_dict = original_mol.to_dict()

        # Deserialize back
        restored_mol = MoleculeType.from_dict(mol_dict)

        # Verify all fields match
        self.assertEqual(original_mol.name, restored_mol.name)
        self.assertEqual(original_mol.radius_nm, restored_mol.radius_nm)
        self.assertEqual(original_mol.D_t_nm2_us, restored_mol.D_t_nm2_us)
        self.assertEqual(original_mol.D_r_rad2_us, restored_mol.D_r_rad2_us)
        self.assertEqual(original_mol.interfaces_neighbors_map,
                         restored_mol.interfaces_neighbors_map)

    def test_serialization_roundtrip_interface(self):
        """Test that interface serialization and deserialization preserve basic data."""
        # Create original interface type
        original_interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([1.0, 2.0, 3.0]),
            local_coord=np.array([0.5, 1.0, 1.5]),
            energy=-5.0,
            required_free=["test1", "test2"],
            signature={"test": "data"}
        )

        # Create interface from basic data (simulating deserialization)
        restored_interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([1.0, 2.0, 3.0]),
            local_coord=np.array([0.5, 1.0, 1.5]),
            energy=-5.0,
            required_free=["test1", "test2"],
            signature={"test": "data"}
        )

        # Verify basic fields match
        self.assertEqual(original_interface.this_mol_type_name,
                         restored_interface.this_mol_type_name)
        self.assertEqual(original_interface.partner_mol_type_name,
                         restored_interface.partner_mol_type_name)
        self.assertEqual(original_interface.interface_index,
                         restored_interface.interface_index)
        np.testing.assert_array_equal(
            original_interface.absolute_coord, restored_interface.absolute_coord)
        np.testing.assert_array_equal(
            original_interface.local_coord, restored_interface.local_coord)
        self.assertEqual(original_interface.energy, restored_interface.energy)
        self.assertEqual(original_interface.required_free,
                         restored_interface.required_free)
        self.assertEqual(original_interface.signature,
                         restored_interface.signature)

    def test_bidirectional_interface_binding(self):
        """Test that two interface types can reference each other."""
        # Create interface types
        interface1 = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([0.0, 0.0, 0.0]),
            local_coord=np.array([1.0, 0.0, 0.0]),
            partner_interface_type=None  # Set later
        )

        interface2 = InterfaceType(
            this_mol_type_name="B",
            partner_mol_type_name="A",
            interface_index=1,
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            local_coord=np.array([-1.0, 0.0, 0.0]),
            partner_interface_type=interface1
        )

        # Set up bidirectional relationship
        interface1.partner_interface_type = interface2

        # Test bidirectional references
        self.assertEqual(interface1.partner_interface_type, interface2)
        self.assertEqual(interface2.partner_interface_type, interface1)

    def test_interface_required_free_constraints(self):
        """Test interface required_free constraints."""
        interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([0.0, 0.0, 0.0]),
            local_coord=np.array([1.0, 0.0, 0.0]),
            required_free=["A_C_1", "A_D_1", "A_E_1"]
        )

        # Test that required_free is properly stored
        self.assertEqual(len(interface.required_free), 3)
        self.assertIn("A_C_1", interface.required_free)
        self.assertIn("A_D_1", interface.required_free)
        self.assertIn("A_E_1", interface.required_free)

    def test_molecule_type_diffusion_integration(self):
        """Test integration of molecule type with diffusion calculations."""
        with patch('ionerdss.model.components.types.compute_diffusion_constants_nm_us') as mock_compute:
            mock_compute.return_value = (2.0, 1.0)

            mol_type = MoleculeType(name="TestMol", radius_nm=5.0)

            # Initially diffusion constants should be 0
            self.assertEqual(mol_type.D_t_nm2_us, 0.0)
            self.assertEqual(mol_type.D_r_rad2_us, 0.0)

            # Calculate diffusion constants
            mol_type.set_diffusion_constants_from_radius()

            # Verify they were updated
            self.assertEqual(mol_type.D_t_nm2_us, 2.0)
            self.assertEqual(mol_type.D_r_rad2_us, 1.0)
            mock_compute.assert_called_once_with(radius_nm=5.0)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions."""

    def test_interface_zero_index(self):
        """Test interface with zero index."""
        interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=0,
            absolute_coord=np.array([0.0, 0.0, 0.0]),
            local_coord=np.array([1.0, 0.0, 0.0])
        )

        self.assertEqual(interface.get_name(), "A_B_0")

    def test_interface_negative_energy(self):
        """Test interface with negative energy (favorable binding)."""
        interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([0.0, 0.0, 0.0]),
            local_coord=np.array([1.0, 0.0, 0.0]),
            energy=-100.0
        )

        self.assertEqual(interface.energy, -100.0)

    def test_interface_positive_energy(self):
        """Test interface with positive energy (unfavorable binding)."""
        interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([0.0, 0.0, 0.0]),
            local_coord=np.array([1.0, 0.0, 0.0]),
            energy=50.0
        )

        self.assertEqual(interface.energy, 50.0)

    def test_molecule_zero_radius(self):
        """Test molecule with zero radius."""
        mol_type = MoleculeType(name="ZeroRadius", radius_nm=0.0)

        self.assertEqual(mol_type.radius_nm, 0.0)

    def test_empty_signature_and_required_free(self):
        """Test interface with empty signature and required_free lists."""
        interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([0.0, 0.0, 0.0]),
            local_coord=np.array([1.0, 0.0, 0.0]),
            signature={},
            required_free=[]
        )

        self.assertEqual(interface.signature, {})
        self.assertEqual(interface.required_free, [])

    def test_large_coordinate_values(self):
        """Test interface with large coordinate values."""
        large_coord = np.array([1000.0, -1000.0, 500.0])

        interface = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=large_coord,
            local_coord=large_coord * 0.5
        )

        np.testing.assert_array_equal(interface.absolute_coord, large_coord)
        np.testing.assert_array_equal(interface.local_coord, large_coord * 0.5)


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)

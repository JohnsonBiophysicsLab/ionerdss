"""
test_reactions.py

Unit tests for the ionerdss.model.components.reactions module.

This test suite provides comprehensive coverage for the molecular binding reaction
classes, including ReactionGeometrySet and ReactionRule. Tests verify proper
initialization, serialization, BNGL expression generation, and geometric
constraint handling.

Test Classes:
    TestReactionGeometrySet: Tests for angular and geometric constraint definitions
    TestReactionRule: Tests for complete reaction rule functionality including
        BNGL expression generation, serialization, and interface handling

Test Coverage:
    - Proper initialization with default and custom values
    - BNGL reaction expression auto-generation
    - Serialization to/from dictionary format
    - Geometric constraint calculations
    - Interface name handling and validation
    - Required free interface constraint processing
    - Kinetic parameter storage and retrieval

Dependencies:
    - unittest: Standard Python testing framework
    - numpy: For coordinate array testing
    - unittest.mock: For mocking external dependencies

Example Usage:
    ```bash
    python -m pytest test_reactions.py -v
    # or
    python -m unittest test_reactions.py
    ```

See Also:
    ionerdss.model.components.reactions: Module under test
    ionerdss.model.components.types: Supporting type definitions
"""

import unittest
from unittest.mock import patch
import numpy as np

from ionerdss.model.components.reactions import ReactionGeometrySet, ReactionRule
from ionerdss.model.components.types import InterfaceType, MoleculeType


class TestReactionGeometrySet(unittest.TestCase):
    """Test suite for ReactionGeometrySet class.

    Tests the angular and geometric constraint definitions used in molecular
    binding reactions, including NERDSS angle conventions and coordinate
    transformations.
    """

    def setUp(self):
        """Set up test fixtures with sample geometric parameters."""
        # Standard test angles in radians
        self.theta1 = 0.5
        self.theta2 = 1.0
        self.phi1 = 0.25
        self.phi2 = 1.5
        self.omega = 0.75
        self.sigma_nm = 2.5

        # Sample normal vectors
        self.norm1 = np.array([1.0, 0.0, 0.0])
        self.norm2 = np.array([0.0, 1.0, 0.0])

        # Create test geometry set
        self.geometry = ReactionGeometrySet(
            theta1=self.theta1,
            theta2=self.theta2,
            phi1=self.phi1,
            phi2=self.phi2,
            omega=self.omega,
            sigma_nm=self.sigma_nm,
            norm1=self.norm1,
            norm2=self.norm2
        )

    def test_initialization(self):
        """Test proper initialization of ReactionGeometrySet."""
        # Verify all angles are set correctly
        self.assertEqual(self.geometry.theta1, self.theta1)
        self.assertEqual(self.geometry.theta2, self.theta2)
        self.assertEqual(self.geometry.phi1, self.phi1)
        self.assertEqual(self.geometry.phi2, self.phi2)
        self.assertEqual(self.geometry.omega, self.omega)
        self.assertEqual(self.geometry.sigma_nm, self.sigma_nm)

        # Verify normal vectors are stored correctly
        np.testing.assert_array_equal(self.geometry.norm1, self.norm1)
        np.testing.assert_array_equal(self.geometry.norm2, self.norm2)

    def test_as_list(self):
        """Test conversion of angles to list format."""
        angle_list = self.geometry.as_list()

        # Verify list contains all angles in correct order
        expected_list = [self.theta1, self.theta2,
                         self.phi1, self.phi2, self.omega, self.sigma_nm]
        self.assertEqual(angle_list, expected_list)
        self.assertEqual(len(angle_list), 6)

    @patch('ionerdss.model.components.reactions.compute_bond_angles_and_length_auto')
    def test_from_absolute_coords(self, mock_compute):
        """Test calculation of angles from absolute coordinates."""
        # Mock the external geometry calculation function
        mock_return_values = (0.1, 0.2, 0.3, 0.4, 0.5, 1.8,
                              np.array([0.0, 0.0, 1.0]), np.array([1.0, 1.0, 0.0]))
        mock_compute.return_value = mock_return_values

        # Sample coordinate inputs
        com1 = np.array([0.0, 0.0, 0.0])
        com2 = np.array([1.0, 1.0, 1.0])
        bind_site1 = np.array([0.5, 0.0, 0.0])
        bind_site2 = np.array([1.5, 1.0, 1.0])

        # Calculate new geometry from coordinates
        new_geometry = self.geometry.from_absolute_coords(
            com1, com2, bind_site1, bind_site2)

        # Verify the mock was called with correct arguments
        mock_compute.assert_called_once_with(
            com1, com2, bind_site1, bind_site2)

        # Verify new geometry has expected values from mock
        self.assertEqual(new_geometry.theta1, 0.1)
        self.assertEqual(new_geometry.theta2, 0.2)
        self.assertEqual(new_geometry.phi1, 0.3)
        self.assertEqual(new_geometry.phi2, 0.4)
        self.assertEqual(new_geometry.omega, 0.5)
        self.assertEqual(new_geometry.sigma_nm, 1.8)
        np.testing.assert_array_equal(
            new_geometry.norm1, np.array([0.0, 0.0, 1.0]))
        np.testing.assert_array_equal(
            new_geometry.norm2, np.array([1.0, 1.0, 0.0]))


class TestReactionRule(unittest.TestCase):
    """Test suite for ReactionRule class.

    Tests the complete reaction rule functionality including BNGL expression
    generation, interface handling, serialization, and kinetic parameters.
    """

    def setUp(self):
        """Set up test fixtures with mock molecule and interface types."""
        # Create mock molecule types
        self.mol_type_a = MoleculeType(name="ProteinA", radius_nm=2.0)
        self.mol_type_b = MoleculeType(name="ProteinB", radius_nm=1.5)

        # Create mock interface types
        self.interface_a_b = InterfaceType(
            this_mol_type_name="ProteinA",
            partner_mol_type_name="ProteinB",
            interface_index=1,
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            local_coord=np.array([1.0, 0.0, 0.0]),
            this_mol_type=self.mol_type_a,
            partner_mol_type=self.mol_type_b
        )

        self.interface_b_a = InterfaceType(
            this_mol_type_name="ProteinB",
            partner_mol_type_name="ProteinA",
            interface_index=1,
            absolute_coord=np.array([-1.0, 0.0, 0.0]),
            local_coord=np.array([-1.0, 0.0, 0.0]),
            this_mol_type=self.mol_type_b,
            partner_mol_type=self.mol_type_a
        )

        # Create sample geometry
        self.geometry = ReactionGeometrySet(
            theta1=0.0, theta2=0.0, phi1=0.0, phi2=0.0, omega=0.0,
            sigma_nm=2.5,
            norm1=np.array([1.0, 0.0, 0.0]),
            norm2=np.array([-1.0, 0.0, 0.0])
        )

    def test_initialization_basic(self):
        """Test basic initialization of ReactionRule without optional parameters."""
        reaction = ReactionRule(
            expr="initial_expr",  # This will be overwritten by update_expr
            reactant_interfaces=(self.interface_a_b, self.interface_b_a),
            ka=1e6,
            kb=1e-3
        )

        # Verify basic attributes
        self.assertEqual(reaction.reactant_interfaces,
                         (self.interface_a_b, self.interface_b_a))
        self.assertEqual(reaction.ka, 1e6)
        self.assertEqual(reaction.kb, 1e-3)
        self.assertEqual(reaction.required_free, ([], [])
                         )  # Default empty lists
        self.assertIsNone(reaction.geometry)

        # Verify BNGL expression was auto-generated
        expected_expr = "ProteinA(ProteinA_ProteinB_1) + ProteinB(ProteinB_ProteinA_1) <-> ProteinA(ProteinA_ProteinB_1!1).ProteinB(ProteinB_ProteinA_1!1)"
        self.assertEqual(reaction.expr, expected_expr)

    def test_initialization_with_geometry(self):
        """Test initialization with geometric constraints."""
        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a),
            geometry=self.geometry,
            ka=5e5,
            kb=2e-2
        )

        # Verify geometry is stored
        self.assertEqual(reaction.geometry, self.geometry)
        self.assertEqual(reaction.ka, 5e5)
        self.assertEqual(reaction.kb, 2e-2)

    def test_initialization_with_required_free(self):
        """Test initialization with required free interface constraints."""
        required_free = (["ProteinA_ProteinC_1"], [
                         "ProteinB_ProteinD_1", "ProteinB_ProteinE_1"])

        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a),
            required_free=required_free
        )

        # Verify required_free constraints are stored
        self.assertEqual(reaction.required_free, required_free)

        # Verify BNGL expression includes required free interfaces
        expected_expr = ("ProteinA(ProteinA_ProteinB_1,ProteinA_ProteinC_1) + "
                         "ProteinB(ProteinB_ProteinA_1,ProteinB_ProteinD_1,ProteinB_ProteinE_1) <-> "
                         "ProteinA(ProteinA_ProteinB_1!1,ProteinA_ProteinC_1).ProteinB(ProteinB_ProteinA_1!1,ProteinB_ProteinD_1,ProteinB_ProteinE_1)")
        self.assertEqual(reaction.expr, expected_expr)

    def test_update_expr_simple(self):
        """Test BNGL expression generation for simple binding reaction."""
        reaction = ReactionRule(
            expr="placeholder",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a)
        )

        # Manually call update_expr to test
        reaction.update_expr()

        expected_expr = "ProteinA(ProteinA_ProteinB_1) + ProteinB(ProteinB_ProteinA_1) <-> ProteinA(ProteinA_ProteinB_1!1).ProteinB(ProteinB_ProteinA_1!1)"
        self.assertEqual(reaction.expr, expected_expr)

    def test_build_molecule_expression_free(self):
        """Test building BNGL molecule expression for free state."""
        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a)
        )

        # Test free state without required free interfaces
        expr = reaction.build_molecule_expression(
            "TestMol", "TestMol_Partner_1", "free", [])
        self.assertEqual(expr, "TestMol(TestMol_Partner_1)")

        # Test free state with required free interfaces
        expr = reaction.build_molecule_expression(
            "TestMol", "TestMol_Partner_1", "free",
            ["TestMol_Other_1", "TestMol_Third_1"]
        )
        self.assertEqual(
            expr, "TestMol(TestMol_Partner_1,TestMol_Other_1,TestMol_Third_1)")

    def test_build_molecule_expression_bound(self):
        """Test building BNGL molecule expression for bound state."""
        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a)
        )

        # Test bound state
        expr = reaction.build_molecule_expression(
            "TestMol", "TestMol_Partner_1", "bound", [])
        self.assertEqual(expr, "TestMol(TestMol_Partner_1!1)")

        # Test bound state with required free interfaces
        expr = reaction.build_molecule_expression(
            "TestMol", "TestMol_Partner_1", "bound",
            ["TestMol_Other_1"]
        )
        self.assertEqual(expr, "TestMol(TestMol_Partner_1!1,TestMol_Other_1)")

    def test_build_molecule_expression_no_duplicate_interfaces(self):
        """Test that binding interface is not duplicated in required_free list."""
        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a)
        )

        # Include binding interface in required_free (should be ignored)
        expr = reaction.build_molecule_expression(
            "TestMol", "TestMol_Partner_1", "free",
            # Duplicate binding interface
            ["TestMol_Partner_1", "TestMol_Other_1"]
        )
        self.assertEqual(expr, "TestMol(TestMol_Partner_1,TestMol_Other_1)")

    def test_reactant_molecule_types_property(self):
        """Test access to reactant molecule types through property."""
        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a)
        )

        mol_types = reaction.reactant_molecule_types
        self.assertEqual(mol_types, (self.mol_type_a, self.mol_type_b))

    def test_get_reactant_interface_names(self):
        """Test retrieval of reactant interface names."""
        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a)
        )

        interface_names = reaction.get_reactant_interface_names()
        expected_names = ("ProteinA_ProteinB_1", "ProteinB_ProteinA_1")
        self.assertEqual(interface_names, expected_names)

    def test_to_dict_without_geometry(self):
        """Test serialization to dictionary without geometry."""
        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a),
            required_free=(["InterfaceA"], ["InterfaceB", "InterfaceC"]),
            ka=1e5,
            kb=1e-4
        )

        result_dict = reaction.to_dict()

        expected_dict = {
            'expr': reaction.expr,  # Auto-generated BNGL expression
            'reactant_interfaces': ["ProteinA_ProteinB_1", "ProteinB_ProteinA_1"],
            'required_free': [["InterfaceA"], ["InterfaceB", "InterfaceC"]],
            'ka': 1e5,
            'kb': 1e-4,
            'geometry': None
        }

        self.assertEqual(result_dict, expected_dict)

    def test_to_dict_with_geometry(self):
        """Test serialization to dictionary with geometry."""
        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a),
            geometry=self.geometry,
            ka=2e6,
            kb=5e-3
        )

        result_dict = reaction.to_dict()

        # Verify basic fields
        self.assertEqual(result_dict['ka'], 2e6)
        self.assertEqual(result_dict['kb'], 5e-3)
        self.assertIsNotNone(result_dict['geometry'])

        # Verify geometry fields
        geometry_dict = result_dict['geometry']
        self.assertEqual(geometry_dict['theta1'], 0.0)
        self.assertEqual(geometry_dict['theta2'], 0.0)
        self.assertEqual(geometry_dict['phi1'], 0.0)
        self.assertEqual(geometry_dict['phi2'], 0.0)
        self.assertEqual(geometry_dict['omega'], 0.0)
        self.assertEqual(geometry_dict['sigma_nm'], 2.5)

        # Verify normal vectors are converted to lists
        self.assertEqual(geometry_dict['norm1'], [1.0, 0.0, 0.0])
        self.assertEqual(geometry_dict['norm2'], [-1.0, 0.0, 0.0])

    def test_to_dict_preserves_list_types(self):
        """Test that to_dict properly converts tuples to lists for JSON compatibility."""
        required_free = (["A_C_1"], ["B_D_1"])
        reaction = ReactionRule(
            expr="test",
            reactant_interfaces=(self.interface_a_b, self.interface_b_a),
            required_free=required_free
        )

        result_dict = reaction.to_dict()

        # Verify required_free is converted to lists (not tuples)
        self.assertIsInstance(result_dict['required_free'], list)
        self.assertIsInstance(result_dict['required_free'][0], list)
        self.assertIsInstance(result_dict['required_free'][1], list)

        # Verify content is preserved
        self.assertEqual(result_dict['required_free'], [["A_C_1"], ["B_D_1"]])


if __name__ == '__main__':
    # Run the test suite
    unittest.main()

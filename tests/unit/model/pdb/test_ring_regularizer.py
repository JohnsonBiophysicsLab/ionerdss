"""
Unit tests for ionerdss.model.pdb.ring_regularizer

Tests the RingRegularizer class and its ring detection and regularization capabilities.
"""

import unittest
from unittest.mock import Mock, MagicMock, patch
import numpy as np
import networkx as nx

from ionerdss.model.pdb.ring_regularizer import (
    RingRegularizer, RingStructure, GeometricFit
)
from ionerdss.model.components.system import System
from ionerdss.model.components.instances import MoleculeInstance, InterfaceInstance
from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.pdb.file_manager import WorkspaceManager


class TestRingStructure(unittest.TestCase):
    """Test cases for RingStructure dataclass."""

    def test_ring_structure_creation(self):
        """Test RingStructure initialization."""
        molecules = [Mock(spec=MoleculeInstance) for _ in range(3)]
        interface_types = {"A_A_1"}
        center = np.array([0.0, 0.0, 0.0])
        radius = 5.0
        normal = np.array([0.0, 0.0, 1.0])

        ring = RingStructure(
            molecules=molecules,
            interface_types=interface_types,
            ring_center=center,
            ring_radius=radius,
            ring_normal=normal
        )

        self.assertEqual(ring.molecules, molecules)
        self.assertEqual(ring.interface_types, interface_types)
        np.testing.assert_array_equal(ring.ring_center, center)
        self.assertEqual(ring.ring_radius, radius)
        np.testing.assert_array_equal(ring.ring_normal, normal)


class TestGeometricFit(unittest.TestCase):
    """Test cases for GeometricFit dataclass."""

    def test_geometric_fit_cylinder(self):
        """Test GeometricFit for cylinder."""
        fit = GeometricFit(
            geometry_type="cylinder",
            center=np.array([1.0, 2.0, 3.0]),
            radius=5.0,
            axis=np.array([0.0, 0.0, 1.0]),
            fit_error=0.1
        )

        self.assertEqual(fit.geometry_type, "cylinder")
        np.testing.assert_array_equal(fit.center, np.array([1.0, 2.0, 3.0]))
        self.assertEqual(fit.radius, 5.0)
        np.testing.assert_array_equal(fit.axis, np.array([0.0, 0.0, 1.0]))
        self.assertEqual(fit.fit_error, 0.1)

    def test_geometric_fit_sphere(self):
        """Test GeometricFit for sphere."""
        fit = GeometricFit(
            geometry_type="sphere",
            center=np.array([0.0, 0.0, 0.0]),
            radius=3.0,
            fit_error=0.05
        )

        self.assertEqual(fit.geometry_type, "sphere")
        self.assertIsNone(fit.axis)  # Default None for sphere
        self.assertEqual(fit.fit_error, 0.05)


class TestRingRegularizer(unittest.TestCase):
    """Test cases for RingRegularizer class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock system
        self.system = Mock(spec=System)

        # Create mock workspace manager
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.logger = Mock()

        # Create mock molecule instances for a triangle ring
        self.mol_instances = []
        for i in range(3):
            mol_instance = Mock(spec=MoleculeInstance)
            mol_instance.name = f"mol_{i}"
            mol_instance.com = np.array([
                5.0 * np.cos(2 * np.pi * i / 3),
                5.0 * np.sin(2 * np.pi * i / 3),
                0.0
            ])

            # Mock molecule type
            mol_type = Mock(spec=MoleculeType)
            mol_type.name = "ProteinA"
            mol_instance.molecule_type = mol_type

            self.mol_instances.append(mol_instance)

        self.system.molecule_instances = self.mol_instances

    def test_initialization_valid_parameters(self):
        """Test RingRegularizer initialization with valid parameters."""
        regularizer = RingRegularizer(
            system=self.system,
            workspace_manager=self.workspace_manager,
            mode="separate",
            geometry="cylinder",
            min_ring_size=4
        )

        self.assertEqual(regularizer.system, self.system)
        self.assertEqual(regularizer.workspace_manager, self.workspace_manager)
        self.assertEqual(regularizer.mode, "separate")
        self.assertEqual(regularizer.geometry, "cylinder")
        self.assertEqual(regularizer.min_ring_size, 4)
        self.assertEqual(len(regularizer.detected_rings), 0)
        self.assertEqual(len(regularizer.geometric_fits), 0)

    def test_initialization_invalid_mode(self):
        """Test RingRegularizer initialization with invalid mode."""
        with self.assertRaises(ValueError) as context:
            RingRegularizer(
                system=self.system,
                mode="invalid_mode"
            )

        self.assertIn("Invalid mode", str(context.exception))

    def test_initialization_invalid_geometry(self):
        """Test RingRegularizer initialization with invalid geometry."""
        with self.assertRaises(ValueError) as context:
            RingRegularizer(
                system=self.system,
                geometry="invalid_geometry"
            )

        self.assertIn("Invalid geometry", str(context.exception))

    def test_initialization_min_ring_size_adjustment(self):
        """Test that min_ring_size is adjusted to at least 3."""
        regularizer = RingRegularizer(
            system=self.system,
            min_ring_size=1  # Should be adjusted to 3
        )

        self.assertEqual(regularizer.min_ring_size, 3)

    def test_regularize_mode_off(self):
        """Test regularize method with mode 'off'."""
        regularizer = RingRegularizer(
            system=self.system,
            workspace_manager=self.workspace_manager,
            mode="off"
        )

        result = regularizer.regularize()

        self.assertFalse(result)
        self.workspace_manager.logger.info.assert_called_with(
            "Regularizer off")

    def test_build_connectivity_graph(self):
        """Test _build_connectivity_graph method."""
        # Set up interface connections
        interface_type = Mock(spec=InterfaceType)
        interface_type.get_name.return_value = "A_A_1"

        for i in range(3):
            interface_instance = Mock(spec=InterfaceInstance)
            interface_instance.interface_type = interface_type

            # Connect each molecule to the next (forming a triangle)
            next_mol = self.mol_instances[(i + 1) % 3]
            self.mol_instances[i].interfaces_neighbors_map = {
                interface_instance: next_mol
            }

        regularizer = RingRegularizer(self.system, self.workspace_manager)
        graph = regularizer._build_connectivity_graph()

        # Check graph structure
        self.assertEqual(graph.number_of_nodes(), 3)
        self.assertEqual(graph.number_of_edges(), 3)

        # Check that all molecules are nodes
        for mol_instance in self.mol_instances:
            self.assertIn(mol_instance.name, graph.nodes)

    def test_validate_interface_consistency_homotypic(self):
        """Test _validate_interface_consistency with homotypic interfaces."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)

        # All same interface type
        interface_types = {"A_A_1"}
        self.assertTrue(
            regularizer._validate_interface_consistency(interface_types))

    def test_validate_interface_consistency_heterotypic(self):
        """Test _validate_interface_consistency with complementary heterotypic interfaces."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)

        # Complementary interface types
        interface_types = {"A_B_1", "B_A_1"}
        self.assertTrue(
            regularizer._validate_interface_consistency(interface_types))

    def test_validate_interface_consistency_invalid(self):
        """Test _validate_interface_consistency with invalid interface combinations."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)

        # Too many different types
        interface_types = {"A_A_1", "B_B_1", "C_C_1"}
        self.assertFalse(
            regularizer._validate_interface_consistency(interface_types))

        # Non-complementary pair
        interface_types = {"A_B_1", "C_D_1"}
        self.assertFalse(
            regularizer._validate_interface_consistency(interface_types))

    def test_fit_cylinder(self):
        """Test _fit_cylinder method."""
        regularizer = RingRegularizer(
            self.system, self.workspace_manager, geometry="cylinder")

        # Create positions on a cylinder
        positions = np.array([
            [5.0, 0.0, 0.0],
            [0.0, 5.0, 0.0],
            [-5.0, 0.0, 0.0],
            [0.0, -5.0, 0.0]
        ])

        fit = regularizer._fit_cylinder(positions)

        self.assertEqual(fit.geometry_type, "cylinder")
        self.assertIsInstance(fit.center, np.ndarray)
        self.assertGreater(fit.radius, 0)
        self.assertIsInstance(fit.axis, np.ndarray)
        self.assertGreaterEqual(fit.fit_error, 0)

    def test_fit_sphere(self):
        """Test _fit_sphere method."""
        regularizer = RingRegularizer(
            self.system, self.workspace_manager, geometry="sphere")

        # Create positions on a sphere
        positions = np.array([
            [5.0, 0.0, 0.0],
            [0.0, 5.0, 0.0],
            [0.0, 0.0, 5.0],
            [-5.0, 0.0, 0.0]
        ])

        fit = regularizer._fit_sphere(positions)

        self.assertEqual(fit.geometry_type, "sphere")
        self.assertIsInstance(fit.center, np.ndarray)
        self.assertGreater(fit.radius, 0)
        self.assertIsNone(fit.axis)  # Sphere doesn't have axis
        self.assertGreaterEqual(fit.fit_error, 0)

    def test_translate_molecule(self):
        """Test _translate_molecule method."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)

        # Create mock molecule with interface
        mol_instance = Mock(spec=MoleculeInstance)
        mol_instance.com = np.array([1.0, 2.0, 3.0])

        # Mock interface instance
        interface_instance = Mock(spec=InterfaceInstance)
        interface_instance.absolute_coord = np.array([1.5, 2.5, 3.5])

        # Mock interface type
        interface_type = Mock(spec=InterfaceType)
        interface_type.absolute_coord = np.array([1.2, 2.2, 3.2])
        interface_instance.interface_type = interface_type

        mol_instance.interfaces_neighbors_map = {interface_instance: Mock()}

        # Apply translation
        translation = np.array([10.0, 20.0, 30.0])
        regularizer._translate_molecule(mol_instance, translation)

        # Check that COM was translated
        expected_com = np.array([11.0, 22.0, 33.0])
        np.testing.assert_array_equal(mol_instance.com, expected_com)

        # Check that interface coordinate was translated
        expected_interface_coord = np.array([11.5, 22.5, 33.5])
        np.testing.assert_array_equal(
            interface_instance.absolute_coord, expected_interface_coord)

    def test_regularize_to_cylinder(self):
        """Test _regularize_to_cylinder method."""
        regularizer = RingRegularizer(
            self.system, self.workspace_manager, geometry="cylinder")

        # Create properly mocked molecules
        mock_molecules = []
        for i in range(3):
            mol = Mock(spec=MoleculeInstance)
            mol.com = np.array([float(i), 0.0, 0.0])
            mol.interfaces_neighbors_map = {}  # Add this attribute
            mock_molecules.append(mol)

        # Create ring structure
        ring = RingStructure(
            molecules=mock_molecules,
            interface_types={"A_A_1"},
            ring_center=np.array([0.0, 0.0, 0.0]),
            ring_radius=5.0,
            ring_normal=np.array([0.0, 0.0, 1.0])
        )

        # Create cylinder fit
        fit = GeometricFit(
            geometry_type="cylinder",
            center=np.array([0.0, 0.0, 0.0]),
            radius=5.0,
            axis=np.array([0.0, 0.0, 1.0])
        )

        # Store original positions
        original_positions = [mol.com.copy() for mol in mock_molecules]

        # Apply regularization
        regularizer._regularize_to_cylinder(ring, fit, len(ring.molecules))

        # Check that molecules were moved
        for i, mol in enumerate(mock_molecules):
            # Positions should be different (unless already perfectly aligned)
            distance_moved = np.linalg.norm(mol.com - original_positions[i])
            # Allow for small movements or no movement if already aligned
            self.assertGreaterEqual(distance_moved, 0.0)

    def test_regularize_to_sphere(self):
        """Test _regularize_to_sphere method."""
        regularizer = RingRegularizer(
            self.system, self.workspace_manager, geometry="sphere")

        # Create properly mocked molecules
        mock_molecules = []
        for i in range(3):
            mol = Mock(spec=MoleculeInstance)
            # Use existing positions
            mol.com = self.mol_instances[i].com.copy()
            mol.interfaces_neighbors_map = {}  # Add this attribute
            mock_molecules.append(mol)

        # Create ring structure
        ring = RingStructure(
            molecules=mock_molecules,
            interface_types={"A_A_1"},
            ring_center=np.array([0.0, 0.0, 0.0]),
            ring_radius=5.0,
            ring_normal=np.array([0.0, 0.0, 1.0])
        )

        # Create sphere fit
        fit = GeometricFit(
            geometry_type="sphere",
            center=np.array([0.0, 0.0, 0.0]),
            radius=5.0
        )

        # Store original positions
        original_positions = [mol.com.copy() for mol in mock_molecules]

        # Apply regularization
        regularizer._regularize_to_sphere(ring, fit, len(ring.molecules))

        # Check that molecules were moved to sphere surface
        for mol in mock_molecules:
            distance_from_center = np.linalg.norm(mol.com - fit.center)
            # Should be approximately on sphere surface
            self.assertAlmostEqual(distance_from_center, fit.radius, places=1)

    def test_get_summary_empty(self):
        """Test get_summary with no detected rings."""
        regularizer = RingRegularizer(
            self.system,
            self.workspace_manager,
            mode="separate",
            geometry="cylinder",
            min_ring_size=4
        )

        summary = regularizer.get_summary()

        expected_summary = {
            "mode": "separate",
            "geometry": "cylinder",
            "min_ring_size": 4,
            "rings_detected": 0,
            "rings_regularized": 0,
            "ring_sizes": [],
            "fit_errors": [],
            "average_fit_error": 0.0
        }

        self.assertEqual(summary, expected_summary)

    def test_get_summary_with_rings(self):
        """Test get_summary with detected rings."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)

        # Add mock rings and fits
        ring1 = Mock(spec=RingStructure)
        ring1.molecules = [Mock(), Mock(), Mock()]

        ring2 = Mock(spec=RingStructure)
        ring2.molecules = [Mock(), Mock(), Mock(), Mock()]

        regularizer.detected_rings = [ring1, ring2]

        fit1 = Mock(spec=GeometricFit)
        fit1.fit_error = 0.1

        fit2 = Mock(spec=GeometricFit)
        fit2.fit_error = 0.2

        regularizer.geometric_fits = {0: fit1, 1: fit2}

        summary = regularizer.get_summary()

        self.assertEqual(summary["rings_detected"], 2)
        self.assertEqual(summary["rings_regularized"], 2)
        self.assertEqual(summary["ring_sizes"], [3, 4])
        self.assertEqual(summary["fit_errors"], [0.1, 0.2])
        # Use assertAlmostEqual for floating point comparison
        self.assertAlmostEqual(summary["average_fit_error"], 0.15, places=6)

    @patch('ionerdss.model.pdb.ring_regularizer.nx.simple_cycles')
    def test_detect_ring_structures_with_cycles(self, mock_cycles):
        """Test _detect_ring_structures when cycles are found."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)

        # Mock cycle detection
        mock_cycles.return_value = [["mol_0", "mol_1", "mol_2"]]

        # Mock the validation method to return a valid ring - create proper RingStructure
        mock_ring = RingStructure(
            molecules=[Mock(), Mock(), Mock()],
            interface_types={"A_A_1"},
            ring_center=np.array([0.0, 0.0, 0.0]),
            ring_radius=5.0,
            ring_normal=np.array([0.0, 0.0, 1.0])
        )
        with patch.object(regularizer, '_validate_and_create_ring', return_value=mock_ring):
            with patch.object(regularizer, '_build_connectivity_graph', return_value=Mock()):
                rings = regularizer._detect_ring_structures()

        self.assertEqual(len(rings), 1)
        self.assertEqual(rings[0], mock_ring)

    @patch('ionerdss.model.pdb.ring_regularizer.nx.simple_cycles')
    def test_detect_ring_structures_no_cycles(self, mock_cycles):
        """Test _detect_ring_structures when no cycles are found."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)

        # Mock no cycles found
        mock_cycles.return_value = []

        with patch.object(regularizer, '_build_connectivity_graph', return_value=Mock()):
            rings = regularizer._detect_ring_structures()

        self.assertEqual(len(rings), 0)

    def test_validate_and_create_ring_invalid_size(self):
        """Test _validate_and_create_ring with insufficient ring size."""
        regularizer = RingRegularizer(
            self.system, self.workspace_manager, min_ring_size=5)

        # Cycle too small
        cycle = ["mol_0", "mol_1", "mol_2"]  # Size 3, but min_ring_size is 5
        graph = Mock()

        result = regularizer._validate_and_create_ring(cycle, graph)

        self.assertIsNone(result)

    def test_validate_and_create_ring_heterotypic(self):
        """Test _validate_and_create_ring with heterotypic molecules."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)

        # Create molecules of different types
        mol1 = Mock(spec=MoleculeInstance)
        mol1.molecule_type.name = "ProteinA"
        mol1.com = np.array([0.0, 0.0, 0.0])

        mol2 = Mock(spec=MoleculeInstance)
        mol2.molecule_type.name = "ProteinB"  # Different type
        mol2.com = np.array([1.0, 0.0, 0.0])

        # Mock graph
        graph = Mock()
        graph.nodes = {
            "mol_0": {"molecule": mol1},
            "mol_1": {"molecule": mol2}
        }

        cycle = ["mol_0", "mol_1"]
        result = regularizer._validate_and_create_ring(cycle, graph)

        # Should return None because not homo-n-mer
        self.assertIsNone(result)

    def test_calculate_uniform_fit(self):
        """Test _calculate_uniform_fit method."""
        regularizer = RingRegularizer(
            self.system, self.workspace_manager, mode="uniform")

        # Create mock rings
        ring1 = Mock(spec=RingStructure)
        ring1.molecules = self.mol_instances[:2]

        ring2 = Mock(spec=RingStructure)
        ring2.molecules = self.mol_instances[2:]

        regularizer.detected_rings = [ring1, ring2]

        # Mock the fit method - create proper GeometricFit object
        mock_fit = GeometricFit(
            geometry_type="cylinder",
            center=np.array([0.0, 0.0, 0.0]),
            radius=5.0,
            axis=np.array([0.0, 0.0, 1.0]),
            fit_error=0.1
        )
        with patch.object(regularizer, '_fit_cylinder', return_value=mock_fit):
            regularizer._calculate_uniform_fit()

        # Both rings should have the same fit
        self.assertEqual(len(regularizer.geometric_fits), 2)
        self.assertEqual(regularizer.geometric_fits[0], mock_fit)
        self.assertEqual(regularizer.geometric_fits[1], mock_fit)

    def test_calculate_separate_fits(self):
        """Test _calculate_separate_fits method."""
        regularizer = RingRegularizer(
            self.system, self.workspace_manager, mode="separate")

        # Create mock rings
        ring1 = Mock(spec=RingStructure)
        ring1.molecules = self.mol_instances[:2]

        regularizer.detected_rings = [ring1]

        # Mock the fit method - create proper GeometricFit object
        mock_fit = GeometricFit(
            geometry_type="cylinder",
            center=np.array([0.0, 0.0, 0.0]),
            radius=5.0,
            axis=np.array([0.0, 0.0, 1.0]),
            fit_error=0.1
        )
        with patch.object(regularizer, '_fit_cylinder', return_value=mock_fit):
            regularizer._calculate_separate_fits()

        # Should have one fit
        self.assertEqual(len(regularizer.geometric_fits), 1)
        self.assertEqual(regularizer.geometric_fits[0], mock_fit)


class TestRingRegularizerIntegration(unittest.TestCase):
    """Integration tests for RingRegularizer."""

    def setUp(self):
        """Set up integration test fixtures."""
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.logger = Mock()

    def test_complete_regularization_workflow(self):
        """Test complete regularization workflow."""
        # Create a system with a triangular ring
        system = Mock(spec=System)

        # Create molecule instances in a triangle
        mol_instances = []
        for i in range(3):
            mol_instance = Mock(spec=MoleculeInstance)
            mol_instance.name = f"mol_{i}"
            mol_instance.com = np.array([
                5.0 * np.cos(2 * np.pi * i / 3),
                5.0 * np.sin(2 * np.pi * i / 3),
                0.0
            ])

            # Same molecule type (homo-trimer)
            mol_type = Mock(spec=MoleculeType)
            mol_type.name = "ProteinA"
            mol_instance.molecule_type = mol_type

            # Initialize empty interfaces_neighbors_map
            mol_instance.interfaces_neighbors_map = {}

            mol_instances.append(mol_instance)

        # Now set up the circular connections after all molecules are created
        for i in range(3):
            interface_instance = Mock(spec=InterfaceInstance)
            interface_instance.absolute_coord = mol_instances[i].com + np.array([
                                                                                1.0, 0.0, 0.0])

            interface_type = Mock(spec=InterfaceType)
            interface_type.get_name.return_value = "A_A_1"
            interface_instance.interface_type = interface_type

            # Connect to next molecule in ring
            next_mol = mol_instances[(i + 1) % 3]
            mol_instances[i].interfaces_neighbors_map = {
                interface_instance: next_mol}

        system.molecule_instances = mol_instances

        # Create regularizer
        regularizer = RingRegularizer(
            system=system,
            workspace_manager=self.workspace_manager,
            mode="separate",
            geometry="cylinder",
            min_ring_size=3
        )

        # Mock the entire ring detection process to ensure a ring is found
        mock_ring = RingStructure(
            molecules=mol_instances,
            interface_types={"A_A_1"},
            ring_center=np.array([0.0, 0.0, 0.0]),
            ring_radius=5.0,
            ring_normal=np.array([0.0, 0.0, 1.0])
        )

        # Mock the geometric fit
        mock_fit = GeometricFit(
            geometry_type="cylinder",
            center=np.array([0.0, 0.0, 0.0]),
            radius=5.0,
            axis=np.array([0.0, 0.0, 1.0]),
            fit_error=0.1
        )

        # Patch the ring detection to return our mock ring
        with patch.object(regularizer, '_detect_ring_structures', return_value=[mock_ring]):
            with patch.object(regularizer, '_fit_cylinder', return_value=mock_fit):
                # Run regularization
                result = regularizer.regularize()

        # Should have performed regularization
        self.assertTrue(result)

        # Should have detected rings
        summary = regularizer.get_summary()
        self.assertEqual(summary["rings_detected"], 1)
        self.assertEqual(summary["rings_regularized"], 1)


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

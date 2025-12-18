"""
Unit tests for ionerdss.model.pdb.template_builder

Tests the TemplateBuilder class and its template generation capabilities.
"""

import unittest
from unittest.mock import Mock, patch
import numpy as np

from ionerdss.model.pdb.template_builder import (
    TemplateBuilder, GeometricSignature
)
from ionerdss.model.pdb.parser import PDBParser
from ionerdss.model.pdb.coarse_graining import CoarseGrainer, CoarseGrainedChain, InterfaceString
from ionerdss.model.pdb.chain_grouping import ChainGrouper, ChainGroup
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.file_manager import WorkspaceManager
from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.components.units import Units


class TestGeometricSignature(unittest.TestCase):
    """Test cases for GeometricSignature dataclass."""

    def test_geometric_signature_creation(self):
        """Test GeometricSignature initialization."""
        signature = GeometricSignature(
            d_i=5.0,
            d_j=6.0,
            theta_i=1.2,
            theta_j=1.5
        )

        self.assertEqual(signature.d_i, 5.0)
        self.assertEqual(signature.d_j, 6.0)
        self.assertEqual(signature.theta_i, 1.2)
        self.assertEqual(signature.theta_j, 1.5)

    def test_normalize(self):
        """Test signature normalization."""
        signature = GeometricSignature(
            d_i=5.123456789,
            d_j=6.987654321,
            theta_i=1.234567,
            theta_j=1.567890
        )

        normalized = signature.normalize(precision=3)
        expected = (5.123, 6.988, 1.235, 1.568)

        self.assertEqual(normalized, expected)

    def test_is_similar_to(self):
        """Test signature similarity comparison."""
        signature1 = GeometricSignature(5.0, 6.0, 1.2, 1.5)
        signature2 = GeometricSignature(5.1, 6.1, 1.25, 1.55)
        signature3 = GeometricSignature(7.0, 8.0, 2.0, 2.5)

        # Similar signatures
        self.assertTrue(signature1.is_similar_to(
            signature2, distance_threshold=0.2, angle_threshold=0.1
        ))

        # Dissimilar signatures
        self.assertFalse(signature1.is_similar_to(
            signature3, distance_threshold=0.2, angle_threshold=0.1
        ))

    def test_is_homotypic(self):
        """Test homotypic signature detection."""
        # Homotypic signature (symmetric)
        homotypic = GeometricSignature(5.0, 5.1, 1.2, 1.25)
        self.assertTrue(homotypic.is_homotypic(
            distance_threshold=0.2, angle_threshold=0.1
        ))

        # Heterotypic signature (asymmetric)
        heterotypic = GeometricSignature(5.0, 7.0, 1.2, 2.0)
        self.assertFalse(heterotypic.is_homotypic(
            distance_threshold=0.2, angle_threshold=0.1
        ))


class TestTemplateBuilder(unittest.TestCase):
    """Test cases for TemplateBuilder class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock components
        self.parser = Mock(spec=PDBParser)
        self.coarse_grainer = Mock(spec=CoarseGrainer)
        self.chain_grouper = Mock(spec=ChainGrouper)
        self.hyperparams = Mock(spec=PDBModelHyperparameters)
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.logger = Mock()

        # Set up mock data
        self._setup_mock_data()

    def _setup_mock_data(self):
        """Set up mock data for testing."""
        # Mock hyperparameters
        self.hyperparams.homodimer_distance_threshold = 1.0
        self.hyperparams.homodimer_angle_threshold = 0.2
        self.hyperparams.signature_precision = 6
        self.hyperparams.steric_clash_mode = "off"

        # Mock coarse-grained chains
        chain_a = Mock(spec=CoarseGrainedChain)
        chain_a.chain_id = "A"
        chain_a.com = np.array([0.0, 0.0, 0.0])
        chain_a.radius = 15.0  # Angstroms

        chain_b = Mock(spec=CoarseGrainedChain)
        chain_b.chain_id = "B"
        chain_b.com = np.array([10.0, 0.0, 0.0])
        chain_b.radius = 12.0  # Angstroms

        self.coarse_grainer.get_coarse_grained_chains.return_value = {
            "A": chain_a,
            "B": chain_b
        }

        # Mock interfaces
        interface_ab = Mock(spec=InterfaceString)
        interface_ab.chain_i = "A"
        interface_ab.chain_j = "B"
        interface_ab.coord_i = np.array([5.0, 0.0, 0.0])
        interface_ab.coord_j = np.array([5.0, 0.0, 0.0])
        interface_ab.residues_i = {1, 2, 3}
        interface_ab.residues_j = {4, 5, 6}
        interface_ab.energy = -5.0

        self.coarse_grainer.get_interfaces.return_value = [interface_ab]

        # Mock chain groups
        group_a = Mock(spec=ChainGroup)
        group_a.representative = "A"
        group_a.members = ["A"]
        group_a.grouping_method = "sequence_similarity"

        group_b = Mock(spec=ChainGroup)
        group_b.representative = "B"
        group_b.members = ["B"]
        group_b.grouping_method = "sequence_similarity"

        self.chain_grouper.get_groups.return_value = [group_a, group_b]
        self.chain_grouper.get_group_for_chain.side_effect = lambda chain_id: {
            "A": group_a,
            "B": group_b
        }.get(chain_id)

        # Mock parser methods
        self.parser.convert_coords_to_nm.side_effect = lambda coords: coords / 10.0
        self.parser.get_chain_data.return_value = {
            'ca_coords': np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
        }

    def test_initialization(self):
        """Test TemplateBuilder initialization."""
        with patch.object(TemplateBuilder, '_build_templates'):
            with patch.object(TemplateBuilder, '_regularize_templates'):
                with patch.object(TemplateBuilder, '_detect_steric_clashes'):
                    builder = TemplateBuilder(
                        parser=self.parser,
                        coarse_grainer=self.coarse_grainer,
                        chain_grouper=self.chain_grouper,
                        hyperparams=self.hyperparams,
                        workspace_manager=self.workspace_manager
                    )

        self.assertEqual(builder.parser, self.parser)
        self.assertEqual(builder.coarse_grainer, self.coarse_grainer)
        self.assertEqual(builder.chain_grouper, self.chain_grouper)
        self.assertEqual(builder.hyperparams, self.hyperparams)
        self.assertEqual(builder.workspace_manager, self.workspace_manager)
        self.assertIsInstance(builder.units, Units)

    def test_initialization_with_custom_units(self):
        """Test TemplateBuilder initialization with custom units."""
        custom_units = Mock(spec=Units)

        with patch.object(TemplateBuilder, '_build_templates'):
            with patch.object(TemplateBuilder, '_regularize_templates'):
                with patch.object(TemplateBuilder, '_detect_steric_clashes'):
                    builder = TemplateBuilder(
                        parser=self.parser,
                        coarse_grainer=self.coarse_grainer,
                        chain_grouper=self.chain_grouper,
                        hyperparams=self.hyperparams,
                        units=custom_units,
                        workspace_manager=self.workspace_manager
                    )

        self.assertEqual(builder.units, custom_units)

    def test_generate_template_name_unique(self):
        """Test template name generation for unique representative."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.used_template_names = set()

        group = Mock(spec=ChainGroup)
        group.representative = "A"
        group.members = ["A"]

        name = builder._generate_template_name(group)

        self.assertEqual(name, "A")
        self.assertIn("A", builder.used_template_names)

    def test_generate_template_name_conflict(self):
        """Test template name generation with conflicts."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.used_template_names = {"A"}  # A is already used

        group = Mock(spec=ChainGroup)
        group.representative = "A"
        group.members = ["A", "C"]  # Multiple members

        name = builder._generate_template_name(group)

        self.assertEqual(name, "A0")
        self.assertIn("A0", builder.used_template_names)

    def test_generate_template_name_numeric_suffix(self):
        """Test template name generation with numeric suffix."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.used_template_names = {"A", "A_group"}  # Both variants used

        group = Mock(spec=ChainGroup)
        group.representative = "A"
        group.members = ["A"]

        name = builder._generate_template_name(group)

        self.assertEqual(name, "A1")
        self.assertIn("A1", builder.used_template_names)

    def test_build_molecule_template(self):
        """Test molecule template building."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.coarse_grainer = self.coarse_grainer
        builder.parser = self.parser
        builder.workspace_manager = self.workspace_manager
        builder.used_template_names = set()
        builder.group_to_template = {}
        builder.molecule_templates = {}

        group = Mock(spec=ChainGroup)
        group.representative = "A"
        group.members = ["A"]
        group.grouping_method = "sequence_similarity"

        builder._build_molecule_template(group)

        # Check that template was created
        self.assertIn("A", builder.molecule_templates)
        self.assertEqual(builder.group_to_template["A"], "A")

        # Check template properties
        template = builder.molecule_templates["A"]
        self.assertIsInstance(template, MoleculeType)
        self.assertEqual(template.name, "A")
        self.assertEqual(template.radius_nm, 1.5)  # 15.0 Å / 10

    def test_calculate_geometric_signature(self):
        """Test geometric signature calculation."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.coarse_grainer = self.coarse_grainer

        interface = Mock(spec=InterfaceString)
        interface.chain_i = "A"
        interface.chain_j = "B"
        interface.coord_i = np.array([5.0, 0.0, 0.0])
        interface.coord_j = np.array([5.0, 0.0, 0.0])

        signature = builder._calculate_geometric_signature(interface)

        self.assertIsInstance(signature, GeometricSignature)
        self.assertGreater(signature.d_i, 0)
        self.assertGreater(signature.d_j, 0)
        self.assertGreaterEqual(signature.theta_i, 0)
        self.assertGreaterEqual(signature.theta_j, 0)

    def test_create_homotypic_interface_template(self):
        """Test homodimer interface template creation."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.coarse_grainer = self.coarse_grainer
        builder.parser = self.parser
        builder.workspace_manager = self.workspace_manager
        builder.hyperparams = self.hyperparams  # Add this line
        builder.interface_templates = {}
        builder.interface_signatures = {}
        builder.molecule_templates = {
            "A": Mock(spec=MoleculeType, interfaces_neighbors_map={})
        }

        interface = Mock(spec=InterfaceString)
        interface.chain_i = "A"
        interface.chain_j = "A"
        interface.coord_i = np.array([5.0, 0.0, 0.0])
        interface.residues_i = {1, 2}
        interface.residues_j = {3, 4}
        interface.energy = -3.0

        signature = GeometricSignature(5.0, 5.0, 1.0, 1.0)

        interface_name = builder._create_homotypic_interface_template(
            interface, "A", signature, 1
        )

        # Check that template was created
        expected_name = "A_A_1"
        self.assertEqual(interface_name, expected_name)
        self.assertIn(expected_name, builder.interface_templates)
        self.assertIn(expected_name, builder.interface_signatures)

        # Check template properties
        template = builder.interface_templates[expected_name]
        self.assertEqual(template.this_mol_type_name, "A")
        self.assertEqual(template.partner_mol_type_name, "A")
        self.assertEqual(template.interface_index, 1)

    def test_create_heterotypic_interface_templates(self):
        """Test heterotypic interface template creation."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.coarse_grainer = self.coarse_grainer
        builder.parser = self.parser
        builder.workspace_manager = self.workspace_manager
        builder.hyperparams = self.hyperparams  # Add this line
        builder.interface_templates = {}
        builder.interface_signatures = {}
        builder.molecule_templates = {
            "A": Mock(spec=MoleculeType, interfaces_neighbors_map={}),
            "B": Mock(spec=MoleculeType, interfaces_neighbors_map={})
        }

        interface = Mock(spec=InterfaceString)
        interface.chain_i = "A"
        interface.chain_j = "B"
        interface.coord_i = np.array([5.0, 0.0, 0.0])
        interface.coord_j = np.array([5.0, 0.0, 0.0])
        interface.residues_i = {1, 2}
        interface.residues_j = {3, 4}
        interface.energy = -3.0

        signature = GeometricSignature(5.0, 5.0, 1.0, 1.2)

        interface_names = builder._create_heterotypic_interface_templates(
            interface, "A", "B", signature, 1
        )

        # Check that both templates were created
        expected_names = ["A_B_1", "B_A_1"]
        self.assertEqual(len(interface_names), 2)
        self.assertCountEqual(interface_names, expected_names)

        # Check templates exist
        for name in expected_names:
            self.assertIn(name, builder.interface_templates)
            self.assertIn(name, builder.interface_signatures)

        # Check cross-references
        template_a = builder.interface_templates["A_B_1"]
        template_b = builder.interface_templates["B_A_1"]
        self.assertEqual(template_a.partner_interface_type, template_b)
        self.assertEqual(template_b.partner_interface_type, template_a)

    def test_find_matching_interface_type_found(self):
        """Test finding matching interface type."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.interface_templates = {}
        builder.interface_signatures = {}
        builder.workspace_manager = self.workspace_manager

        # Create existing interface template
        existing_template = Mock(spec=InterfaceType)
        existing_template.this_mol_type_name = "A"
        existing_template.partner_mol_type_name = "B"

        existing_signature = GeometricSignature(5.0, 6.0, 1.0, 1.2)

        builder.interface_templates["A_B_1"] = existing_template
        builder.interface_signatures["A_B_1"] = existing_signature

        # Test with similar signature
        test_signature = GeometricSignature(5.1, 6.1, 1.05, 1.25)

        match = builder._find_matching_interface_type("A", "B", test_signature)

        self.assertEqual(match, "A_B_1")

    def test_find_matching_interface_type_not_found(self):
        """Test finding matching interface type when none exists."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.interface_templates = {}
        builder.interface_signatures = {}
        builder.workspace_manager = self.workspace_manager

        test_signature = GeometricSignature(5.0, 6.0, 1.0, 1.2)

        match = builder._find_matching_interface_type("A", "B", test_signature)

        self.assertIsNone(match)

    def test_process_interface_with_signature_existing(self):
        """Test processing interface with existing signature match."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.workspace_manager = self.workspace_manager

        # Mock finding existing match
        with patch.object(builder, '_find_matching_interface_type', return_value="A_B_1"):
            interface = Mock(spec=InterfaceString)
            interface.chain_i = "A"  # Add this line
            interface.chain_j = "B"  # Add this line
            signature = GeometricSignature(5.0, 6.0, 1.0, 1.2)

            builder._process_interface_with_signature(
                interface, "A", "B", signature
            )

            # Check that interface was assigned existing type
            self.assertEqual(interface.interface_type, "A_B_1")

    def test_process_interface_with_signature_new(self):
        """Test processing interface with new signature."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.workspace_manager = self.workspace_manager

        # Mock no existing match and new type creation
        with patch.object(builder, '_find_matching_interface_type', return_value=None):
            with patch.object(builder, '_create_new_interface_type', return_value=["A_B_1"]):
                interface = Mock(spec=InterfaceString)
                signature = GeometricSignature(5.0, 6.0, 1.0, 1.2)

                builder._process_interface_with_signature(
                    interface, "A", "B", signature
                )

                # Check that interface was assigned new type
                self.assertEqual(interface.interface_type, "A_B_1")

    def test_assign_interface_to_heterotypic_type(self):
        """Test assigning interface to heterotypic type."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.interface_templates = {}

        # Create mock interface templates
        template1 = Mock(spec=InterfaceType)
        template1.this_mol_type_name = "A"
        template1.signature = {"this_side": "i"}

        template2 = Mock(spec=InterfaceType)
        template2.this_mol_type_name = "B"
        template2.signature = {"this_side": "j"}

        builder.interface_templates["A_B_1"] = template1
        builder.interface_templates["B_A_1"] = template2

        interface = Mock(spec=InterfaceString)
        interface_type_names = ["A_B_1", "B_A_1"]

        assigned_type = builder._assign_interface_to_heterotypic_type(
            interface, interface_type_names, "A", "B"
        )

        self.assertEqual(assigned_type, "A_B_1")

    def test_get_interface_type_for_interface(self):
        """Test getting interface type for interface."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.interface_to_type_mapping = {
            "A_B_5.000_0.000_0.000": "A_B_1"
        }

        interface = Mock(spec=InterfaceString)
        interface.chain_i = "A"
        interface.chain_j = "B"
        interface.coord_i = np.array([5.0, 0.0, 0.0])

        interface_type = builder.get_interface_type_for_interface(interface)

        self.assertEqual(interface_type, "A_B_1")

    def test_get_interface_type_for_interface_not_found(self):
        """Test getting interface type for interface when not found."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.interface_to_type_mapping = {}

        interface = Mock(spec=InterfaceString)
        interface.chain_i = "A"
        interface.chain_j = "B"
        interface.coord_i = np.array([5.0, 0.0, 0.0])

        interface_type = builder.get_interface_type_for_interface(interface)

        self.assertIsNone(interface_type)

    def test_compute_rigid_transform(self):
        """Test rigid transform computation."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.parser = self.parser

        ref_data = Mock(spec=CoarseGrainedChain)
        ref_data.chain_id = "A"

        mem_data = Mock(spec=CoarseGrainedChain)
        mem_data.chain_id = "B"

        # Mock parser to return same coordinates (identity transform expected)
        self.parser.get_chain_data.return_value = {
            'ca_coords': np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        }

        transform = builder._compute_rigid_transform(ref_data, mem_data)

        # Should return 4x4 transformation matrix
        self.assertEqual(transform.shape, (4, 4))
        # For identical coordinates, should be close to identity
        np.testing.assert_allclose(transform, np.eye(4), atol=1e-10)

    def test_regularize_group(self):
        """Test group regularization."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.coarse_grainer = self.coarse_grainer
        builder.parser = self.parser

        # Create group with multiple members
        group = Mock(spec=ChainGroup)
        group.representative = "A"
        group.members = ["A", "C"]  # A is reference, C needs regularization

        # Add chain C to mock data
        chain_c = Mock(spec=CoarseGrainedChain)
        chain_c.chain_id = "C"
        chain_c.com = np.array([20.0, 0.0, 0.0])
        self.coarse_grainer.get_coarse_grained_chains.return_value["C"] = chain_c

        # Mock compute_rigid_transform
        with patch.object(builder, '_compute_rigid_transform', return_value=np.eye(4)):
            builder._regularize_group(group)

        # Check that transform was computed and stored
        self.assertTrue(hasattr(chain_c, 'transform_from_reference'))

    def test_detect_steric_clashes(self):
        """Test steric clash detection."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.interface_templates = {}

        # Create mock interface templates on same molecule
        template1 = Mock(spec=InterfaceType)
        template1.this_mol_type = Mock(name="A")
        template1.local_coord = np.array([0.0, 0.0, 0.0])
        template1.required_free = []

        template2 = Mock(spec=InterfaceType)
        template2.this_mol_type = template1.this_mol_type  # Same molecule
        template2.local_coord = np.array([0.2, 0.0, 0.0])  # Close position
        template2.required_free = []

        builder.interface_templates["intf1"] = template1
        builder.interface_templates["intf2"] = template2

        builder._detect_steric_clashes()

        # Check that mutual exclusion was set up
        self.assertIn("intf2", template1.required_free)
        self.assertIn("intf1", template2.required_free)

    def test_get_molecule_templates(self):
        """Test getting molecule templates."""
        builder = TemplateBuilder.__new__(TemplateBuilder)

        template_a = Mock(spec=MoleculeType)
        template_b = Mock(spec=MoleculeType)

        builder.molecule_templates = {
            "A": template_a,
            "B": template_b
        }

        templates = builder.get_molecule_templates()

        self.assertEqual(len(templates), 2)
        self.assertIn("A", templates)
        self.assertIn("B", templates)
        self.assertEqual(templates["A"], template_a)
        self.assertEqual(templates["B"], template_b)

    def test_get_interface_templates(self):
        """Test getting interface templates."""
        builder = TemplateBuilder.__new__(TemplateBuilder)

        template_ab = Mock(spec=InterfaceType)
        template_ba = Mock(spec=InterfaceType)

        builder.interface_templates = {
            "A_B_1": template_ab,
            "B_A_1": template_ba
        }

        templates = builder.get_interface_templates()

        self.assertEqual(len(templates), 2)
        self.assertIn("A_B_1", templates)
        self.assertIn("B_A_1", templates)

    def test_get_template_name_for_group(self):
        """Test getting template name for group."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.group_to_template = {
            "A": "ProteinA",
            "B": "ProteinB"
        }

        name = builder.get_template_name_for_group("A")
        self.assertEqual(name, "ProteinA")

        name = builder.get_template_name_for_group("C")
        self.assertIsNone(name)

    def test_get_chain_name_mapping(self):
        """Test getting chain name mapping."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.chain_grouper = self.chain_grouper
        builder.group_to_template = {
            "A": "ProteinA",
            "B": "ProteinB"
        }

        mapping = builder.get_chain_name_mapping()

        expected = {
            "A": "ProteinA",
            "B": "ProteinB"
        }

        self.assertEqual(mapping, expected)

    def test_get_summary(self):
        """Test getting summary."""
        builder = TemplateBuilder.__new__(TemplateBuilder)
        builder.molecule_templates = {
            "A": Mock(spec=MoleculeType),
            "B": Mock(spec=MoleculeType)
        }

        # Mock interface templates
        template_ab = Mock(spec=InterfaceType)
        template_ab.this_mol_type_name = "A"
        template_ab.partner_mol_type_name = "B"

        template_ba = Mock(spec=InterfaceType)
        template_ba.this_mol_type_name = "B"
        template_ba.partner_mol_type_name = "A"

        builder.interface_templates = {
            "A_B_1": template_ab,
            "B_A_1": template_ba
        }

        builder.group_to_template = {"A": "A", "B": "B"}

        # Mock chain grouper for get_chain_name_mapping
        with patch.object(builder, 'get_chain_name_mapping', return_value={"A": "A", "B": "B"}):
            summary = builder.get_summary()

        self.assertEqual(summary["num_molecule_templates"], 2)
        self.assertEqual(summary["num_interface_templates"], 2)
        self.assertIn("A", summary["molecule_templates"])
        self.assertIn("B", summary["molecule_templates"])
        self.assertIn("A_B_1", summary["interface_templates"])
        self.assertIn("B_A_1", summary["interface_templates"])


class TestTemplateBuilderIntegration(unittest.TestCase):
    """Integration tests for TemplateBuilder."""

    def setUp(self):
        """Set up integration test fixtures."""
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.logger = Mock()

    def test_complete_template_building_workflow(self):
        """Test complete template building workflow."""
        # Create comprehensive mock setup
        parser = Mock(spec=PDBParser)
        coarse_grainer = Mock(spec=CoarseGrainer)
        chain_grouper = Mock(spec=ChainGrouper)
        hyperparams = Mock(spec=PDBModelHyperparameters)

        # Set up realistic mock data
        self._setup_integration_mocks(
            parser, coarse_grainer, chain_grouper, hyperparams)

        # Build templates
        builder = TemplateBuilder(
            parser=parser,
            coarse_grainer=coarse_grainer,
            chain_grouper=chain_grouper,
            hyperparams=hyperparams,
            workspace_manager=self.workspace_manager
        )

        # Verify templates were created
        molecule_templates = builder.get_molecule_templates()
        interface_templates = builder.get_interface_templates()

        self.assertGreater(len(molecule_templates), 0)
        self.assertGreater(len(interface_templates), 0)

        # Verify summary
        summary = builder.get_summary()
        self.assertGreater(summary["num_molecule_templates"], 0)
        self.assertGreater(summary["num_interface_templates"], 0)

    def _setup_integration_mocks(self, parser, coarse_grainer, chain_grouper, hyperparams):
        """Set up comprehensive mocks for integration testing."""
        # Mock hyperparameters
        hyperparams.homodimer_distance_threshold = 1.0
        hyperparams.homodimer_angle_threshold = 0.2
        hyperparams.signature_precision = 6
        hyperparams.steric_clash_mode = "off"

        # Mock coarse-grained chains
        chain_a = Mock(spec=CoarseGrainedChain)
        chain_a.chain_id = "A"
        chain_a.com = np.array([0.0, 0.0, 0.0])
        chain_a.radius = 15.0

        coarse_grainer.get_coarse_grained_chains.return_value = {"A": chain_a}

        # Mock interfaces (self-interaction)
        interface = Mock(spec=InterfaceString)
        interface.chain_i = "A"
        interface.chain_j = "A"
        interface.coord_i = np.array([5.0, 0.0, 0.0])
        interface.coord_j = np.array([-5.0, 0.0, 0.0])
        interface.residues_i = {1, 2}
        interface.residues_j = {3, 4}
        interface.energy = -3.0

        coarse_grainer.get_interfaces.return_value = [interface]

        # Mock chain group
        group = Mock(spec=ChainGroup)
        group.representative = "A"
        group.members = ["A"]
        group.grouping_method = "sequence_similarity"

        chain_grouper.get_groups.return_value = [group]
        chain_grouper.get_group_for_chain.return_value = group

        # Mock parser
        parser.convert_coords_to_nm.side_effect = lambda coords: coords / 10.0
        parser.get_chain_data.return_value = {
            'ca_coords': np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        }


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

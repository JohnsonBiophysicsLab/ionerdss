"""
Unit tests for ionerdss.model.pdb.system_builder

Tests the SystemBuilder class and its system assembly capabilities.
"""

import unittest
from unittest.mock import Mock, MagicMock, patch
from pathlib import Path
import numpy as np

from ionerdss.model.pdb.system_builder import SystemBuilder
from ionerdss.model.pdb.parser import PDBParser
from ionerdss.model.pdb.coarse_graining import CoarseGrainer, CoarseGrainedChain
from ionerdss.model.pdb.chain_grouping import ChainGrouper, ChainGroup
from ionerdss.model.pdb.template_builder import TemplateBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.file_manager import WorkspaceManager
from ionerdss.model.components.system import System
from ionerdss.model.components.instances import MoleculeInstance, InterfaceInstance
from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.components.units import Units


class MockRegistry:
    """Mock registry that supports len() and add()."""

    def __init__(self):
        self.items = []

    def add(self, item):
        self.items.append(item)

    def __len__(self):
        return len(self.items)

    def __iter__(self):
        return iter(self.items)


class TestSystemBuilder(unittest.TestCase):
    """Test cases for SystemBuilder class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock components
        self.parser = Mock(spec=PDBParser)
        self.coarse_grainer = Mock(spec=CoarseGrainer)
        self.chain_grouper = Mock(spec=ChainGrouper)
        self.template_builder = Mock(spec=TemplateBuilder)
        self.hyperparams = Mock(spec=PDBModelHyperparameters)
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.logger = Mock()

        # Set up workspace path and PDB ID
        self.workspace_path = "/test/workspace"
        self.pdb_id = "1ABC"

        # Configure hyperparameters with default values to avoid ring regularization
        self._setup_hyperparameters()

        # Set up mock data
        self._setup_mock_data()

    def _setup_hyperparameters(self):
        """Set up hyperparameters to avoid ring regularization by default."""
        # Don't set ring_regularization_mode attribute by default
        # This will make hasattr() return False and skip ring regularization
        pass

    def _setup_mock_data(self):
        """Set up mock data for testing."""
        # Mock coarse-grained chains
        chain_a = Mock(spec=CoarseGrainedChain)
        chain_a.com = np.array([10.0, 20.0, 30.0])  # Angstroms

        chain_b = Mock(spec=CoarseGrainedChain)
        chain_b.com = np.array([40.0, 50.0, 60.0])  # Angstroms

        self.coarse_grainer.get_coarse_grained_chains.return_value = {
            "A": chain_a,
            "B": chain_b
        }

        # Mock interfaces
        mock_interface = Mock()
        mock_interface.chain_i = "A"
        mock_interface.chain_j = "B"
        mock_interface.coord_i = np.array([15.0, 25.0, 35.0])
        mock_interface.coord_j = np.array([35.0, 45.0, 55.0])
        mock_interface.residues_i = {1, 2, 3}
        mock_interface.residues_j = {4, 5, 6}
        mock_interface.energy = -5.0
        # Add interface_type attribute that the code checks for
        mock_interface.interface_type = None  # Will trigger fallback lookup

        self.coarse_grainer.get_interfaces.return_value = [mock_interface]

        # Mock chain groups
        group_a = Mock(spec=ChainGroup)
        group_a.representative = "A"
        group_a.chains = ["A"]

        group_b = Mock(spec=ChainGroup)
        group_b.representative = "B"
        group_b.chains = ["B"]

        self.chain_grouper.get_group_for_chain.side_effect = lambda chain_id: {
            "A": group_a,
            "B": group_b
        }.get(chain_id)

        # Mock molecule templates with all required attributes
        mol_type_a = Mock(spec=MoleculeType)
        mol_type_a.name = "ProteinA"
        mol_type_a.this_mol_type_name = "ProteinA"  # Add required attribute

        mol_type_b = Mock(spec=MoleculeType)
        mol_type_b.name = "ProteinB"
        mol_type_b.this_mol_type_name = "ProteinB"  # Add required attribute

        self.template_builder.get_template_name_for_group.side_effect = lambda group: {
            "A": "ProteinA",
            "B": "ProteinB"
        }.get(group)

        self.template_builder.molecule_templates = {
            "ProteinA": mol_type_a,
            "ProteinB": mol_type_b
        }

        self.template_builder.get_molecule_templates.return_value = {
            "ProteinA": mol_type_a,
            "ProteinB": mol_type_b
        }

        # Mock interface templates with all required attributes
        interface_type = Mock(spec=InterfaceType)
        interface_type.get_name.return_value = "A_B_1"
        interface_type.interface_index = 1
        interface_type.partner_interface_type = None
        interface_type.this_mol_type_name = "ProteinA"  # Add required attribute
        interface_type.partner_mol_type_name = "ProteinB"  # Add required attribute

        self.template_builder.get_interface_type_for_interface.return_value = "A_B_1"
        self.template_builder.interface_templates = {
            "A_B_1": interface_type
        }
        self.template_builder.get_interface_templates.return_value = {
            "A_B_1": interface_type
        }

        self.template_builder.group_to_template = {
            "A": "ProteinA",
            "B": "ProteinB"
        }

        # Mock parser coordinate conversion
        self.parser.convert_coords_to_nm.side_effect = lambda coords: coords / 10.0

    def _create_mock_system(self):
        """Create a properly mocked system with registries that support len()."""
        mock_system = Mock(spec=System)
        mock_system.molecule_types = MockRegistry()
        mock_system.interface_types = MockRegistry()
        mock_system.molecule_instances = MockRegistry()
        mock_system.interface_instances = MockRegistry()
        mock_system._rebuild_cross_references = Mock()
        mock_system.get_summary.return_value = {
            "molecule_types": 2, "interface_types": 1}
        mock_system.validate_system.return_value = {
            "errors": [], "warnings": []}
        return mock_system

    def test_initialization(self):
        """Test SystemBuilder initialization."""
        with patch.object(SystemBuilder, '_build_system'):
            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                pdb_id=self.pdb_id,
                workspace_manager=self.workspace_manager
            )

        self.assertEqual(builder.parser, self.parser)
        self.assertEqual(builder.coarse_grainer, self.coarse_grainer)
        self.assertEqual(builder.chain_grouper, self.chain_grouper)
        self.assertEqual(builder.template_builder, self.template_builder)
        self.assertEqual(builder.hyperparams, self.hyperparams)
        self.assertEqual(builder.workspace_path, self.workspace_path)
        self.assertEqual(builder.pdb_id, self.pdb_id)
        self.assertEqual(builder.workspace_manager, self.workspace_manager)
        self.assertIsInstance(builder.units, Units)

    def test_initialization_with_custom_units(self):
        """Test SystemBuilder initialization with custom units."""
        custom_units = Mock(spec=Units)

        with patch.object(SystemBuilder, '_build_system'):
            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                units=custom_units,
                workspace_manager=self.workspace_manager
            )

        self.assertEqual(builder.units, custom_units)

    def test_create_molecule_instances(self):
        """Test _create_molecule_instances method."""
        builder = SystemBuilder.__new__(SystemBuilder)
        builder.parser = self.parser
        builder.coarse_grainer = self.coarse_grainer
        builder.chain_grouper = self.chain_grouper
        builder.template_builder = self.template_builder
        builder.workspace_manager = self.workspace_manager

        instances = builder._create_molecule_instances()

        # Should create instances for both chains
        self.assertEqual(len(instances), 2)

        # Check instance properties
        instance_names = [inst.name for inst in instances]
        self.assertIn("A_ProteinA", instance_names)
        self.assertIn("B_ProteinB", instance_names)

        # Check that coordinates were converted to nm
        for instance in instances:
            self.assertIsInstance(instance, MoleculeInstance)
            self.assertIsInstance(instance.com, np.ndarray)
            # COM should be in nm (original / 10)
            # Should be much smaller after conversion
            self.assertTrue(np.all(instance.com < 10.0))

    def test_create_molecule_instances_missing_group(self):
        """Test _create_molecule_instances with missing group."""
        builder = SystemBuilder.__new__(SystemBuilder)
        builder.parser = self.parser
        builder.coarse_grainer = self.coarse_grainer
        builder.chain_grouper = self.chain_grouper
        builder.template_builder = self.template_builder
        builder.workspace_manager = self.workspace_manager

        # Mock missing group for chain A
        self.chain_grouper.get_group_for_chain.side_effect = lambda chain_id: {
            "B": Mock(representative="B")
        }.get(chain_id)  # Chain A returns None

        instances = builder._create_molecule_instances()

        # Should only create instance for chain B
        self.assertEqual(len(instances), 1)
        self.assertEqual(instances[0].name, "B_ProteinB")

    def test_create_interface_instances(self):
        """Test _create_interface_instances method."""
        builder = SystemBuilder.__new__(SystemBuilder)
        builder.parser = self.parser
        builder.coarse_grainer = self.coarse_grainer
        builder.chain_grouper = self.chain_grouper
        builder.template_builder = self.template_builder
        builder.workspace_manager = self.workspace_manager

        instances = builder._create_interface_instances()

        # Should create bidirectional interface instances (2 per interface)
        self.assertEqual(len(instances), 2)

        # Check instance properties
        for instance in instances:
            self.assertIsInstance(instance, InterfaceInstance)
            self.assertIsInstance(instance.absolute_coord, np.ndarray)
            self.assertIsNotNone(instance.interface_type)
            self.assertIsNotNone(instance.this_mol_name)
            self.assertIsNotNone(instance.partner_mol_name)

        # Check bidirectional linking
        instance_i, instance_j = instances
        self.assertEqual(instance_i.partner_interface, instance_j)
        self.assertEqual(instance_j.partner_interface, instance_i)

    def test_create_interface_instances_missing_template(self):
        """Test _create_interface_instances with missing interface template."""
        builder = SystemBuilder.__new__(SystemBuilder)
        builder.parser = self.parser
        builder.coarse_grainer = self.coarse_grainer
        builder.chain_grouper = self.chain_grouper
        builder.template_builder = self.template_builder
        builder.workspace_manager = self.workspace_manager

        # Mock missing interface template
        self.template_builder.get_interface_type_for_interface.return_value = None

        instances = builder._create_interface_instances()

        # Should create no instances due to missing template
        self.assertEqual(len(instances), 0)

    def test_establish_cross_references(self):
        """Test _establish_cross_references method."""
        builder = SystemBuilder.__new__(SystemBuilder)
        builder.workspace_manager = self.workspace_manager

        # Create mock molecule instances
        mol_a = Mock(spec=MoleculeInstance)
        mol_a.name = "A_ProteinA"
        mol_a.interfaces_neighbors_map = {}

        mol_b = Mock(spec=MoleculeInstance)
        mol_b.name = "B_ProteinB"
        mol_b.interfaces_neighbors_map = {}

        builder.molecule_instances = [mol_a, mol_b]

        # Create mock interface instances
        intf_a = Mock(spec=InterfaceInstance)
        intf_a.this_mol_name = "A_ProteinA"
        intf_a.partner_mol_name = "B_ProteinB"
        intf_a.interface_index = 1
        intf_a.get_name.return_value = "A_ProteinA_B_ProteinB_1"

        intf_b = Mock(spec=InterfaceInstance)
        intf_b.this_mol_name = "B_ProteinB"
        intf_b.partner_mol_name = "A_ProteinA"
        intf_b.interface_index = 1
        intf_b.get_name.return_value = "B_ProteinB_A_ProteinA_1"

        builder.interface_instances = [intf_a, intf_b]

        # Establish cross-references
        builder._establish_cross_references()

        # Check that this_mol references were set
        self.assertEqual(intf_a.this_mol, mol_a)
        self.assertEqual(intf_b.this_mol, mol_b)

        # Check that partner interfaces were linked
        self.assertEqual(intf_a.partner_interface, intf_b)
        self.assertEqual(intf_b.partner_interface, intf_a)

        # Check that interfaces_neighbors_map was populated
        self.assertIn(intf_a, mol_a.interfaces_neighbors_map)
        self.assertEqual(mol_a.interfaces_neighbors_map[intf_a], mol_b)
        self.assertIn(intf_b, mol_b.interfaces_neighbors_map)
        self.assertEqual(mol_b.interfaces_neighbors_map[intf_b], mol_a)

    def test_create_system(self):
        """Test _create_system method."""
        builder = SystemBuilder.__new__(SystemBuilder)
        builder.workspace_path = self.workspace_path
        builder.pdb_id = self.pdb_id
        builder.units = Units()
        builder.template_builder = self.template_builder
        builder.workspace_manager = self.workspace_manager

        # Mock instances with required attributes
        mock_mol_instance = Mock(spec=MoleculeInstance)
        mock_mol_instance.name = "A_ProteinA"

        mock_interface_instance = Mock(spec=InterfaceInstance)
        mock_interface_instance.name = "A_B_1_instance"

        builder.molecule_instances = [mock_mol_instance]
        builder.interface_instances = [mock_interface_instance]

        # Mock the system with proper registries
        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            # Create system
            builder._create_system()

            # Check that system was created
            self.assertIsInstance(builder.system, Mock)
            mock_system._rebuild_cross_references.assert_called_once()

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    def test_build_system_with_ring_regularization(self, mock_ring_regularizer_class):
        """Test _build_system with ring regularization enabled."""
        # Mock hyperparameters with ring regularization - use actual string values
        self.hyperparams.ring_regularization_mode = "separate"
        self.hyperparams.ring_geometry = "cylinder"

        # Mock ring regularizer
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        # Mock the System class with proper registries
        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            # Create builder (this will call _build_system)
            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                workspace_manager=self.workspace_manager
            )

            # Check that ring regularizer was created and called
            mock_ring_regularizer_class.assert_called_once_with(
                system=builder.system,
                workspace_manager=self.workspace_manager,
                mode="separate",
                geometry="cylinder"
            )
            mock_regularizer.regularize.assert_called_once()

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    def test_build_system_without_ring_regularization(self, mock_ring_regularizer_class):
        """Test _build_system without ring regularization."""
        # Mock ring regularizer to avoid validation issues
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        # Mock the System class with proper registries
        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            # Create builder (should not raise error)
            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                workspace_manager=self.workspace_manager
            )

            # Should complete successfully without ring regularization
            self.assertIsInstance(builder.system, Mock)

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    def test_get_system(self, mock_ring_regularizer_class):
        """Test get_system method."""
        # Mock ring regularizer to avoid validation issues
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                workspace_manager=self.workspace_manager
            )

            system = builder.get_system()

            self.assertEqual(system, builder.system)

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    def test_validate_system(self, mock_ring_regularizer_class):
        """Test validate_system method."""
        # Mock ring regularizer to avoid validation issues
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                workspace_manager=self.workspace_manager
            )

            # Mock system validation
            mock_validation = {"errors": [], "warnings": ["test warning"]}
            builder.system.validate_system.return_value = mock_validation

            validation = builder.validate_system()

            self.assertEqual(validation, mock_validation)
            #builder.system.validate_system.assert_called_once()

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    def test_get_summary(self, mock_ring_regularizer_class):
        """Test get_summary method."""
        # Mock ring regularizer to avoid validation issues
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                workspace_manager=self.workspace_manager
            )

            # Mock hyperparams dict
            mock_hyperparams_dict = {"distance_cutoff": 0.6}
            self.hyperparams.to_dict.return_value = mock_hyperparams_dict

            summary = builder.get_summary()

            # Check that summary contains all expected components
            self.assertEqual(summary["molecule_types"], 2)
            self.assertEqual(summary["interface_types"], 1)
            self.assertIn("validation", summary)
            self.assertEqual(summary["hyperparameters"], mock_hyperparams_dict)

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    @patch('ionerdss.model.pdb.system_builder.PDBVisualizer')
    def test_generate_visualizations(self, mock_visualizer_class, mock_ring_regularizer_class):
        """Test generate_visualizations method."""
        # Mock ring regularizer to avoid validation issues
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                workspace_manager=self.workspace_manager
            )

            # Mock visualizer
            mock_visualizer = Mock()
            mock_viz_outputs = {
                "structure": Path("/test/structure.png"),
                "interfaces": Path("/test/interfaces.png")
            }
            mock_visualizer.visualize_all.return_value = mock_viz_outputs
            mock_visualizer_class.return_value = mock_visualizer

            viz_outputs = builder.generate_visualizations()

            # Check that visualizer was created and called
            mock_visualizer_class.assert_called_once_with(
                self.workspace_manager)
            mock_visualizer.visualize_all.assert_called_once_with(
                self.parser, self.coarse_grainer, self.chain_grouper, self.template_builder
            )

            self.assertEqual(viz_outputs, mock_viz_outputs)

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    def test_generate_visualizations_no_workspace_manager(self, mock_ring_regularizer_class):
        """Test generate_visualizations without workspace manager."""
        # Mock ring regularizer to avoid validation issues
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                workspace_manager=None  # No workspace manager
            )

            viz_outputs = builder.generate_visualizations()

            # Should return empty dict
            self.assertEqual(viz_outputs, {})

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    @patch('ionerdss.model.pdb.system_builder.NERDSSExporter')
    def test_export_nerdss_files(self, mock_exporter_class, mock_ring_regularizer_class):
        """Test export_nerdss_files method."""
        # Mock ring regularizer to avoid validation issues
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                workspace_manager=self.workspace_manager
            )

            # Mock exporter
            mock_exporter = Mock()
            mock_export_outputs = {
                "ProteinA_mol": Path("/test/ProteinA.mol"),
                "parms": Path("/test/parms.inp")
            }
            mock_exporter.export_all.return_value = mock_export_outputs
            mock_exporter_class.return_value = mock_exporter

            # Test export with custom parameters
            molecule_counts = {"ProteinA": 10, "ProteinB": 5}
            box_nm = (200.0, 200.0, 200.0)
            parms_overrides = {"timestep": 0.1}

            export_outputs = builder.export_nerdss_files(
                molecule_counts=molecule_counts,
                box_nm=box_nm,
                parms_overrides=parms_overrides
            )

            # Check that exporter was created and called correctly
            mock_exporter_class.assert_called_once_with(
                builder.system, self.workspace_manager)
            mock_exporter.export_all.assert_called_once_with(
                molecule_counts=molecule_counts,
                box_nm=box_nm,
                parms_overrides=parms_overrides
            )

            self.assertEqual(export_outputs, mock_export_outputs)

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    def test_export_nerdss_files_default_parameters(self, mock_ring_regularizer_class):
        """Test export_nerdss_files with default parameters."""
        # Mock ring regularizer to avoid validation issues
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            builder = SystemBuilder(
                parser=self.parser,
                coarse_grainer=self.coarse_grainer,
                chain_grouper=self.chain_grouper,
                template_builder=self.template_builder,
                hyperparams=self.hyperparams,
                workspace_path=self.workspace_path,
                workspace_manager=self.workspace_manager
            )

            with patch('ionerdss.model.pdb.system_builder.NERDSSExporter') as mock_exporter_class:
                mock_exporter = Mock()
                mock_exporter_class.return_value = mock_exporter
                mock_exporter.export_all.return_value = {}

                builder.export_nerdss_files()

                # Should call with default parameters
                mock_exporter.export_all.assert_called_once_with(
                    molecule_counts=None,
                    box_nm=(100.0, 100.0, 100.0),
                    parms_overrides=None
                )


class TestSystemBuilderIntegration(unittest.TestCase):
    """Integration tests for SystemBuilder."""

    def setUp(self):
        """Set up integration test fixtures."""
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.logger = Mock()

    def _create_mock_system(self):
        """Create a properly mocked system with registries that support len()."""
        mock_system = Mock(spec=System)
        mock_system.molecule_types = MockRegistry()
        mock_system.interface_types = MockRegistry()
        mock_system.molecule_instances = MockRegistry()
        mock_system.interface_instances = MockRegistry()
        mock_system._rebuild_cross_references = Mock()
        mock_system.get_summary.return_value = {
            "molecule_types": 1, "interface_types": 1}
        mock_system.validate_system.return_value = {
            "errors": [], "warnings": []}
        return mock_system

    @patch('ionerdss.model.pdb.system_builder.RingRegularizer')
    def test_complete_system_building_workflow(self, mock_ring_regularizer_class):
        """Test complete system building workflow."""
        # Mock ring regularizer to avoid validation issues
        mock_regularizer = Mock()
        mock_ring_regularizer_class.return_value = mock_regularizer

        # Create comprehensive mock setup
        parser = Mock(spec=PDBParser)
        coarse_grainer = Mock(spec=CoarseGrainer)
        chain_grouper = Mock(spec=ChainGrouper)
        template_builder = Mock(spec=TemplateBuilder)
        hyperparams = Mock(spec=PDBModelHyperparameters)

        # Set up realistic mock data
        self._setup_integration_mocks(
            parser, coarse_grainer, chain_grouper, template_builder
        )

        # Mock the System class with proper registries
        with patch('ionerdss.model.pdb.system_builder.System') as mock_system_class:
            mock_system = self._create_mock_system()
            mock_system_class.return_value = mock_system

            # Build system
            builder = SystemBuilder(
                parser=parser,
                coarse_grainer=coarse_grainer,
                chain_grouper=chain_grouper,
                template_builder=template_builder,
                hyperparams=hyperparams,
                workspace_path="/test/workspace",
                pdb_id="1ABC",
                workspace_manager=self.workspace_manager
            )

            # Verify system was built
            system = builder.get_system()
            self.assertIsInstance(system, Mock)

            # Verify components were created
            self.assertGreater(len(builder.molecule_instances), 0)
            self.assertGreater(len(builder.interface_instances), 0)

            # Verify summary contains expected information
            hyperparams.to_dict.return_value = {"distance_cutoff": 0.6}
            summary = builder.get_summary()
            self.assertIn("validation", summary)
            self.assertIn("hyperparameters", summary)

    def _setup_integration_mocks(self, parser, coarse_grainer, chain_grouper, template_builder):
        """Set up comprehensive mocks for integration testing."""
        # Mock coarse-grained chains
        chain_data = Mock(spec=CoarseGrainedChain)
        chain_data.com = np.array([10.0, 20.0, 30.0])

        coarse_grainer.get_coarse_grained_chains.return_value = {
            "A": chain_data}

        # Mock interfaces
        interface = Mock()
        interface.chain_i = "A"
        interface.chain_j = "A"  # Self-interaction
        interface.coord_i = np.array([15.0, 25.0, 35.0])
        interface.coord_j = np.array([25.0, 35.0, 45.0])
        interface.residues_i = {1, 2}
        interface.residues_j = {3, 4}
        interface.energy = -3.0
        interface.interface_type = None  # Will trigger fallback lookup

        coarse_grainer.get_interfaces.return_value = [interface]

        # Mock chain group
        group = Mock(spec=ChainGroup)
        group.representative = "A"
        chain_grouper.get_group_for_chain.return_value = group

        # Mock templates
        mol_type = Mock(spec=MoleculeType)
        mol_type.name = "ProteinA"
        mol_type.this_mol_type_name = "ProteinA"

        interface_type = Mock(spec=InterfaceType)
        interface_type.get_name.return_value = "A_A_1"
        interface_type.interface_index = 1
        interface_type.partner_interface_type = None
        interface_type.this_mol_type_name = "ProteinA"
        interface_type.partner_mol_type_name = "ProteinA"

        template_builder.get_template_name_for_group.return_value = "ProteinA"
        template_builder.molecule_templates = {"ProteinA": mol_type}
        template_builder.get_molecule_templates.return_value = {
            "ProteinA": mol_type}
        template_builder.get_interface_type_for_interface.return_value = "A_A_1"
        template_builder.interface_templates = {"A_A_1": interface_type}
        template_builder.get_interface_templates.return_value = {
            "A_A_1": interface_type}
        template_builder.group_to_template = {"A": "ProteinA"}

        # Mock parser
        parser.convert_coords_to_nm.side_effect = lambda coords: coords / 10.0


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

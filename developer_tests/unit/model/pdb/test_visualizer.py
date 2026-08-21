"""
Unit tests for ionerdss.model.pdb.visualizer

Tests the PDBVisualizer class and its visualization generation capabilities.

"""

import unittest
from unittest.mock import Mock, patch, mock_open
from pathlib import Path
import tempfile
import shutil
import numpy as np

from ionerdss.model.pdb.visualizer import PDBVisualizer
from ionerdss.model.pdb.file_manager import WorkspaceManager
from ionerdss.model.pdb.coarse_graining import CoarseGrainer, CoarseGrainedChain, InterfaceString
from ionerdss.model.pdb.chain_grouping import ChainGrouper, ChainGroup
from ionerdss.model.pdb.template_builder import TemplateBuilder
from ionerdss.model.components.types import MoleculeType, InterfaceType


class TestPDBVisualizer(unittest.TestCase):
    """Test cases for PDBVisualizer class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create temporary directory for testing
        self.temp_dir = Path(tempfile.mkdtemp())

        # Create mock workspace manager
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.workspace_path = self.temp_dir
        self.workspace_manager.pdb_id = "1ABC"
        self.workspace_manager.logger = Mock()

        # Set up mock data
        self._setup_mock_data()

    def tearDown(self):
        """Clean up test fixtures."""
        # Remove temporary directory
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    def _setup_mock_data(self):
        """Set up mock data for testing."""
        # Mock coarse-grained chains
        self.chain_a = Mock(spec=CoarseGrainedChain)
        self.chain_a.chain_id = "A"
        self.chain_a.com = np.array([0.0, 0.0, 0.0])
        self.chain_a.radius = 15.0
        self.chain_a.sequence = "MKLAVQN"

        self.chain_b = Mock(spec=CoarseGrainedChain)
        self.chain_b.chain_id = "B"
        self.chain_b.com = np.array([20.0, 0.0, 0.0])
        self.chain_b.radius = 12.0
        self.chain_b.sequence = "AKLVQNC"

        # Mock coarse grainer
        self.coarse_grainer = Mock(spec=CoarseGrainer)
        self.coarse_grainer.get_coarse_grained_chains.return_value = {
            "A": self.chain_a,
            "B": self.chain_b
        }

        # Mock interfaces
        self.interface_ab = Mock(spec=InterfaceString)
        self.interface_ab.chain_i = "A"
        self.interface_ab.chain_j = "B"
        self.interface_ab.coord_i = np.array([10.0, 0.0, 0.0])
        self.interface_ab.coord_j = np.array([10.0, 0.0, 0.0])
        self.interface_ab.residues_i = {1, 2, 3}
        self.interface_ab.residues_j = {4, 5, 6}
        self.interface_ab.energy = -5.0

        self.coarse_grainer.get_interfaces.return_value = [self.interface_ab]

        # Mock chain groups
        self.group_a = Mock(spec=ChainGroup)
        self.group_a.representative = "A"
        self.group_a.members = ["A"]
        self.group_a.grouping_method = "sequence_similarity"

        self.group_b = Mock(spec=ChainGroup)
        self.group_b.representative = "B"
        self.group_b.members = ["B"]
        self.group_b.grouping_method = "sequence_similarity"

        # Mock chain grouper
        self.chain_grouper = Mock(spec=ChainGrouper)
        self.chain_grouper.get_groups.return_value = [
            self.group_a, self.group_b]

        # Mock molecule templates
        self.mol_template_a = Mock(spec=MoleculeType)
        self.mol_template_a.name = "ProteinA"
        self.mol_template_a.radius_nm = 1.5
        self.mol_template_a.D_t_nm2_us = 1.0e-3
        self.mol_template_a.D_r_rad2_us = 1.0e-4

        self.mol_template_b = Mock(spec=MoleculeType)
        self.mol_template_b.name = "ProteinB"
        self.mol_template_b.radius_nm = 1.2
        self.mol_template_b.D_t_nm2_us = 1.2e-3
        self.mol_template_b.D_r_rad2_us = 1.2e-4

        # Mock interface templates
        self.intf_template = Mock(spec=InterfaceType)
        self.intf_template.energy = -5.0

        # Mock template builder
        self.template_builder = Mock(spec=TemplateBuilder)
        self.template_builder.get_molecule_templates.return_value = {
            "ProteinA": self.mol_template_a,
            "ProteinB": self.mol_template_b
        }
        self.template_builder.get_interface_templates.return_value = {
            "A_B_1": self.intf_template
        }
        self.template_builder.get_template_name_for_group.side_effect = lambda x: {
            "A": "ProteinA",
            "B": "ProteinB"
        }.get(x)

        # Mock parser
        self.parser = Mock()
        self.parser.filepath = Path("/test/structure.pdb")

    def test_initialization(self):
        """Test PDBVisualizer initialization."""
        visualizer = PDBVisualizer(self.workspace_manager)

        self.assertEqual(visualizer.workspace_manager, self.workspace_manager)
        self.assertEqual(visualizer.logger, self.workspace_manager.logger)

        # Check that visualization directory was created
        viz_dir = self.temp_dir / 'visualizations'
        self.assertTrue(viz_dir.exists())
        self.assertEqual(visualizer.viz_dir, viz_dir)

    @patch('ionerdss.model.pdb.visualizer.plt')
    def test_plot_basic_coarse_grained_structure(self, mock_plt):
        """Test basic coarse-grained structure plotting."""
        # Mock matplotlib components
        mock_fig = Mock()
        mock_ax = Mock()
        mock_plt.figure.return_value = mock_fig
        mock_fig.add_subplot.return_value = mock_ax

        visualizer = PDBVisualizer(self.workspace_manager)

        output_path = visualizer.plot_basic_coarse_grained_structure(
            self.coarse_grainer)

        # Check that plot was created
        expected_path = visualizer.viz_dir / 'basic_coarse_grained_structure.png'
        self.assertEqual(output_path, expected_path)

        # Check that matplotlib methods were called
        mock_plt.figure.assert_called_once()
        mock_fig.add_subplot.assert_called_once_with(111, projection='3d')
        mock_ax.scatter.assert_called()
        mock_ax.plot.assert_called()
        mock_plt.savefig.assert_called_once()
        mock_plt.close.assert_called_once()

    @patch('ionerdss.model.pdb.visualizer.plt')
    def test_plot_basic_coarse_grained_structure_empty_data(self, mock_plt):
        """Test basic coarse-grained structure plotting with empty data."""
        # Mock matplotlib components
        mock_fig = Mock()
        mock_ax = Mock()
        mock_plt.figure.return_value = mock_fig
        mock_fig.add_subplot.return_value = mock_ax

        # Mock empty coarse grainer
        empty_coarse_grainer = Mock(spec=CoarseGrainer)
        empty_coarse_grainer.get_coarse_grained_chains.return_value = {}
        empty_coarse_grainer.get_interfaces.return_value = []

        visualizer = PDBVisualizer(self.workspace_manager)

        output_path = visualizer.plot_basic_coarse_grained_structure(
            empty_coarse_grainer)

        # Should still create plot even with empty data
        expected_path = visualizer.viz_dir / 'basic_coarse_grained_structure.png'
        self.assertEqual(output_path, expected_path)

        mock_plt.figure.assert_called_once()
        mock_plt.savefig.assert_called_once()

    @patch('ionerdss.model.pdb.visualizer.plt')
    def test_plot_interface_connections(self, mock_plt):
        """Test interface connections plotting."""
        # Mock matplotlib components
        mock_fig = Mock()
        mock_ax = Mock()
        mock_plt.figure.return_value = mock_fig
        mock_fig.add_subplot.return_value = mock_ax

        visualizer = PDBVisualizer(self.workspace_manager)

        output_path = visualizer.plot_interface_connections(
            self.coarse_grainer)

        # Check that plot was created
        expected_path = visualizer.viz_dir / 'interface_connections.png'
        self.assertEqual(output_path, expected_path)

        # Check that matplotlib methods were called
        mock_plt.figure.assert_called_once()
        mock_ax.scatter.assert_called()
        mock_ax.plot.assert_called()
        mock_plt.savefig.assert_called_once()

    @patch('ionerdss.model.pdb.visualizer.plt')
    def test_plot_chain_groups(self, mock_plt):
        """Test chain groups plotting."""
        # Mock matplotlib components and colormaps
        mock_fig = Mock()
        mock_ax = Mock()
        mock_plt.figure.return_value = mock_fig
        mock_fig.add_subplot.return_value = mock_ax
        mock_plt.cm.Set1.return_value = np.array([[1, 0, 0], [0, 1, 0]])

        visualizer = PDBVisualizer(self.workspace_manager)

        output_path = visualizer.plot_chain_groups(
            self.coarse_grainer, self.chain_grouper)

        # Check that plot was created
        expected_path = visualizer.viz_dir / 'chain_groups.png'
        self.assertEqual(output_path, expected_path)

        # Check that matplotlib methods were called
        mock_plt.figure.assert_called_once()
        mock_ax.scatter.assert_called()
        mock_plt.savefig.assert_called_once()

    @patch('ionerdss.model.pdb.visualizer.plt')
    def test_plot_template_overview(self, mock_plt):
        """Test template overview plotting."""
        # Mock matplotlib components
        mock_fig = Mock()
        mock_axes = [[Mock(), Mock()], [Mock(), Mock()]]
        mock_plt.subplots.return_value = (mock_fig, mock_axes)

        visualizer = PDBVisualizer(self.workspace_manager)

        output_path = visualizer.plot_template_overview(self.template_builder)

        # Check that plot was created
        expected_path = visualizer.viz_dir / 'template_overview.png'
        self.assertEqual(output_path, expected_path)

        # Check that matplotlib methods were called
        mock_plt.subplots.assert_called_once_with(2, 2, figsize=(14, 10))
        mock_plt.tight_layout.assert_called_once()
        mock_plt.savefig.assert_called_once()

    @patch('ionerdss.model.pdb.visualizer.plt')
    def test_plot_template_overview_empty_interfaces(self, mock_plt):
        """Test template overview plotting with empty interface data."""
        # Mock matplotlib components
        mock_fig = Mock()
        mock_axes = [[Mock(), Mock()], [Mock(), Mock()]]
        mock_plt.subplots.return_value = (mock_fig, mock_axes)

        # Mock template builder with empty interfaces
        empty_template_builder = Mock(spec=TemplateBuilder)
        empty_template_builder.get_molecule_templates.return_value = {
            "ProteinA": self.mol_template_a
        }
        empty_template_builder.get_interface_templates.return_value = {}

        visualizer = PDBVisualizer(self.workspace_manager)

        output_path = visualizer.plot_template_overview(empty_template_builder)

        # Should still create plot
        expected_path = visualizer.viz_dir / 'template_overview.png'
        self.assertEqual(output_path, expected_path)

    def test_save_coarse_grained_structure(self):
        """Test saving coarse-grained structure to CIF."""
        visualizer = PDBVisualizer(self.workspace_manager)

        with patch('builtins.open', mock_open()) as mock_file:
            output_path = visualizer.save_coarse_grained_structure(
                self.coarse_grainer)

        # Check that file was created
        expected_path = visualizer.viz_dir / '1ABC_coarse_grained.cif'
        self.assertEqual(output_path, expected_path)

        # Check that file was opened for writing
        mock_file.assert_called_once_with(expected_path, 'w', encoding='utf-8')

        # Check that content was written
        handle = mock_file()
        handle.write.assert_called()

        # Verify some expected content
        written_content = ''.join(call.args[0]
                                  for call in handle.write.call_args_list)
        self.assertIn('data_coarse_grained', written_content)
        self.assertIn('_atom_site.group_PDB', written_content)
        self.assertIn('COM', written_content)
        self.assertIn('INT', written_content)

    def test_generate_pymol_script(self):
        """Test PyMOL script generation."""
        visualizer = PDBVisualizer(self.workspace_manager)

        with patch('builtins.open', mock_open()) as mock_file:
            output_path = visualizer.generate_pymol_script(
                self.parser, self.coarse_grainer, self.template_builder
            )

        # Check that file was created
        expected_path = visualizer.viz_dir / '1ABC_visualization.pml'
        self.assertEqual(output_path, expected_path)

        # Check that file was opened for writing
        mock_file.assert_called_once_with(expected_path, 'w', encoding='utf-8')

        # Check that content was written
        handle = mock_file()
        handle.write.assert_called()

        # Verify some expected content
        written_content = ''.join(call.args[0]
                                  for call in handle.write.call_args_list)
        self.assertIn('load', written_content)
        self.assertIn('pseudoatom', written_content)
        self.assertIn('distance', written_content)
        # self.assertIn('png', written_content)  # Deprecated: png generation muted by default

    def test_set_equal_aspect_3d(self):
        """Test setting equal aspect ratio for 3D plots."""
        visualizer = PDBVisualizer(self.workspace_manager)

        # Mock 3D axis
        mock_ax = Mock()

        # Test with valid coordinates
        coords = np.array([[0, 0, 0], [10, 5, 3], [-5, 8, -2]])

        visualizer._set_equal_aspect_3d(mock_ax, coords)

        # Check that limits were set
        mock_ax.set_xlim.assert_called_once()
        mock_ax.set_ylim.assert_called_once()
        mock_ax.set_zlim.assert_called_once()

    def test_set_equal_aspect_3d_empty_coords(self):
        """Test setting equal aspect ratio with empty coordinates."""
        visualizer = PDBVisualizer(self.workspace_manager)

        # Mock 3D axis
        mock_ax = Mock()

        # Test with empty coordinates
        coords = np.array([]).reshape(0, 3)

        # Should not raise error and not call set methods
        visualizer._set_equal_aspect_3d(mock_ax, coords)

        mock_ax.set_xlim.assert_not_called()
        mock_ax.set_ylim.assert_not_called()
        mock_ax.set_zlim.assert_not_called()

    def test_generate_summary_report(self):
        """Test generating summary report."""
        visualizer = PDBVisualizer(self.workspace_manager)

        # Create some test files in viz directory
        test_file1 = visualizer.viz_dir / 'test1.png'
        test_file2 = visualizer.viz_dir / 'test2.pml'
        test_file1.touch()
        test_file2.touch()

        with patch('builtins.open', mock_open()) as mock_file:
            output_path = visualizer.generate_summary_report(
                self.coarse_grainer, self.chain_grouper, self.template_builder
            )

        # Check that file was created
        expected_path = visualizer.viz_dir / 'visualization_summary.txt'
        self.assertEqual(output_path, expected_path)

        # Check that file was opened for writing
        mock_file.assert_called_once_with(expected_path, 'w', encoding='utf-8')

        # Check that content was written
        handle = mock_file()
        handle.write.assert_called()

        # Verify some expected content
        written_content = ''.join(call.args[0]
                                  for call in handle.write.call_args_list)
        self.assertIn('Visualization Summary Report', written_content)
        self.assertIn('1ABC', written_content)
        self.assertIn('Chains: 2', written_content)
        self.assertIn('Interfaces: 1', written_content)
        self.assertIn('test1.png', written_content)
        self.assertIn('test2.pml', written_content)

    @patch('ionerdss.model.pdb.visualizer.plt')
    def test_visualize_all(self, mock_plt):
        """Test generating all visualizations."""
        # Mock matplotlib components
        mock_fig = Mock()
        mock_ax = Mock()
        mock_plt.figure.return_value = mock_fig
        mock_fig.add_subplot.return_value = mock_ax
        mock_plt.subplots.return_value = (
            mock_fig, [[mock_ax, mock_ax], [mock_ax, mock_ax]])
        mock_plt.cm.Set1.return_value = np.array([[1, 0, 0], [0, 1, 0]])

        visualizer = PDBVisualizer(self.workspace_manager)

        with patch('builtins.open', mock_open()):
            outputs = visualizer.visualize_all(
                self.parser, self.coarse_grainer, self.chain_grouper, self.template_builder
            )

        # Check that all expected outputs were generated
        expected_keys = [
            'basic_cg', 'interfaces', 'groups', 'templates',
            'pymol', 'cg_structure'
        ]

        for key in expected_keys:
            self.assertIn(key, outputs)
            self.assertIsInstance(outputs[key], Path)

    @patch('ionerdss.model.pdb.visualizer.plt')
    def test_visualize_all_with_error(self, mock_plt):
        """Test visualize_all with error handling."""
        # Mock matplotlib to raise an error
        mock_plt.figure.side_effect = Exception("Test error")

        visualizer = PDBVisualizer(self.workspace_manager)

        # Should raise the exception
        with self.assertRaises(Exception) as context:
            visualizer.visualize_all(
                self.parser, self.coarse_grainer, self.chain_grouper, self.template_builder
            )

        self.assertEqual(str(context.exception), "Test error")

        # Check that error was logged
        visualizer.logger.error.assert_called()

    def test_plot_with_custom_figsize(self):
        """Test plotting with custom figure size."""
        with patch('ionerdss.model.pdb.visualizer.plt') as mock_plt:
            mock_fig = Mock()
            mock_ax = Mock()
            mock_plt.figure.return_value = mock_fig
            mock_fig.add_subplot.return_value = mock_ax

            visualizer = PDBVisualizer(self.workspace_manager)

            # Test with custom figure size
            custom_figsize = (16, 12)
            visualizer.plot_basic_coarse_grained_structure(
                self.coarse_grainer, figsize=custom_figsize
            )

            # Check that custom figsize was used
            mock_plt.figure.assert_called_with(figsize=custom_figsize)


class TestPDBVisualizerIntegration(unittest.TestCase):
    """Integration tests for PDBVisualizer."""

    def setUp(self):
        """Set up integration test fixtures."""
        self.temp_dir = Path(tempfile.mkdtemp())

        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.workspace_path = self.temp_dir
        self.workspace_manager.pdb_id = "TEST"
        self.workspace_manager.logger = Mock()

    def tearDown(self):
        """Clean up test fixtures."""
        if self.temp_dir.exists():
            shutil.rmtree(self.temp_dir)

    @patch('ionerdss.model.pdb.visualizer.plt')
    def test_complete_visualization_workflow(self, mock_plt):
        """Test complete visualization workflow."""
        # Mock matplotlib components
        mock_fig = Mock()
        mock_ax = Mock()
        mock_plt.figure.return_value = mock_fig
        mock_fig.add_subplot.return_value = mock_ax
        mock_plt.subplots.return_value = (
            mock_fig, [[mock_ax, mock_ax], [mock_ax, mock_ax]])
        mock_plt.cm.Set1.return_value = np.array([[1, 0, 0]])

        # Create comprehensive mock data
        parser, coarse_grainer, chain_grouper, template_builder = self._create_integration_mocks()

        # Create visualizer and generate all outputs
        visualizer = PDBVisualizer(self.workspace_manager)

        with patch('builtins.open', mock_open()):
            outputs = visualizer.visualize_all(
                parser, coarse_grainer, chain_grouper, template_builder
            )

        # Verify all outputs were generated
        self.assertGreater(len(outputs), 0)

        # Verify visualization directory was created
        viz_dir = self.temp_dir / 'visualizations'
        self.assertTrue(viz_dir.exists())

    def _create_integration_mocks(self):
        """Create comprehensive mocks for integration testing."""
        # Mock parser
        parser = Mock()
        parser.filepath = Path("/test/structure.pdb")

        # Mock coarse-grained chain
        chain = Mock(spec=CoarseGrainedChain)
        chain.chain_id = "A"
        chain.com = np.array([0.0, 0.0, 0.0])
        chain.radius = 15.0
        chain.sequence = "MKLA"

        # Mock coarse grainer
        coarse_grainer = Mock(spec=CoarseGrainer)
        coarse_grainer.get_coarse_grained_chains.return_value = {"A": chain}
        coarse_grainer.get_interfaces.return_value = []

        # Mock chain group
        group = Mock(spec=ChainGroup)
        group.representative = "A"
        group.members = ["A"]
        group.grouping_method = "single_chain"

        # Mock chain grouper
        chain_grouper = Mock(spec=ChainGrouper)
        chain_grouper.get_groups.return_value = [group]

        # Mock molecule template
        mol_template = Mock(spec=MoleculeType)
        mol_template.name = "ProteinA"
        mol_template.radius_nm = 1.5
        mol_template.D_t_nm2_us = 1.0e-3
        mol_template.D_r_rad2_us = 1.0e-4

        # Mock template builder
        template_builder = Mock(spec=TemplateBuilder)
        template_builder.get_molecule_templates.return_value = {
            "ProteinA": mol_template}
        template_builder.get_interface_templates.return_value = {}
        template_builder.get_template_name_for_group.return_value = "ProteinA"

        return parser, coarse_grainer, chain_grouper, template_builder


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

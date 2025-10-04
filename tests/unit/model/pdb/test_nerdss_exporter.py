"""
Unit tests for ionerdss.model.pdb.nerdss_exporter

Tests the NERDSSExporter class and its file generation capabilities.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch
import numpy as np

from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.components.instances import MoleculeInstance, InterfaceInstance
from ionerdss.model.pdb.file_manager import WorkspaceManager


class TestNERDSSExporter(unittest.TestCase):
    """Test cases for NERDSSExporter class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create temporary directory
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)

        # Create mock system
        self.system = self._create_mock_system()

        # Create mock workspace manager
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.workspace_path = self.temp_path
        self.workspace_manager.logger = Mock()

    def tearDown(self):
        """Clean up test fixtures."""
        if self.temp_path.exists():
            shutil.rmtree(self.temp_path)

    def _create_mock_system(self):
        """Create a mock system with molecule types and instances."""
        system = Mock(spec=System)

        # Create molecule type A
        mol_type_a = Mock(spec=MoleculeType)
        mol_type_a.name = "A"
        mol_type_a.D_t_nm2_us = 1.0
        mol_type_a.D_r_rad2_us = 0.1

        # Create interface types
        interface_type_1 = Mock(spec=InterfaceType)
        interface_type_1.get_name.return_value = "A_A_1"
        interface_type_1.this_mol_type_name = "A"
        interface_type_1.partner_mol_type_name = "A"
        interface_type_1.interface_index = 1
        interface_type_1.local_coord = np.array([1.0, 0.0, 0.0])
        interface_type_1.absolute_coord = np.array([1.0, 0.0, 0.0])

        interface_type_2 = Mock(spec=InterfaceType)
        interface_type_2.get_name.return_value = "A_A_2"
        interface_type_2.this_mol_type_name = "A"
        interface_type_2.partner_mol_type_name = "A"
        interface_type_2.interface_index = 2
        interface_type_2.local_coord = np.array([0.0, 1.0, 0.0])
        interface_type_2.absolute_coord = np.array([0.0, 1.0, 0.0])

        # Create molecule instances
        mol_instance_1 = Mock(spec=MoleculeInstance)
        mol_instance_1.molecule_type = mol_type_a
        mol_instance_1.com = np.array([0.0, 0.0, 0.0])

        mol_instance_2 = Mock(spec=MoleculeInstance)
        mol_instance_2.molecule_type = mol_type_a
        mol_instance_2.com = np.array([5.0, 0.0, 0.0])

        # Create interface instances
        interface_instance_1 = Mock(spec=InterfaceInstance)
        interface_instance_1.interface_type = interface_type_1
        interface_instance_1.absolute_coord = np.array([1.0, 0.0, 0.0])

        interface_instance_2 = Mock(spec=InterfaceInstance)
        interface_instance_2.interface_type = interface_type_2
        interface_instance_2.absolute_coord = np.array([0.0, 1.0, 0.0])

        # Set up partner interfaces
        interface_instance_1.partner_interface = interface_instance_2
        interface_instance_2.partner_interface = interface_instance_1

        # Set up interfaces_neighbors_map
        mol_instance_1.interfaces_neighbors_map = {
            interface_instance_1: mol_instance_2,
            interface_instance_2: mol_instance_2
        }

        mol_instance_2.interfaces_neighbors_map = {
            interface_instance_1: mol_instance_1,
            interface_instance_2: mol_instance_1
        }

        # Set up system
        system.molecule_types = [mol_type_a]
        system.interface_types = [interface_type_1, interface_type_2]
        system.molecule_instances = [mol_instance_1, mol_instance_2]

        return system

    def test_initialization_with_workspace_manager(self):
        """Test NERDSSExporter initialization with workspace manager."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Check basic attributes
        self.assertEqual(exporter.system, self.system)
        self.assertEqual(exporter.workspace_manager, self.workspace_manager)

        # Check output directory creation
        expected_output_dir = self.temp_path / 'nerdss_files'
        self.assertEqual(exporter.output_dir, expected_output_dir)
        self.assertTrue(expected_output_dir.exists())

        # Check empty mappings
        self.assertEqual(len(exporter.interface_to_site_map), 0)
        self.assertEqual(len(exporter.reaction_metadata), 0)
        self.assertEqual(len(exporter.homotypic_interface_map), 0)
        self.assertEqual(len(exporter.reaction_params_cache), 0)

    def test_initialization_without_workspace_manager(self):
        """Test NERDSSExporter initialization without workspace manager."""
        exporter = NERDSSExporter(self.system)

        self.assertEqual(exporter.system, self.system)
        self.assertIsNone(exporter.workspace_manager)

        # Should create default output directory
        expected_output_dir = Path("nerdss_files")
        self.assertEqual(exporter.output_dir, expected_output_dir)

    def test_get_base_site_label(self):
        """Test _get_base_site_label method."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Test with standard interface type name
        label = exporter._get_base_site_label("A", "A_B_1")
        self.assertEqual(label, "a1")

        # Test with different molecule name
        label = exporter._get_base_site_label("Protein", "A_B_2")
        self.assertEqual(label, "p2")

        # Test with interface type without index
        label = exporter._get_base_site_label("A", "simple_name")
        self.assertEqual(label, "a1")

        # Test with empty molecule name
        label = exporter._get_base_site_label("", "A_B_1")
        self.assertEqual(label, "x1")

    def test_get_unique_site_label_with_base(self):
        """Test _get_unique_site_label_with_base method."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Test with unused label
        used_labels = {"a1", "b1"}
        label = exporter._get_unique_site_label_with_base("a2", used_labels)
        self.assertEqual(label, "a2")

        # Test with used label
        used_labels = {"a1", "a2"}
        label = exporter._get_unique_site_label_with_base("a1", used_labels)
        self.assertEqual(label, "a3")  # Should increment

        # Test with multiple conflicts
        used_labels = {"a1", "a2", "a3"}
        label = exporter._get_unique_site_label_with_base("a1", used_labels)
        self.assertEqual(label, "a4")

    def test_write_mol_file(self):
        """Test _write_mol_file method."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        mol_type = self.system.molecule_types[0]
        mol_file_path = exporter._write_mol_file(mol_type)

        # Check file was created
        expected_path = exporter.output_dir / f"{mol_type.name}.mol"
        self.assertEqual(mol_file_path, expected_path)
        self.assertTrue(mol_file_path.exists())

        # Check file content
        content = mol_file_path.read_text()
        self.assertIn(f"Name = {mol_type.name}", content)
        self.assertIn("# translational diffusion constants", content)
        self.assertIn("# rotational diffusion constants", content)
        self.assertIn("# Coordinates", content)
        self.assertIn("COM", content)
        self.assertIn("# bonds", content)

        # Check that interface mapping was created
        self.assertGreater(len(exporter.interface_to_site_map), 0)
        self.assertGreater(len(exporter.homotypic_interface_map), 0)

    def test_write_mol_file_no_instances(self):
        """Test _write_mol_file with no molecule instances."""
        # Create system with no instances
        system_no_instances = Mock(spec=System)
        system_no_instances.molecule_types = [self.system.molecule_types[0]]
        system_no_instances.molecule_instances = []

        exporter = NERDSSExporter(system_no_instances, self.workspace_manager)

        mol_type = system_no_instances.molecule_types[0]
        mol_file_path = exporter._write_mol_file(mol_type)

        # Should still create file path but log warning
        expected_path = exporter.output_dir / f"{mol_type.name}.mol"
        self.assertEqual(mol_file_path, expected_path)

        # Check that warning was logged
        exporter.workspace_manager.logger.warning.assert_called()

    def test_generate_reactions_homotypic(self):
        """Test _generate_reactions for homotypic case."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # First write mol file to build mappings
        mol_type = self.system.molecule_types[0]
        exporter._write_mol_file(mol_type)

        # Generate reactions
        reactions = exporter._generate_reactions()

        # Should generate homotypic reactions
        self.assertGreater(len(reactions), 0)

        # Check reaction format
        for reaction in reactions:
            self.assertIn("<->", reaction)
            self.assertIn("!1", reaction)

        # Check metadata was created
        self.assertEqual(len(exporter.reaction_metadata), len(reactions))

    def test_format_vec(self):
        """Test _format_vec method."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Test with default precision
        result = exporter._format_vec([1.0, 2.5, -3.14159])
        expected = "1.0000000   2.5000000   -3.1415900"
        self.assertEqual(result, expected)

        # Test with custom precision
        result = exporter._format_vec([1.0, 2.5, -3.14159], precision=2)
        expected = "1.00   2.50   -3.14"
        self.assertEqual(result, expected)

        # Test with numpy array
        vec = np.array([1.0, 2.0, 3.0])
        result = exporter._format_vec(vec, precision=1)
        expected = "1.0   2.0   3.0"
        self.assertEqual(result, expected)

    def test_get_coms_interfaces(self):
        """Test _get_coms_interfaces method."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Set up interface mapping
        exporter.interface_to_site_map = {
            "A_A_1": "a1",
            "A_A_2": "a2"
        }

        # Test getting COMs and interfaces
        com1, com2, coord1, coord2 = exporter._get_coms_interfaces(
            "A", "a1", "A", "a2")

        # Should return valid coordinates
        self.assertIsNotNone(com1)
        self.assertIsNotNone(com2)
        self.assertIsNotNone(coord1)
        self.assertIsNotNone(coord2)

    def test_get_coms_interfaces_not_found(self):
        """Test _get_coms_interfaces with non-existent sites."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Test with empty mapping
        com1, com2, coord1, coord2 = exporter._get_coms_interfaces(
            "A", "nonexistent", "A", "also_nonexistent")

        # Should return None values
        self.assertIsNone(com1)
        self.assertIsNone(com2)
        self.assertIsNone(coord1)
        self.assertIsNone(coord2)

    @patch('ionerdss.model.pdb.nerdss_exporter.angles_from_points')
    @patch('ionerdss.model.pdb.nerdss_exporter.dihedrals_from_points')
    def test_generate_reaction_angles(self, mock_dihedrals, mock_angles):
        """Test _generate_reaction_angles method."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Mock angle functions
        mock_angles.return_value = 1.5
        mock_dihedrals.return_value = 0.5

        intf1 = np.array([0.0, 0.0, 0.0])
        intf2 = np.array([1.0, 0.0, 0.0])
        com1 = np.array([0.0, 1.0, 0.0])
        com2 = np.array([1.0, 1.0, 0.0])

        bond_length, angles = exporter._generate_reaction_angles(
            intf1, intf2, com1, com2)

        # Check bond length calculation
        expected_length = 1.0
        self.assertAlmostEqual(bond_length, expected_length, places=6)

        # Check angles tuple
        self.assertEqual(len(angles), 5)
        # All angles should be 1.5 (theta1, theta2) or 0.5 (phi1, phi2, omega)
        for angle in angles[:2]:
            self.assertEqual(angle, 1.5)
        for angle in angles[2:]:
            self.assertEqual(angle, 0.5)

    def test_calculate_reaction_parameters_caching(self):
        """Test _calculate_reaction_parameters with caching."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Set up homotypic mapping
        exporter.homotypic_interface_map = {
            "a1": "a1",
            "a2": "a1",  # Maps to same representative
        }

        reactions = [
            "A(a1) + A(a1) <-> A(a1!1).A(a1!1)",
            "A(a2) + A(a2) <-> A(a2!1).A(a2!1)"  # Should use cached parameters
        ]

        # Mock the coordinate lookup to return None (will use defaults)
        with patch.object(exporter, '_get_coms_interfaces', return_value=(None, None, None, None)):
            sigma_list, angles_list = exporter._calculate_reaction_parameters(
                reactions)

        # Should have parameters for both reactions
        self.assertEqual(len(sigma_list), 2)
        self.assertEqual(len(angles_list), 2)

        # Both should use same cached values (defaults in this case)
        self.assertEqual(sigma_list[0], sigma_list[1])
        self.assertEqual(angles_list[0], angles_list[1])

    def test_write_parms_file(self):
        """Test _write_parms_file method."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Set up reaction metadata
        exporter.reaction_metadata = [
            {'is_cross_reaction': False},
            {'is_cross_reaction': True}  # Cross-reaction should get doubled rate
        ]

        reactions = [
            "A(a1) + A(a1) <-> A(a1!1).A(a1!1)",
            "A(a1) + A(a2) <-> A(a1!1).A(a2!1)"
        ]
        molecule_counts = {"A": 10}
        box_nm = (100.0, 100.0, 100.0)
        sigma_list = [1.0, 1.5]
        angles_list = [(0.1, 0.2, 0.3, 0.4, 0.5), (0.6, 0.7, 0.8, 0.9, 1.0)]

        parms_path = exporter._write_parms_file(
            reactions, molecule_counts, box_nm, sigma_list, angles_list
        )

        # Check file was created
        expected_path = exporter.output_dir / "parms.inp"
        self.assertEqual(parms_path, expected_path)
        self.assertTrue(parms_path.exists())

        # Check file content
        content = parms_path.read_text()
        self.assertIn("start parameters", content)
        self.assertIn("start boundaries", content)
        self.assertIn("start molecules", content)
        self.assertIn("start reactions", content)
        self.assertIn("A : 10", content)
        self.assertIn("WaterBox = [100.0, 100.0, 100.0]", content)

        # Check that cross-reaction gets doubled rate
        lines = content.split('\n')
        on_rate_lines = [line for line in lines if 'onRate3Dka' in line]
        self.assertEqual(len(on_rate_lines), 2)

        # First reaction should have base rate (100.0)
        self.assertIn("onRate3Dka = 100.0", on_rate_lines[0])
        # Second reaction should have doubled rate (200.0)
        self.assertIn("onRate3Dka = 200.0", on_rate_lines[1])

    def test_write_parms_file_with_overrides(self):
        """Test _write_parms_file with parameter overrides."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        reactions = ["A(a1) + A(a1) <-> A(a1!1).A(a1!1)"]
        molecule_counts = {"A": 5}
        box_nm = (50.0, 50.0, 50.0)
        sigma_list = [2.0]
        angles_list = [(1.0, 1.0, 1.0, 1.0, 1.0)]

        overrides = {
            'nItr': 2e5,
            'timestep': 1.0,
            'onRate3Dka': 50.0
        }

        parms_path = exporter._write_parms_file(
            reactions, molecule_counts, box_nm, sigma_list, angles_list, overrides
        )

        # Check overrides were applied
        content = parms_path.read_text()
        self.assertIn("nItr = 200000.0", content)
        self.assertIn("timestep = 1.0", content)
        self.assertIn("onRate3Dka = 50.0", content)

    def test_export_all(self):
        """Test complete export_all workflow."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Mock angle calculation functions to avoid complex setup
        with patch('ionerdss.model.pdb.nerdss_exporter.angles_from_points', return_value=1.0):
            with patch('ionerdss.model.pdb.nerdss_exporter.dihedrals_from_points', return_value=0.5):
                output_files = exporter.export_all()

        # Check that files were created
        self.assertIn('A_mol', output_files)
        self.assertIn('parms', output_files)

        # Check files exist
        for file_path in output_files.values():
            self.assertTrue(file_path.exists())

        # Check mol file content
        mol_file = output_files['A_mol']
        mol_content = mol_file.read_text()
        self.assertIn("Name = A", mol_content)

        # Check parms file content
        parms_file = output_files['parms']
        parms_content = parms_file.read_text()
        self.assertIn("A : 10", parms_content)  # Default count

    def test_export_all_custom_parameters(self):
        """Test export_all with custom parameters."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        molecule_counts = {"A": 20}
        box_nm = (200.0, 200.0, 200.0)
        parms_overrides = {"timestep": 0.1}

        with patch('ionerdss.model.pdb.nerdss_exporter.angles_from_points', return_value=1.0):
            with patch('ionerdss.model.pdb.nerdss_exporter.dihedrals_from_points', return_value=0.5):
                output_files = exporter.export_all(
                    molecule_counts=molecule_counts,
                    box_nm=box_nm,
                    parms_overrides=parms_overrides
                )

        # Check custom parameters were used
        parms_file = output_files['parms']
        parms_content = parms_file.read_text()
        self.assertIn("A : 20", parms_content)
        self.assertIn("WaterBox = [200.0, 200.0, 200.0]", parms_content)
        self.assertIn("timestep = 0.1", parms_content)

    def test_homotypic_mapping_consistency(self):
        """Test that homotypic mapping ensures parameter consistency."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Write mol file to build mappings
        mol_type = self.system.molecule_types[0]
        exporter._write_mol_file(mol_type)

        # Check that homotypic mapping was created
        self.assertGreater(len(exporter.homotypic_interface_map), 0)

        # All sites of the same type should map to the same representative
        representatives = set(exporter.homotypic_interface_map.values())
        self.assertGreater(len(representatives), 0)

        # Each representative should map to itself
        for site, rep in exporter.homotypic_interface_map.items():
            if site == rep:  # This is a representative
                self.assertIn(rep, exporter.homotypic_interface_map)

    def test_reaction_metadata_cross_reaction_detection(self):
        """Test that cross-reactions are properly detected."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Write mol file and generate reactions
        mol_type = self.system.molecule_types[0]
        exporter._write_mol_file(mol_type)
        reactions = exporter._generate_reactions()

        # Check that metadata was created
        self.assertEqual(len(exporter.reaction_metadata), len(reactions))

        # Check for cross-reaction detection
        has_cross_reaction = any(meta['is_cross_reaction']
                                 for meta in exporter.reaction_metadata)
        has_non_cross_reaction = any(
            not meta['is_cross_reaction'] for meta in exporter.reaction_metadata)

        # Should have both types if multiple sites exist
        if len(exporter.interface_to_site_map) > 1:
            self.assertTrue(has_cross_reaction or has_non_cross_reaction)


class TestNERDSSExporterIntegration(unittest.TestCase):
    """Integration tests for NERDSSExporter."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)

    def tearDown(self):
        """Clean up test fixtures."""
        if self.temp_path.exists():
            shutil.rmtree(self.temp_path)

    def test_complete_export_workflow(self):
        """Test complete export workflow with realistic system."""
        # Create a more complex system
        system = self._create_complex_system()

        # Create workspace manager
        workspace_manager = Mock(spec=WorkspaceManager)
        workspace_manager.workspace_path = self.temp_path
        workspace_manager.logger = Mock()

        # Create exporter
        exporter = NERDSSExporter(system, workspace_manager)

        # Mock angle calculations
        with patch('ionerdss.model.pdb.nerdss_exporter.angles_from_points', return_value=1.2):
            with patch('ionerdss.model.pdb.nerdss_exporter.dihedrals_from_points', return_value=0.8):
                output_files = exporter.export_all(
                    molecule_counts={"ProteinA": 15, "ProteinB": 10},
                    box_nm=(150.0, 150.0, 150.0)
                )

        # Verify all expected files were created
        expected_files = ['ProteinA_mol', 'ProteinB_mol', 'parms']
        for file_key in expected_files:
            self.assertIn(file_key, output_files)
            self.assertTrue(output_files[file_key].exists())

        # Verify file contents
        parms_content = output_files['parms'].read_text()
        self.assertIn("ProteinA : 15", parms_content)
        self.assertIn("ProteinB : 10", parms_content)
        self.assertIn("WaterBox = [150.0, 150.0, 150.0]", parms_content)

    def _create_complex_system(self):
        """Create a more complex system for integration testing."""
        system = Mock(spec=System)

        # Create two different molecule types
        mol_type_a = Mock(spec=MoleculeType)
        mol_type_a.name = "ProteinA"
        mol_type_a.D_t_nm2_us = 1.5
        mol_type_a.D_r_rad2_us = 0.15

        mol_type_b = Mock(spec=MoleculeType)
        mol_type_b.name = "ProteinB"
        mol_type_b.D_t_nm2_us = 1.2
        mol_type_b.D_r_rad2_us = 0.12

        # Create interface types for heterotypic interaction
        interface_type_ab = Mock(spec=InterfaceType)
        interface_type_ab.get_name.return_value = "ProteinA_ProteinB_1"
        interface_type_ab.this_mol_type_name = "ProteinA"
        interface_type_ab.partner_mol_type_name = "ProteinB"
        interface_type_ab.interface_index = 1
        interface_type_ab.local_coord = np.array([2.0, 0.0, 0.0])
        interface_type_ab.absolute_coord = np.array([2.0, 0.0, 0.0])

        interface_type_ba = Mock(spec=InterfaceType)
        interface_type_ba.get_name.return_value = "ProteinB_ProteinA_1"
        interface_type_ba.this_mol_type_name = "ProteinB"
        interface_type_ba.partner_mol_type_name = "ProteinA"
        interface_type_ba.interface_index = 1
        interface_type_ba.local_coord = np.array([-2.0, 0.0, 0.0])
        interface_type_ba.absolute_coord = np.array([-2.0, 0.0, 0.0])

        # Create molecule instances
        mol_instance_a = Mock(spec=MoleculeInstance)
        mol_instance_a.molecule_type = mol_type_a
        mol_instance_a.com = np.array([0.0, 0.0, 0.0])

        mol_instance_b = Mock(spec=MoleculeInstance)
        mol_instance_b.molecule_type = mol_type_b
        mol_instance_b.com = np.array([4.0, 0.0, 0.0])

        # Create interface instances
        interface_instance_ab = Mock(spec=InterfaceInstance)
        interface_instance_ab.interface_type = interface_type_ab
        interface_instance_ab.absolute_coord = np.array([2.0, 0.0, 0.0])

        interface_instance_ba = Mock(spec=InterfaceInstance)
        interface_instance_ba.interface_type = interface_type_ba
        interface_instance_ba.absolute_coord = np.array([2.0, 0.0, 0.0])

        # Set up partner interfaces
        interface_instance_ab.partner_interface = interface_instance_ba
        interface_instance_ba.partner_interface = interface_instance_ab

        # Set up interfaces_neighbors_map
        mol_instance_a.interfaces_neighbors_map = {
            interface_instance_ab: mol_instance_b}
        mol_instance_b.interfaces_neighbors_map = {
            interface_instance_ba: mol_instance_a}

        # Set up system
        system.molecule_types = [mol_type_a, mol_type_b]
        system.interface_types = [interface_type_ab, interface_type_ba]
        system.molecule_instances = [mol_instance_a, mol_instance_b]

        return system


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

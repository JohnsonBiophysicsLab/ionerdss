"""
Unit tests for ionerdss.model.pdb.nerdss_exporter

Complete rewrite to match current API and implementation.
Tests the NERDSSExporter class and its file generation capabilities.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, MagicMock, patch
import numpy as np

from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.components.instances import MoleculeInstance, InterfaceInstance
from ionerdss.model.pdb.file_manager import WorkspaceManager


class TestNERDSSExporterBasics(unittest.TestCase):
    """Test basic NERDSSExporter functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)

        # Create minimal system
        self.system = Mock(spec=System)
        self.system.molecule_types = []
        self.system.interface_types = []
        self.system.molecule_instances = []

        # Create workspace manager
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.workspace_path = self.temp_path
        self.workspace_manager.logger = Mock()

    def tearDown(self):
        """Clean up test fixtures."""
        if self.temp_path.exists():
            shutil.rmtree(self.temp_path)

    def test_initialization_with_workspace_manager(self):
        """Test NERDSSExporter initialization with workspace manager."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        self.assertEqual(exporter.system, self.system)
        self.assertEqual(exporter.workspace_manager, self.workspace_manager)

        # Check output directory creation
        expected_output_dir = self.temp_path / 'nerdss_files'
        self.assertEqual(exporter.output_dir, expected_output_dir)
        self.assertTrue(expected_output_dir.exists())

    def test_initialization_without_workspace_manager(self):
        """Test NERDSSExporter initialization without workspace manager."""
        exporter = NERDSSExporter(self.system)

        self.assertEqual(exporter.system, self.system)
        self.assertIsNone(exporter.workspace_manager)

        # Should create default output directory
        expected_output_dir = Path("nerdss_files")
        self.assertEqual(exporter.output_dir, expected_output_dir)

    def test_format_vec(self):
        """Test _format_vec method."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Test with default precision
        result = exporter._format_vec([1.0, 2.5, -3.14159])
        self.assertIn("1.0000000", result)
        self.assertIn("2.5000000", result)
        self.assertIn("-3.1415900", result)

        # Test with numpy array
        vec = np.array([1.0, 2.0, 3.0])
        result = exporter._format_vec(vec, precision=1)
        expected = "1.0   2.0   3.0"
        self.assertEqual(result, expected)

    def test_circular_mean_signed(self):
        """Test circular mean calculation."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Test with angles around 0
        angles = [0.1, -0.1, 0.05]
        mean = exporter._circular_mean_signed(angles)
        self.assertAlmostEqual(mean, 0.0166, places=3)

        # Test with angles around π
        angles = [np.pi - 0.1, -np.pi + 0.1]
        mean = exporter._circular_mean_signed(angles)
        self.assertAlmostEqual(abs(mean), np.pi, places=1)

    def test_circular_var_signed(self):
        """Test circular variance calculation."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Test with very similar angles (low variance)
        angles = [0.0, 0.01, -0.01]
        var = exporter._circular_var_signed(angles)
        self.assertLess(var, 0.01)

        # Test with spread angles (higher variance)
        angles = [-1.0, 0.0, 1.0]
        var = exporter._circular_var_signed(angles)
        self.assertGreater(var, 0.01)


class TestNERDSSExporterWithMolecules(unittest.TestCase):
    """Test NERDSSExporter with realistic molecule setup."""

    def setUp(self):
        """Set up test fixtures with molecules."""
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)

        # Create workspace manager
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.workspace_path = self.temp_path
        self.workspace_manager.logger = Mock()

        # Create real molecule type
        self.mol_type_a = MoleculeType(
            name="A",
            radius_nm=2.0,
            D_t_nm2_us=1.0,
            D_r_rad2_us=0.1,
            ref1_local=np.array([1.0, 0.0, 0.0]),
            ref2_local=np.array([0.0, 0.0, 1.0])
        )

        # Create interface types
        self.interface_type_1 = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="A",
            interface_index=1,
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            local_coord=np.array([0.5, 0.0, 0.0]),
            energy=-5.0,
            this_mol_type=self.mol_type_a,
            partner_mol_type=self.mol_type_a
        )

        # Create molecule instances
        self.mol_instance_1 = MoleculeInstance(
            name="A_1",
            molecule_type=self.mol_type_a,
            com=np.array([0.0, 0.0, 0.0]),
            norm=np.array([0.0, 0.0, 1.0]),
            ref1=np.array([1.0, 0.0, 0.0]),
            ref2=np.array([0.0, 0.0, 1.0])
        )

        self.mol_instance_2 = MoleculeInstance(
            name="A_2",
            molecule_type=self.mol_type_a,
            com=np.array([5.0, 0.0, 0.0]),
            norm=np.array([0.0, 0.0, 1.0]),
            ref1=np.array([1.0, 0.0, 0.0]),
            ref2=np.array([0.0, 0.0, 1.0])
        )

        # Create interface instances
        self.interface_instance_1 = InterfaceInstance(
            interface_type=self.interface_type_1,
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            this_mol=self.mol_instance_1,
            this_mol_name="A",
            partner_mol_name="A",
            interface_index=1
        )

        self.interface_instance_2 = InterfaceInstance(
            interface_type=self.interface_type_1,
            absolute_coord=np.array([4.0, 0.0, 0.0]),
            this_mol=self.mol_instance_2,
            this_mol_name="A",
            partner_mol_name="A",
            interface_index=1
        )

        # Set up bindings
        self.interface_instance_1.partner_interface = self.interface_instance_2
        self.interface_instance_2.partner_interface = self.interface_instance_1

        self.mol_instance_1.interfaces_neighbors_map = {
            self.interface_instance_1: self.mol_instance_2
        }
        self.mol_instance_2.interfaces_neighbors_map = {
            self.interface_instance_2: self.mol_instance_1
        }

        # Create system (use mock to avoid workspace_path requirement)
        self.system = Mock(spec=System)
        self.system.molecule_types = [self.mol_type_a]
        self.system.interface_types = [self.interface_type_1]
        self.system.molecule_instances = [self.mol_instance_1, self.mol_instance_2]

    def tearDown(self):
        """Clean up test fixtures."""
        if self.temp_path.exists():
            shutil.rmtree(self.temp_path)

    def test_write_mol_file_creates_file(self):
        """Test that _write_mol_file creates a .mol file."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        mol_file_path = exporter._write_mol_file(self.mol_type_a)

        # Check file was created
        self.assertTrue(mol_file_path.exists())
        self.assertEqual(mol_file_path.suffix, ".mol")
        self.assertEqual(mol_file_path.stem, "A")

    def test_write_mol_file_content(self):
        """Test that .mol file contains expected content."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        mol_file_path = exporter._write_mol_file(self.mol_type_a)
        content = mol_file_path.read_text()

        # Check essential sections
        self.assertIn("Name = A", content)
        self.assertIn("# translational diffusion constants", content)
        self.assertIn("# rotational diffusion constants", content)
        self.assertIn("# Coordinates", content)
        self.assertIn("COM", content)
        self.assertIn("# bonds", content)

    def test_get_base_site_label(self):
        """Test _get_base_site_label method."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        # Current implementation includes partner initial
        label = exporter._get_base_site_label("A", "A_A_1")
        # Expected format: first letter of mol + first letter of partner + index
        self.assertTrue(label.startswith("a"))  # First letter lowercase
        self.assertIn("1", label)  # Contains index

    def test_write_parms_file_creates_file(self):
        """Test that _write_parms_file creates parms.inp."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        reactions = ["A(a1) + A(a1) <-> A(a1!1).A(a1!1)"]
        molecule_counts = {"A": 10}
        box_nm = (100.0, 100.0, 100.0)
        sigma_list = [1.0]
        angles_list = [(1.57, 1.57, 0.0, 0.0, 0.0)]

        parms_path = exporter._write_parms_file(
            reactions, molecule_counts, box_nm, sigma_list, angles_list
        )

        # Check file was created
        self.assertTrue(parms_path.exists())
        self.assertEqual(parms_path.name, "parms.inp")

    def test_write_parms_file_content(self):
        """Test that parms.inp contains expected content."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        reactions = ["A(a1) + A(a1) <-> A(a1!1).A(a1!1)"]
        molecule_counts = {"A": 15}
        box_nm = (150.0, 150.0, 150.0)
        sigma_list = [1.5]
        angles_list = [(1.0, 1.0, 0.5, 0.5, 0.0)]

        parms_path = exporter._write_parms_file(
            reactions, molecule_counts, box_nm, sigma_list, angles_list
        )

        content = parms_path.read_text()

        # Check essential sections
        self.assertIn("start parameters", content)
        self.assertIn("start boundaries", content)
        self.assertIn("start molecules", content)
        self.assertIn("start reactions", content)
        self.assertIn("A : 15", content)
        self.assertIn("WaterBox = [150.0, 150.0, 150.0]", content)

    def test_export_all_generates_files(self):
        """Test that export_all generates all required files."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        output_files = exporter.export_all(
            molecule_counts={"A": 20},
            box_nm=(200.0, 200.0, 200.0)
        )

        # Check that files were created
        self.assertIn('A_mol', output_files)
        self.assertIn('parms', output_files)

        # Check files exist
        for file_key, file_path in output_files.items():
            self.assertTrue(file_path.exists(), f"{file_key} file should exist")

    def test_export_all_with_custom_parameters(self):
        """Test export_all with custom parameters."""
        exporter = NERDSSExporter(self.system, self.workspace_manager)

        output_files = exporter.export_all(
            molecule_counts={"A": 25},
            box_nm=(250.0, 250.0, 250.0),
            parms_overrides={"timestep": 0.5, "nItr": 1e5}
        )

        # Check parms file contains overrides
        parms_content = output_files['parms'].read_text()
        self.assertIn("timestep = 0.5", parms_content)
        self.assertIn("nItr = 100000.0", parms_content)
        self.assertIn("A : 25", parms_content)
        self.assertIn("WaterBox = [250.0, 250.0, 250.0]", parms_content)


class TestAngleCalculation(unittest.TestCase):
    """Test angle calculation methods."""

    def setUp(self):
        """Set up test fixtures."""
        self.system = Mock(spec=System)
        self.system.molecule_types = []
        self.system.interface_types = []
        self.system.molecule_instances = []

        self.exporter = NERDSSExporter(self.system)

    def test_generate_reaction_angles_signature(self):
        """Test that _generate_reaction_angles has expected signature."""
        # The method signature changed to require mol names and site labels
        # This is a minimal test to ensure the method exists
        import inspect
        
        sig = inspect.signature(self.exporter._generate_reaction_angles)
        params = list(sig.parameters.keys())
        
        # Check expected parameters exist
        self.assertIn('intf1', params)
        self.assertIn('intf2', params)
        self.assertIn('com1', params)
        self.assertIn('com2', params)
        self.assertIn('mol1_name', params)
        self.assertIn('mol2_name', params)
        self.assertIn('site1', params)
        self.assertIn('site2', params)


class TestNormalVectorCalculation(unittest.TestCase):
    """Test normal vector calculation methods."""

    def setUp(self):
        """Set up test system with molecules."""
        # Create minimal system for normal vector testing
        self.mol_type = MoleculeType(
            name="TestMol",
            radius_nm=1.0,
            ref1_local=np.array([1.0, 0.0, 0.0]),
            ref2_local=np.array([0.0, 0.0, 1.0])
        )

        self.mol_instance = MoleculeInstance(
            name="TestMol_1",
            molecule_type=self.mol_type,
            com=np.array([0.0, 0.0, 0.0]),
            norm=np.array([0.0, 0.0, 1.0]),
            ref1=np.array([1.0, 0.0, 0.0]),
            ref2=np.array([0.0, 0.0, 1.0])
        )

        self.system = Mock(spec=System)
        self.system.molecule_types = [self.mol_type]
        self.system.molecule_instances = [self.mol_instance]
        self.system.interface_types = []

        self.exporter = NERDSSExporter(self.system)

    def test_calculate_normal_vectors_runs(self):
        """Test that _calculate_normal_vectors executes without error."""
        # Should not raise an exception
        try:
            self.exporter._calculate_normal_vectors()
            success = True
        except Exception as e:
            success = False
            print(f"Error: {e}")

        self.assertTrue(success)


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

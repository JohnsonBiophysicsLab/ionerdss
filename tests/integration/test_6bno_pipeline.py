"""
Integration tests for the complete ionerdss pipeline.

Tests the full workflow from PDB/CIF parsing through NERDSS file generation,
using the 6bno structure as a realistic test case.
"""

import unittest
import tempfile
import shutil
import json
from pathlib import Path
import numpy as np

from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter


class Test6BNOPipelineIntegration(unittest.TestCase):
    """Integration test for complete 6bno processing pipeline."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures that are expensive to create."""
        # Use the local test data copy
        # This ensures tests run in clean environments (like CI)
        project_root = Path(__file__).parent.parent.parent
        cls.cif_path = project_root / "tests/data/test_6BNO.cif"
        
        if not cls.cif_path.exists():
            raise FileNotFoundError(
                f"test_6BNO.cif not found at {cls.cif_path}. "
            )
        
        # Create temporary workspace for testing
        cls.temp_dir = tempfile.mkdtemp(prefix="ionerdss_test_")
        cls.workspace_path = Path(cls.temp_dir) / "test_6bno_workspace"

    @classmethod
    def tearDownClass(cls):
        """Clean up temporary files."""
        if Path(cls.temp_dir).exists():
            shutil.rmtree(cls.temp_dir)

    def build_test_system(self):
        """Helper to build a system for testing."""
        model = PDBModelBuilder(source=str(self.cif_path))
        
        hyperparams = PDBModelHyperparameters(
            interface_detect_distance_cutoff=1.0,
            ring_regularization_mode="off"
        )
        
        system = model.build_system(
            workspace_path=str(self.workspace_path),
            hyperparams=hyperparams
        )
        
        return model, system

    def test_01_pdb_model_builder_initialization(self):
        """Test that PDBModelBuilder can initialize with 6bno CIF file."""
        model = PDBModelBuilder(source=str(self.cif_path))
        
        # Verify model was created
        self.assertIsNotNone(model)
        self.assertIsNotNone(model.source)

    def test_02_build_system_with_hyperparameters(self):
        """Test building the complete system with custom hyperparameters."""
        model, system = self.build_test_system()
        
        # Verify workspace was created
        self.assertTrue(self.workspace_path.exists())
        
        # Verify system was created
        self.assertIsNotNone(system)

    def test_03_system_structure_validation(self):
        """Test that the generated system has expected structure."""
        model, system = self.build_test_system()
        
        # Verify system has molecules
        self.assertGreater(len(system.molecule_types), 0,
                          "System should have at least one molecule type")
        
        # Verify system has instances
        self.assertGreater(len(system.molecule_instances), 0,
                          "System should have molecule instances")
        
        # Verify molecule types have required attributes
        for mol_type in system.molecule_types:
            self.assertIsNotNone(mol_type.name)
            self.assertGreater(mol_type.radius_nm, 0)
            self.assertIsNotNone(mol_type.ref1_local)
            self.assertIsNotNone(mol_type.ref2_local)
            
            # Verify reference vectors are unit vectors
            ref1_mag = np.linalg.norm(mol_type.ref1_local)
            ref2_mag = np.linalg.norm(mol_type.ref2_local)
            self.assertAlmostEqual(ref1_mag, 1.0, places=5,
                                  msg=f"ref1 should be unit vector for {mol_type.name}")
            self.assertAlmostEqual(ref2_mag, 1.0, places=5,
                                  msg=f"ref2 should be unit vector for {mol_type.name}")
            
            # Verify reference vectors are orthogonal
            dot_product = np.dot(mol_type.ref1_local, mol_type.ref2_local)
            self.assertAlmostEqual(dot_product, 0.0, places=5,
                                  msg=f"ref1 and ref2 should be orthogonal for {mol_type.name}")

    def test_04_molecule_instances_have_coordinates(self):
        """Test that molecule instances have valid coordinates."""
        model, system = self.build_test_system()
        
        # Check molecule instances
        for mol_instance in system.molecule_instances:
            # Verify COM exists and is valid
            self.assertIsNotNone(mol_instance.com)
            self.assertEqual(len(mol_instance.com), 3)
            self.assertFalse(np.any(np.isnan(mol_instance.com)),
                           f"COM should not contain NaN for {mol_instance.name}")
            
            # Verify reference vectors exist
            self.assertIsNotNone(mol_instance.ref1)
            self.assertIsNotNone(mol_instance.ref2)
            self.assertEqual(len(mol_instance.ref1), 3)
            self.assertEqual(len(mol_instance.ref2), 3)
            
            # Verify norm vector exists
            self.assertIsNotNone(mol_instance.norm)
            self.assertEqual(len(mol_instance.norm), 3)

    def test_05_interface_detection(self):
        """Test that interfaces are detected between molecules."""
        model, system = self.build_test_system()
        
        # Verify interfaces were detected
        self.assertGreater(len(system.interface_types), 0,
                          "System should have interface types")
        
        # Check interface properties
        for intf_type in system.interface_types:
            self.assertIsNotNone(intf_type.this_mol_type_name)
            self.assertIsNotNone(intf_type.partner_mol_type_name)
            self.assertGreaterEqual(intf_type.interface_index, 1)
            
            # Verify coordinates
            self.assertIsNotNone(intf_type.absolute_coord)
            self.assertIsNotNone(intf_type.local_coord)
            self.assertEqual(len(intf_type.absolute_coord), 3)
            self.assertEqual(len(intf_type.local_coord), 3)

    def test_06_nerdss_export_generates_files(self):
        """Test that NERDSS exporter generates required files."""
        model, system = self.build_test_system()
        workspace_manager = model.workspace_manager
        
        # Create exporter
        exporter = NERDSSExporter(system, workspace_manager)
        
        # Export NERDSS files
        molecule_counts = {mol_type.name: 10 for mol_type in system.molecule_types}
        
        output_files = exporter.export_all(
            molecule_counts=molecule_counts,
            box_nm=(100.0, 100.0, 100.0)
        )
        
        # Verify files were generated
        self.assertIn('parms', output_files)
        self.assertTrue(output_files['parms'].exists(),
                       "parms.inp should be created")
        
        # Verify .mol files for each molecule type
        for mol_type in system.molecule_types:
            mol_file_key = f'{mol_type.name}_mol'
            self.assertIn(mol_file_key, output_files,
                         f"{mol_type.name}.mol should be in output")
            self.assertTrue(output_files[mol_file_key].exists(),
                           f"{mol_type.name}.mol should exist")

    def test_07_parms_file_content_validation(self):
        """Test that parms.inp contains expected content and structure."""
        model, system = self.build_test_system()
        workspace_manager = model.workspace_manager
        
        # Export
        exporter = NERDSSExporter(system, workspace_manager)
        molecule_counts = {mol_type.name: 15 for mol_type in system.molecule_types}
        
        output_files = exporter.export_all(
            molecule_counts=molecule_counts,
            box_nm=(150.0, 150.0, 150.0)
        )
        
        # Read parms.inp
        parms_content = output_files['parms'].read_text()
        
        # Verify essential sections
        self.assertIn("start parameters", parms_content)
        self.assertIn("end parameters", parms_content)
        self.assertIn("start boundaries", parms_content)
        self.assertIn("end boundaries", parms_content)
        self.assertIn("start molecules", parms_content)
        self.assertIn("end molecules", parms_content)
        self.assertIn("start reactions", parms_content)
        self.assertIn("end reactions", parms_content)
        
        # Verify box dimensions
        self.assertIn("WaterBox = [150.0, 150.0, 150.0]", parms_content)
        
        # Verify molecule counts
        for mol_type in system.molecule_types:
            self.assertIn(f"{mol_type.name} : 15", parms_content)

    def test_08_mol_file_content_validation(self):
        """Test that .mol files contain expected structure."""
        model, system = self.build_test_system()
        workspace_manager = model.workspace_manager
        
        # Export
        exporter = NERDSSExporter(system, workspace_manager)
        molecule_counts = {mol_type.name: 10 for mol_type in system.molecule_types}
        
        output_files = exporter.export_all(
            molecule_counts=molecule_counts,
            box_nm=(100.0, 100.0, 100.0)
        )
        
        # Check each mol file
        for mol_type in system.molecule_types:
            mol_file_key = f'{mol_type.name}_mol'
            mol_content = output_files[mol_file_key].read_text()
            
            # Verify essential sections
            self.assertIn(f"Name = {mol_type.name}", mol_content)
            self.assertIn("# translational diffusion constants", mol_content)
            self.assertIn("# rotational diffusion constants", mol_content)
            self.assertIn("# Coordinates", mol_content)
            self.assertIn("COM", mol_content)
            
            # Verify diffusion constants are present and numeric
            lines = mol_content.split('\n')
            found_D_t = False
            found_D_r = False
            
            for line in lines:
                if line.strip().startswith("D ="):
                    found_D_t = True
                    # Extract and verify it's a number
                    # Find the first '[' and the first ',' after it
                    start_index = line  .find('[') + 1
                    end_index = line.find(',', start_index)
                    if start_index > 0 and end_index != -1:
                        value = line[start_index:end_index]
                    self.assertGreater(float(value), 0,
                                     f"D_t should be positive for {mol_type.name}")
                elif line.strip().startswith("Dr ="):
                    found_D_r = True
                    # Extract and verify it's a number
                    # Find the first '[' and the first ',' after it
                    start_index = line.find('[') + 1
                    end_index = line.find(',', start_index)
                    if start_index > 0 and end_index != -1:
                        value = line[start_index:end_index]
                    self.assertGreater(float(value), 0,
                                     f"D_r should be positive for {mol_type.name}")
            
            self.assertTrue(found_D_t, f"D_t not found in {mol_type.name}.mol")
            self.assertTrue(found_D_r, f"D_r not found in {mol_type.name}.mol")

    def test_09_reaction_angles_are_valid(self):
        """Test that generated reaction angles are within valid ranges."""
        model, system = self.build_test_system()
        workspace_manager = model.workspace_manager
        
        # Export
        exporter = NERDSSExporter(system, workspace_manager)
        molecule_counts = {mol_type.name: 10 for mol_type in system.molecule_types}
        
        output_files = exporter.export_all(
            molecule_counts=molecule_counts,
            box_nm=(100.0, 100.0, 100.0)
        )
        
        # Read parms.inp
        parms_content = output_files['parms'].read_text()
        lines = parms_content.split('\n')
        
        # Find reaction lines and extract angles
        in_reactions = False
        for line in lines:
            if "start reactions" in line:
                in_reactions = True
                continue
            if "end reactions" in line:
                in_reactions = False
                continue
            
            if in_reactions and "onRate3Dka" in line:
                # Extract angle values (theta1, theta2, phi1, phi2, omega)
                parts = line.split(',')
                for part in parts:
                    if 'theta' in part or 'phi' in part or 'omega' in part:
                        # Extract numeric value
                        value_str = part.split('=')[1].strip()
                        angle_value = float(value_str)
                        
                        # Verify angle is within valid range [-π, π]
                        self.assertGreaterEqual(angle_value, -np.pi - 0.01,
                                              f"Angle {angle_value} too negative")
                        self.assertLessEqual(angle_value, np.pi + 0.01,
                                           f"Angle {angle_value} too positive")
                        
                        # Verify angle is not NaN or Inf
                        self.assertFalse(np.isnan(angle_value),
                                       "Angle should not be NaN")
                        self.assertFalse(np.isinf(angle_value),
                                       "Angle should not be Inf")

    def test_10_coordinate_system_consistency(self):
        """Test that coordinate systems are consistent throughout pipeline."""
        model, system = self.build_test_system()
        
        # Verify all molecule types have consistent coordinate frames
        for mol_type in system.molecule_types:
            # Verify orthogonality
            dot = np.dot(mol_type.ref1_local, mol_type.ref2_local)
            self.assertAlmostEqual(dot, 0.0, places=6,
                                  msg=f"ref1⊥ref2 for {mol_type.name}")
            
            # Verify cross product gives a valid third axis
            ref3 = np.cross(mol_type.ref1_local, mol_type.ref2_local)
            ref3_mag = np.linalg.norm(ref3)
            self.assertAlmostEqual(ref3_mag, 1.0, places=5,
                                  msg=f"ref1×ref2 should be unit for {mol_type.name}")

    def test_11_end_to_end_pipeline_no_exceptions(self):
        """Test that complete pipeline runs without exceptions."""
        try:
            model, system = self.build_test_system()
            workspace_manager = model.workspace_manager
            
            # Export NERDSS files
            exporter = NERDSSExporter(system, workspace_manager)
            molecule_counts = {mol_type.name: 20 for mol_type in system.molecule_types}
            
            output_files = exporter.export_all(
                molecule_counts=molecule_counts,
                box_nm=(200.0, 200.0, 200.0)
            )
            
            # If we get here, success!
            self.assertTrue(True)
            self.assertGreater(len(output_files), 0)
            
        except Exception as e:
            self.fail(f"Pipeline raised unexpected exception: {e}")

    def test_12_output_files_are_readable(self):
        """Test that all output files can be read and parsed."""
        model, system = self.build_test_system()
        workspace_manager = model.workspace_manager
        
        exporter = NERDSSExporter(system, workspace_manager)
        molecule_counts = {mol_type.name: 5 for mol_type in system.molecule_types}
        
        output_files = exporter.export_all(
            molecule_counts=molecule_counts,
            box_nm=(50.0, 50.0, 50.0)
        )
        
        # Verify all files can be read
        for file_key, file_path in output_files.items():
            self.assertTrue(file_path.exists(), f"{file_key} file should exist")
            
            # Try to read content
            try:
                content = file_path.read_text()
                self.assertGreater(len(content), 0,
                                 f"{file_key} should not be empty")
            except Exception as e:
                self.fail(f"Failed to read {file_key}: {e}")


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

"""Unit tests for affinity prediction in PDBModelBuilder."""

import unittest
import sys
import os
import tempfile
import shutil

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters


class TestAffinityPrediction(unittest.TestCase):
    """Test suite for affinity prediction feature."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.pdb_id = "8erq"
    
    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
    
    def test_build_system_default_energy(self):
        """Test build_system with default fixed energy (predict_affinity=False)."""
        hyperparams = PDBModelHyperparameters(
            interface_detect_distance_cutoff=0.35,
            interface_detect_n_residue_cutoff=3,
            predict_affinity=False
        )
        
        model = PDBModelBuilder(source=self.pdb_id, hyperparams=hyperparams)
        system = model.build_system(
            workspace_path=os.path.join(self.test_dir, self.pdb_id),
            hyperparams=hyperparams
        )
        
        # Check that the system was built successfully
        self.assertIsNotNone(system)
        self.assertGreater(len(system.molecule_types), 0)
        
        # Check that interfaces were detected
        self.assertGreater(len(system.interface_types), 0)
    
    def test_hyperparameters_signature(self):
        """Test that hyperparameters include affinity prediction parameters."""
        hyperparams = PDBModelHyperparameters()
        
        # Check that new parameters exist
        self.assertTrue(hasattr(hyperparams, 'predict_affinity'))
        self.assertTrue(hasattr(hyperparams, 'adfr_path'))
        
        # Check default values
        self.assertEqual(hyperparams.predict_affinity, False)
        self.assertEqual(hyperparams.adfr_path, None)
    
    def test_interface_detection(self):
        """Test that interface detection works correctly."""
        hyperparams = PDBModelHyperparameters(
            interface_detect_distance_cutoff=0.35,
            interface_detect_n_residue_cutoff=3,
            predict_affinity=False
        )
        
        model = PDBModelBuilder(source=self.pdb_id, hyperparams=hyperparams)
        system = model.build_system(
            workspace_path=os.path.join(self.test_dir, self.pdb_id),
            hyperparams=hyperparams
        )
        
        # Check that system components are properly initialized
        summary = system.get_summary()
        self.assertGreater(summary['num_molecule_types'], 0, "No molecule types created")
        self.assertGreater(summary['num_interface_types'], 0, "No interface types detected")


class TestAffinityPredictionWithProAffinity(unittest.TestCase):
    """Test suite for ProAffinity integration (requires ADFR and model)."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.pdb_id = "8erq"
        # Path to ADFR - modify this if running tests
        self.adfr_path = os.environ.get('ADFR_PATH', None)
    
    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
    
    @unittest.skipIf(
        os.environ.get('ADFR_PATH') is None,
        "ADFR_PATH not set, skipping ProAffinity tests"
    )
    def test_build_system_with_prediction(self):
        """Test build_system with ProAffinity prediction enabled."""
        hyperparams = PDBModelHyperparameters(
            interface_detect_distance_cutoff=0.35,
            interface_detect_n_residue_cutoff=3,
            predict_affinity=True,
            adfr_path=self.adfr_path
        )
        
        model = PDBModelBuilder(source=self.pdb_id, hyperparams=hyperparams)
        system = model.build_system(
            workspace_path=os.path.join(self.test_dir, self.pdb_id),
            hyperparams=hyperparams
        )
        
        # Check that the system was built successfully
        self.assertIsNotNone(system)
        self.assertGreater(len(system.molecule_types), 0)
        
        # Check that interfaces were detected
        self.assertGreater(len(system.interface_types), 0)
        
        # With ProAffinity enabled, energies should be predicted
        # (actual values will vary based on structure)
        for interface_type in system.interface_types:
            # Check that binding energy is set
            # The value should be different from the default -16 RT
            self.assertIsNotNone(interface_type.sigma)


if __name__ == '__main__':
    unittest.main()

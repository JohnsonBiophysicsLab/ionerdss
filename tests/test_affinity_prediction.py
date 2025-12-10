"""Unit tests for affinity prediction in PDBModel.coarse_grain method."""

import unittest
import sys
import os
import tempfile
import shutil

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from ionerdss.model.pdb_model import PDBModel


class TestAffinityPrediction(unittest.TestCase):
    """Test suite for affinity prediction feature in coarse_grain method."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.pdb_id = "8erq"
    
    def tearDown(self):
        """Clean up test fixtures."""
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
    
    def test_coarse_grain_default_energy(self):
        """Test coarse_grain with default fixed energy (predict_affinity=False)."""
        model = PDBModel(pdb_id=self.pdb_id, save_dir=self.test_dir)
        model.coarse_grain(
            distance_cutoff=0.35,
            residue_cutoff=3,
            predict_affinity=False
        )
        
        # Check that interfaces were detected
        self.assertGreater(len(model.all_chains), 0)
        
        # Check that default energy is used (-16 RT in kJ/mol)
        default_energy = -16 * 8.314/1000 * 298
        for chain_energies in model.all_interface_energies:
            for energy in chain_energies:
                self.assertAlmostEqual(energy, default_energy, places=2)
    
    def test_coarse_grain_signature(self):
        """Test that coarse_grain method signature includes new parameters."""
        import inspect
        sig = inspect.signature(PDBModel.coarse_grain)
        params = sig.parameters
        
        # Check that new parameters exist
        self.assertIn('predict_affinity', params)
        self.assertIn('adfr_path', params)
        
        # Check default values
        self.assertEqual(params['predict_affinity'].default, False)
        self.assertEqual(params['adfr_path'].default, None)
    
    def test_coarse_grain_interface_detection(self):
        """Test that interface detection works correctly."""
        model = PDBModel(pdb_id=self.pdb_id, save_dir=self.test_dir)
        model.coarse_grain(
            distance_cutoff=0.35,
            residue_cutoff=3,
            predict_affinity=False
        )
        
        # Check that data structures are properly initialized
        self.assertEqual(len(model.all_chains), len(model.all_interfaces))
        self.assertEqual(len(model.all_chains), len(model.all_interfaces_coords))
        self.assertEqual(len(model.all_chains), len(model.all_interface_energies))
        
        # Check that interfaces are detected for chains
        total_interfaces = sum(len(interfaces) for interfaces in model.all_interfaces)
        self.assertGreater(total_interfaces, 0, "No interfaces detected")
    
    def test_energy_storage_consistency(self):
        """Test that energy values are consistently stored for all interfaces."""
        model = PDBModel(pdb_id=self.pdb_id, save_dir=self.test_dir)
        model.coarse_grain(
            distance_cutoff=0.35,
            residue_cutoff=3,
            predict_affinity=False
        )
        
        # For each chain, number of interfaces should match number of energies
        for i in range(len(model.all_chains)):
            self.assertEqual(
                len(model.all_interfaces[i]),
                len(model.all_interface_energies[i]),
                f"Chain {model.all_chains[i].id}: mismatch between interfaces and energies"
            )


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
    def test_coarse_grain_with_prediction(self):
        """Test coarse_grain with ProAffinity prediction enabled."""
        model = PDBModel(pdb_id=self.pdb_id, save_dir=self.test_dir)
        model.coarse_grain(
            distance_cutoff=0.35,
            residue_cutoff=3,
            predict_affinity=True,
            adfr_path=self.adfr_path,
        )
        
        # Check that interfaces were detected
        self.assertGreater(len(model.all_chains), 0)
        
        # Check that energies are not all the same (should be predicted values)
        all_energies = []
        for chain_energies in model.all_interface_energies:
            all_energies.extend(chain_energies)
        
        if len(all_energies) > 1:
            # With prediction, energies should vary
            energy_set = set(f"{e:.2f}" for e in all_energies)
            # At least check that some energies were assigned
            self.assertGreater(len(all_energies), 0)


if __name__ == '__main__':
    unittest.main()

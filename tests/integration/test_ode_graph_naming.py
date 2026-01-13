"""
Integration test for ODE generation with graph-based naming.

Tests that the complete pipeline generates proper ODE models with topology-aware
complex names for various molecular systems.
"""

import unittest
import os
import shutil
from pathlib import Path
from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters


class TestODEGraphBasedNaming(unittest.TestCase):
    """Test ODE generation with graph-based naming."""
    
    def setUp(self):
        """Set up test workspace."""
        self.test_workspace = "test_ode_integration"
        if os.path.exists(self.test_workspace):
            shutil.rmtree(self.test_workspace)
    
    def tearDown(self):
        """Clean up test workspace."""
        if os.path.exists(self.test_workspace):
            shutil.rmtree(self.test_workspace)
    
    def test_ode_with_graph_naming_8y7s(self):
        """Test ODE generation with graph-based naming for 8y7s structure."""
        # Build system with ODE enabled
        model = PDBModelBuilder(source="8y7s")
        hyperparams = PDBModelHyperparameters(
            interface_detect_distance_cutoff=1.0,
            ode_enabled=True,
            ode_time_span=(0.0, 1.0),  # Shorter for faster test
            ode_plot=False,  # Skip plotting for test
            ode_save_csv=True
        )
        
        # Build system
        system = model.build_system(
            workspace_path=self.test_workspace,
            hyperparams=hyperparams
        )
        
        # Check ODE results exist
        ode_csv = Path(self.test_workspace) / "ode_results" / "ode_solution.csv"
        self.assertTrue(ode_csv.exists(), "ODE results CSV should be generated")
        
        # Read CSV and check for graph-based naming
        with open(ode_csv, 'r') as f:
            header = f.readline().strip()
            species_names = header.split(',')[1:]  # Skip time column
        
        # Should have multiple species
        self.assertGreater(len(species_names), 1, "Should generate multiple species")
        
        # Check for topology-aware naming (should contain topology keywords or hashes)
        has_topology_naming = any(
            '_' in name and any(topo in name.lower() for topo in ['linear', 'cyclic', 'complete', 'branched']) or
            any(c.isdigit() for c in name.split('_')[-1] if '_' in name)  # Has hash suffix
            for name in species_names if name != species_names[0]  # Skip monomer
        )
        
        self.assertTrue(
            has_topology_naming,
            f"Should use graph-based naming. Got species: {species_names}"
        )
        
        # Check monomer is just molecule type name
        self.assertEqual(species_names[0], 'A', "Monomer should be named 'A'")
    
    def test_ode_species_names_descriptive(self):
        """Test that ODE species names are descriptive and include topology."""
        model = PDBModelBuilder(source="8y7s")
        hyperparams = PDBModelHyperparameters(
            interface_detect_distance_cutoff=1.0,
            ode_enabled=True,
            ode_time_span=(0.0, 0.5),
            ode_plot=False,
            ode_save_csv=True
        )
        
        system = model.build_system(
            workspace_path=self.test_workspace,
            hyperparams=hyperparams
        )
        
        # Read species names
        ode_csv = Path(self.test_workspace) / "ode_results" / "ode_solution.csv"
        with open(ode_csv, 'r') as f:
            header = f.readline().strip()
            species_names = header.split(',')[1:]
        
        # Check that complex names contain molecule count
        for name in species_names[1:]:  # Skip monomer
            # Should have format like "A2_complete_7937" or similar
            parts = name.split('_')
            self.assertGreaterEqual(len(parts), 2, 
                f"Complex name '{name}' should have at least 2 parts (composition and topology/hash)")
            
            # First part should indicate molecule type and count (e.g., "A2", "A3")
            self.assertTrue(parts[0].startswith('A'), 
                f"Complex name '{name}' should start with molecule type")


if __name__ == '__main__':
    unittest.main()

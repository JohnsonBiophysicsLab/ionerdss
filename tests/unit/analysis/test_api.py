"""
Unit tests for ionerdss.analysis.api

Tests the Analyzer API for loading and computing on simulations.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import pandas as pd
from ionerdss.analysis import Analyzer


class TestAnalyzerAPI(unittest.TestCase):
    """Test cases for Analyzer API."""

    def setUp(self):
        """Set up mock simulation directory for testing."""
        # Create temporary directory structure
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Create simulation directory structure
        # /temp_dir/
        #   /1/  (simulation 1)
        #     /DATA/
        #       histogram_complexes_time.dat
        #       copy_numbers_time.dat
        #       transition_matrix_time.dat
        
        sim_dir = self.temp_path / "1"
        sim_dir.mkdir()
        data_dir = sim_dir / "DATA"
        data_dir.mkdir()
        
        # Create histogram file
        histogram_file = data_dir / "histogram_complexes_time.dat"
        with open(histogram_file, 'w') as f:
            f.write("time: 0.0\n")
            f.write("10\tA: 1.\n")
            f.write("\n")
        
        # Create copy numbers file
        copy_numbers_file = data_dir / "copy_numbers_time.dat"
        with open(copy_numbers_file, 'w') as f:
            f.write("Time,Complex,A\n")
            f.write("0.0,10,10\n")
        
        # Create transition matrix file
        transition_file = data_dir / "transition_matrix_time.dat"
        with open(transition_file, 'w') as f:
            f.write("time: 0.0\n")
            f.write("transition matrix for each mol type:\n")
            f.write("5 0\n")
            f.write("0 5\n")
            f.write("\n")
            f.write("time: 0.1\n")
            f.write("transition matrix for each mol type:\n")
            f.write("3 2\n")
            f.write("1 4\n")
            f.write("\n")
        
        self.mock_simulation_dir = self.temp_path
    
    def tearDown(self):
        """Clean up temporary files."""
        if Path(self.temp_dir).exists():
            shutil.rmtree(self.temp_dir)

    def test_analyzer_loading(self):
        """Test Analyzer initialization and simulation loading."""
        analyzer = Analyzer(self.mock_simulation_dir)
        
        self.assertEqual(len(analyzer.simulations), 1)
        sim = analyzer.get_simulation(0)
        
        self.assertEqual(sim.id, "1")
        
        # Test lazy loading
        self.assertIsNone(sim._data)
        sim.load()
        self.assertIsNotNone(sim._data)
        self.assertEqual(len(sim.data.transitions), 2)
        self.assertIsNotNone(sim.data.copy_numbers)

    def test_analyzer_integration_compute(self):
        """Test Analyzer compute methods integration."""
        analyzer = Analyzer(self.mock_simulation_dir)
        sim = analyzer.get_simulation(0)
        
        # Compute Free Energy (should trigger load)
        df_fe = analyzer.compute_free_energy(sim)
        
        self.assertFalse(df_fe.empty)
        self.assertIn('free_energy', df_fe.columns)
        
        # Check caching
        self.assertIsNotNone(sim.data.df_free_energy)
        self.assertIs(sim.data.df_free_energy, df_fe)


if __name__ == '__main__':
    unittest.main(verbosity=2)

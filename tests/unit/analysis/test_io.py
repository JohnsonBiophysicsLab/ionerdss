"""
Unit tests for ionerdss.analysis.io

Tests file parsing functions for NERDSS simulation output files.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from ionerdss.analysis.io import parser


class TestIOParser(unittest.TestCase):
    """Test cases for IO parser functions."""

    def setUp(self):
        """Set up temporary test files for each test."""
        # Create a temporary directory
        self.temp_dir = tempfile.mkdtemp()
        self.temp_path = Path(self.temp_dir)
        
        # Create sample transition file
        self.sample_transition_file = self.temp_path / "transition.dat"
        with open(self.sample_transition_file, 'w') as f:
            f.write("time: 0.0\n")
            f.write("transition matrix for each mol type:\n")
            f.write("0 0\n")
            f.write("0 0\n")
            f.write("\n")
            f.write("time: 0.1\n")
            f.write("transition matrix for each mol type:\n")
            f.write("0 2\n")
            f.write("1 0\n")
            f.write("lifetime for each mol type:\n")
            f.write("size of the cluster: 2\n")
            f.write("0.5 0.5\n")
            f.write("\n")
        
        # Create sample copy numbers file (CSV format)
        self.sample_copy_numbers_file = self.temp_path / "copy_numbers.dat"
        with open(self.sample_copy_numbers_file, 'w') as f:
            f.write("Time,Complex,A,B\n")
            f.write("0.0,10,10,10\n")
            f.write("0.1,9,8,9\n")
            f.write("0.2,8,7,8\n")
        
        # Create sample complex histogram file
        self.sample_complex_histogram_file = self.temp_path / "complex_histogram.dat"
        with open(self.sample_complex_histogram_file, 'w') as f:
            f.write("time: 0.0\n")
            f.write("10\tA: 1.\n")
            f.write("\n")
            f.write("time: 0.1\n")
            f.write("5\tA: 1. B: 1.\n")
            f.write("3\tA: 2.\n")
            f.write("\n")
    
    def tearDown(self):
        """Clean up temporary files."""
        if Path(self.temp_dir).exists():
            shutil.rmtree(self.temp_dir)

    @unittest.skip("Skipping transition matrix tests as requested")
    def test_parse_transition_file(self):
        """Test transition file parsing."""
        pass
        
    def test_parse_copy_numbers(self):
        """Test copy numbers file parsing."""
        df = parser.parse_copy_numbers(self.sample_copy_numbers_file)
        self.assertFalse(df.empty)
        self.assertEqual(len(df), 3)
        self.assertIn('Complex', df.columns)
        self.assertEqual(df.iloc[1]['A'], 8)

    def test_parse_complex_histogram(self):
        """Test complex histogram file parsing with real data."""
        # Use the provided test file in tests/data
        # Path is relative to this test file: ../../../data/test_6BNO_histogram_complexes_time.dat
        test_file = Path(__file__).resolve().parent.parent.parent.parent / "tests/data/test_6BNO_histogram_complexes_time.dat"
        
        if not test_file.exists():
            self.skipTest(f"Test data file not found: {test_file}")
            
        times, all_comps, hist_matrix = parser.parse_complex_histogram(test_file)
        
        # Based on file inspection:
        # Time 0: 1 line
        # Time 1.35727e-05: 3 lines
        # Time 2.71454e-05: 4 lines
        # Time 4.07181e-05: 4 lines
        # Total = 4 time points
        self.assertEqual(len(times), 4)
        
        # Test first time point (Time 0)
        # Content: "75 A: 1."
        # This matches index 0
        self.assertEqual(times[0], 0.0)
        # Find the row in matrix corresponding to t=0
        # Use _getrow() for sparse matrix to avoid NotImplementedError on 1D slicing
        # Do we still need this since we bumped up to python 3.10? - M. Ying
        row_0 = hist_matrix._getrow(0).toarray().flatten()
        # It should have exactly one non-zero entry (75)
        self.assertEqual(row_0.sum(), 75)
        
        # Test last time point matches counts
        # Content: 3 (A:3), 28 (A:1), 5 (A:4), 9 (A:2) -> Total 3+28+5+9 = 45
        row_last = hist_matrix._getrow(len(times)-1).toarray().flatten()
        self.assertEqual(row_last.sum(), 45)


if __name__ == '__main__':
    unittest.main(verbosity=2)

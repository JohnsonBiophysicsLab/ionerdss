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

    def test_parse_transition_file(self):
        """Test transition file parsing."""
        transitions, lifetimes = parser.parse_transition_file(self.sample_transition_file)
        
        # Check transitions
        self.assertEqual(len(transitions), 2)
        self.assertEqual(transitions[0]['time'], 0.0)
        self.assertTrue(np.all(transitions[0]['matrix'] == np.zeros((2, 2))))
        
        self.assertEqual(transitions[1]['time'], 0.1)
        expected_mat = np.array([[0, 2], [1, 0]])
        self.assertTrue(np.all(transitions[1]['matrix'] == expected_mat))
        
        # Check lifetimes
        self.assertEqual(len(lifetimes), 1)  # Only one time point had lifetimes
        self.assertEqual(lifetimes[0]['time'], 0.1)
        self.assertEqual(lifetimes[0]['lifetimes'][2], [0.5, 0.5])

    def test_parse_copy_numbers(self):
        """Test copy numbers file parsing."""
        df = parser.parse_copy_numbers(self.sample_copy_numbers_file)
        self.assertFalse(df.empty)
        self.assertEqual(len(df), 3)
        self.assertIn('Complex', df.columns)
        self.assertEqual(df.iloc[1]['A'], 8)

    def test_parse_complex_histogram(self):
        """Test complex histogram file parsing."""
        data = parser.parse_complex_histogram(self.sample_complex_histogram_file)
        self.assertEqual(len(data), 2)
        
        t0 = data[0]
        self.assertEqual(t0['time'], 0.0)
        self.assertEqual(len(t0['complexes']), 1)
        self.assertEqual(t0['complexes'][0]['count'], 10)
        self.assertEqual(t0['complexes'][0]['composition'], {'A': 1})
        
        t1 = data[1]
        self.assertEqual(t1['time'], 0.1)
        self.assertEqual(len(t1['complexes']), 2)
        # Check "5 A: 1. B: 1."
        c1 = t1['complexes'][0]
        self.assertEqual(c1['count'], 5)
        self.assertEqual(c1['composition'], {'A': 1, 'B': 1})


if __name__ == '__main__':
    unittest.main(verbosity=2)

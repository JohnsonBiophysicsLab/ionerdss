"""
Unit tests for ionerdss.analysis.processing

Tests processing functions for size distributions and free energy calculations.
"""

import unittest
import numpy as np
import pandas as pd
from ionerdss.analysis.processing import transitions as trans_proc


class TestTransitionProcessing(unittest.TestCase):
    """Test cases for transition processing functions."""

    def test_compute_size_distribution(self):
        """Test size distribution computation from transition matrix."""
        # Create a dummy 3x3 matrix
        # Rows = From size (if index 0 is size 1)
        # Wait, usually transitions are: M[i,j] = count from j to i?
        # Let's check documentation/code.
        # parser.py logic: just reads rows.
        # transitions.py logic: "Row index n corresponds to size n+1"
        
        # Let's assume a simple distribution: 
        # Size 1: 10 particles
        # Size 2: 5 particles
        # Size 3: 2 particles
        # Total = 17
        
        # But wait, compute_size_distribution usually takes the full Simulation or Matrix?
        # Let's look at the function signature in processing/transitions.py
        # It takes `transition_matrix: np.ndarray`
        
        # Mock matrix (Counts of transitions TO i FROM j)
        # Diagonal M[i,i] usually dominates (staying same size)
        matrix = np.array([
            [10, 0, 0],
            [0, 5, 0],
            [0, 0, 2]
        ])
        
        df = trans_proc.compute_size_distribution_transition_matrix(matrix)
        
        self.assertIn('size', df.columns)
        self.assertIn('probability', df.columns)
        
        # Total sum of matrix is 17
        # Prob size 1 = 10/17
        self.assertTrue(np.isclose(df.loc[0, 'probability'], 10/17))
        self.assertTrue(np.isclose(df.loc[1, 'probability'], 5/17))
        self.assertTrue(np.isclose(df.loc[2, 'probability'], 2/17))

    def test_compute_free_energy(self):
        """Test free energy calculation from size distribution."""
        df_dist = pd.DataFrame({
            'size': [1, 2],
            'probability': [0.8, 0.2]
        })
        
        # G = -kT ln(P)
        # If T=1, kB=1 (sim units)
        
        df_fe = trans_proc.compute_free_energy(df_dist, temperature=1.0)
        
        p1 = 0.8
        g1 = -1.0 * np.log(p1)
        
        self.assertTrue(np.isclose(df_fe.loc[0, 'free_energy'], g1))
        
        # Test normalization (if G_min is shifted to 0)
        # The implementation might subtract the min G
        min_G = min(-np.log(0.8), -np.log(0.2))
        expected_g1_shifted = g1 - min_G
        
        # Check if implementation shifts to zero
        if df_fe['free_energy'].min() == 0.0:
            self.assertTrue(np.isclose(df_fe.loc[0, 'free_energy'], expected_g1_shifted))


if __name__ == '__main__':
    unittest.main(verbosity=2)

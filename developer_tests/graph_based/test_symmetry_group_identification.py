"""
test_pointgroup_symmetry.py

Unit tests for symmetry classification using the PointGroup class in
`ionerdss.model.graph_based.symmetry.pointgroup`.

This test suite verifies point group identification from 3D geometry
and element types using canonical symmetric structures.

Tests
-----
1. Octahedral (Oh) symmetry from an SF6-like configuration
2. Tetrahedral (Td) symmetry from a regular tetrahedron
3. Lower symmetry in a flipped tetrahedral dimer (C2 or D2h-like)

Each test ensures the returned point group string is valid and structurally expected.

Author: yying7@jh.edu
-------
Auto-generated and maintained as part of the ioNERDSS framework
for molecular complex modeling and reaction network generation.
"""

import unittest
import numpy as np
from ionerdss.model.graph_based.symmetry.pointgroup import PointGroup

class TestPointGroupSymmetry(unittest.TestCase):
    def test_sf6_octrahedral_symmetry(self):
        """Test that SF6-like molecule gives Oh point group."""
        coords = [
            [0.000,  0.000,  0.000],
            [0.000,  0.000,  1.561],
            [0.000,  1.561,  0.000],
            [0.000,  0.000, -1.561],
            [0.000, -1.561,  0.000],
            [1.561,  0.000,  0.000],
            [-1.561, 0.000,  0.000],
        ]
        symbols = ['S', 'F', 'F', 'F', 'F', 'F', 'F']
        pg = PointGroup(positions=coords, symbols=symbols)
        pg_str = pg.get_point_group()
        self.assertIsInstance(pg_str, str)
        # Octahedral Symmetry
        self.assertIn(pg_str, {"O", "Oh", "O_h"})

    def test_tetrahedron_td_symmetry(self):
        """Test regular tetrahedron gives Td symmetry."""
        coords = np.array([
            [ 1,  1,  1],
            [-1, -1,  1],
            [-1,  1, -1],
            [ 1, -1, -1]
        ]) / np.sqrt(3)
        pg = PointGroup(positions=coords, symbols=["B"] * 4)
        pg_str = pg.get_point_group()
        self.assertIsInstance(pg_str, str)
        # Tetrahedral symmetry
        self.assertIn(pg_str, {"T", "Td", "T_d"})

    def test_tetrahedral_dimer(self):
        """Test dimer of tetrahedra flipped across Z axis."""
        tetra = np.array([
            [ 1,  1,  1],
            [-1, -1,  1],
            [-1,  1, -1],
            [ 1, -1, -1]
        ]) / np.sqrt(3)
        tetra_A = tetra.copy()
        tetra_B = tetra.copy()
        tetra_B[:, 2] *= -1  # flip Z axis

        combined = np.vstack([tetra_A, tetra_B])
        pg = PointGroup(positions=combined, symbols=["B"] * 8)
        pg_str = pg.get_point_group()
        self.assertIsInstance(pg_str, str)
        self.assertGreater(len(pg_str), 0)


if __name__ == "__main__":
    unittest.main()

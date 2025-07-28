"""
test_geometry.py

Unit tests for geometry utilities in NERDSS, particularly rigid chain alignment
and steric clash detection.

This test suite ensures:
- Rigid alignment between chains works with at least 3 residues.
- Steric clash detection behaves as expected for overlapping and non-overlapping spheres.
"""

import unittest
import numpy as np
from Bio.PDB import Chain, Residue, Atom

from ionerdss.model.pdb.geometry import rigid_transform_chains, check_steric_clashes

class TestGeometryUtilities(unittest.TestCase):

    def setUp(self):
        # Create two chains with 3 CA atoms each (requirement)
        self.chain1 = Chain.Chain('A')
        self.chain2 = Chain.Chain('B')
        atoms1 = [(-1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
        atoms2 = [(1.0, 1.0, 0.0), (2.0, 1.0, 0.0), (3.0, 1.0, 0.0), (4.0, 1.0, 0.0)]


        for i, (a1, a2) in enumerate(zip(atoms1, atoms2)):
            res_id = (' ', i + 1, ' ')
            resname = 'ALA'

            res1 = Residue.Residue(res_id, resname, '')
            res1.add(Atom.Atom('CA', np.array(a1), 1.0, 1.0, '', 'CA', i))
            self.chain1.add(res1)

            res2 = Residue.Residue(res_id, resname, '')
            res2.add(Atom.Atom('CA', np.array(a2), 1.0, 1.0, '', 'CA', i))
            self.chain2.add(res2)

    def test_rigid_transform_chains(self):
        """Tests whether rigid_transform_chains returns a valid rotation and translation."""
        R, t = rigid_transform_chains(self.chain1, self.chain2)

        # Check orthogonality: R.T @ R ≈ I
        identity = np.dot(R.T, R)
        np.testing.assert_allclose(identity, np.eye(3), atol=1e-6)

        # Check determinant is 1 (proper rotation)
        self.assertAlmostEqual(np.linalg.det(R), 1.0, places=6)

        # Apply transform to all points from chain1 and verify alignment with chain2
        ca1 = np.array([res['CA'].coord for res in self.chain1])
        ca2 = np.array([res['CA'].coord for res in self.chain2])
        transformed = (R @ ca1.T).T + t
        np.testing.assert_allclose(transformed, ca2, atol=1e-6)

    def test_check_steric_clashes_clashing(self):
        """Test that steric clash is detected correctly."""
        pos1 = np.array([0.0, 0.0, 0.0])
        pos2 = np.array([1.0, 0.0, 0.0])
        r1 = 1.0
        r2 = 1.0
        buffer = 0.0

        self.assertTrue(check_steric_clashes(pos1, pos2, r1, r2, buffer))

    def test_check_steric_clashes_no_clash(self):
        """Test that no clash is detected when spheres are sufficiently far."""
        pos1 = np.array([0.0, 0.0, 0.0])
        pos2 = np.array([3.0, 0.0, 0.0])
        r1 = 1.0
        r2 = 1.0
        buffer = 0.0

        self.assertFalse(check_steric_clashes(pos1, pos2, r1, r2, buffer))


if __name__ == "__main__":
    unittest.main()

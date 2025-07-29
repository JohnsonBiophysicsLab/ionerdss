"""
Unit tests for the coarse-graining and interface detection functions
in `ionerdss.model.pdb.coarse_grain`.

This test suite validates:
- The functionality of `compute_interface` for detecting residue-residue interactions
  between two chains based on spatial proximity and an energy scoring table.
- The correctness of `coarse_grain_structure` in identifying chains, computing
  center-of-mass (COM), estimating radius, and finding valid binding interfaces
  between chains in a Biopython `Structure`.

Test cases include:
- Synthetic PDB-like structure with two interacting chains containing alpha-carbon (CA) atoms.
- Verification of returned COMs, interface partners, interacting residue lists,
  and calculated interaction energies.
- Use of a default energy table from `ionerdss.model.pdb.energy_table`.

The setup avoids I/O and constructs test structures entirely in memory using Biopython objects.
"""
import unittest
import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

from ionerdss.model.pdb.coarse_grain import coarse_grain_structure, compute_interface
from ionerdss.model.pdb.energy_table import get_default_energy_table
from ionerdss.utils.coords import Coords

class TestCoarseGrainStructure(unittest.TestCase):
    def setUp(self):
        # Create synthetic structure with 2 chains, each with 3 residues having CA atoms
        self.structure = Structure("test")
        model = Model(0)
        self.structure.add(model)

        chainA = Chain("A")
        chainB = Chain("B")
        model.add(chainA)
        model.add(chainB)

        def add_residue(chain, res_id, res_name, coord):
            res = Residue((" ", res_id, " "), res_name, "")
            atom = Atom("CA", np.array(coord, dtype=float), 1.0, 1.0, "", "CA", 0, "C")
            res.add(atom)
            chain.add(res)

        # Chain A: clustered around origin
        add_residue(chainA, 1, "ALA", [0.0, 0.0, 0.0])
        add_residue(chainA, 2, "GLY", [0.3, 0.0, 0.0])
        add_residue(chainA, 3, "SER", [0.1, 0.3, 0.0])

        # Chain B: close enough to interact
        add_residue(chainB, 1, "VAL", [0.25, 0.05, 0.2])
        add_residue(chainB, 2, "LEU", [0.6, 0.0, 0.0])
        add_residue(chainB, 3, "THR", [0.5, 0.3, 0.1])

        self.distance_cutoff = 0.35
        self.residue_cutoff = 2
        self.energy_table = get_default_energy_table()

    def test_compute_interface_returns_valid_result(self):
        chainA = self.structure[0]["A"]
        chainB = self.structure[0]["B"]

        residues_i = [(res.id[1], res.get_resname().upper(), res['CA'].coord)
                      for res in chainA]
        residues_j = [(res.id[1], res.get_resname().upper(), res['CA'].coord)
                      for res in chainB]

        result = compute_interface(residues_i, residues_j, self.energy_table,
                                   self.distance_cutoff, self.residue_cutoff)
        self.assertIsNotNone(result, "Interface should be detected with close chains")

        com_i, com_j, ids_i, ids_j, energy = result
        self.assertIsInstance(com_i, Coords)
        self.assertIsInstance(com_j, Coords)
        self.assertGreaterEqual(len(ids_i), self.residue_cutoff)
        self.assertGreaterEqual(len(ids_j), self.residue_cutoff)
        self.assertIsInstance(energy, float)

    def test_coarse_grain_structure_detects_interface(self):
        result = coarse_grain_structure(self.structure, self.distance_cutoff, self.residue_cutoff)

        self.assertEqual(len(result['chains']), 2)
        self.assertEqual(len(result['COMs']), 2)
        self.assertTrue(all(isinstance(c, Coords) or c is None for c in result['COMs']))

        # Should detect 1 interface per chain
        self.assertEqual(len(result['interfaces'][0]), 1)
        self.assertEqual(result['interfaces'][0][0], "B")
        self.assertEqual(result['interfaces'][1][0], "A")

        # Check energy and residue information
        self.assertAlmostEqual(result['interface_energies'][0][0], -1.2)
        self.assertIsInstance(result['interface_coords'][0][0], Coords)
        self.assertIsInstance(result['interface_residues'][0][0], list)

if __name__ == "__main__":
    unittest.main()

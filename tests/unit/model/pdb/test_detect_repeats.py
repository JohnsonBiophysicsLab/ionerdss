"""
test_detect_repeats.py

Unit tests for the repeated chain identification logic in `identify_repeated_chains`.

This module validates:
- Header-based chain grouping (from mmCIF metadata)
- Fallback to structure-based chain grouping using C-alpha superposition
- Fallback to sequence-based grouping using pairwise sequence identity
- Correct parsing and assignment of canonical representative chains
- Mode-switching behavior via the `mode` argument ('default', 'structure', 'sequence')

Synthetic structures are created in-memory to simulate controlled chain duplications
and are compatible with both structural and sequence similarity testing.

Dependencies:
- Biopython
- Python unittest
"""

import unittest
import warnings
import numpy as np
from Bio.PDB import Structure, Model, Chain, Residue, Atom

from ionerdss.model.pdb.detect_repeats import (
    identify_repeated_chains,
)
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

class TestIdentifyRepeatedChains(unittest.TestCase):
    
    def setUp(self):
        self.params = PDBModelHyperparameters()

    def build_chain(self, chain_id, coords, seq_resnames=None):
        """Create a chain with N-CA-C atoms so that PPBuilder can build peptides."""
        chain = Chain.Chain(chain_id)
        for i, coord in enumerate(coords):
            coord = np.array(coord)  # make sure it's an array
            resname = seq_resnames[i] if seq_resnames else "ALA"
            res = Residue.Residue((" ", i + 1, " "), resname, "")
            ca = Atom.Atom("CA", coord, 1.0, 1.0, "", "CA", i, "C")
            n = Atom.Atom("N", coord - np.array([1.2, 0, 0]), 1.0, 1.0, "", "N", i, "N")
            c = Atom.Atom("C", coord + np.array([1.2, 0, 0]), 1.0, 1.0, "", "C", i, "C")
            res.add(n)
            res.add(ca)
            res.add(c)
            chain.add(res)
        return chain

    def build_structure(self, chain_specs):
        """
        chain_specs: list of tuples (chain_id, coords, optional seq)
        """
        structure = Structure.Structure("test")
        model = Model.Model(0)
        structure.add(model)
        for spec in chain_specs:
            chain_id, coords = spec[:2]
            seq = spec[2] if len(spec) > 2 else None
            chain = self.build_chain(chain_id, coords, seq)
            model.add(chain)
        return structure

    def test_structure_mode_identifies_redundant_chains(self):
        """
        The structures provided are translated so they should fully
        overlap with each other and only has 1 chain group based on
        structure mode
        coords[0]: baseline
        coords[1]: exact copy
        coords[2]: translation
        coords[3]: slight deviation
        """
        coords = [[(0,0,0), (1,0,0), (7,9,11), (10,20,30)],
                  [(0,0,0), (1,0,0), (7,9,11), (10,20,30)],
                  [(5,5,5), (6,5,5), (12,14,16), (15,25,35)],
                  [(0,0.01,0.01), (0.99,0,0), (7.98,9.02,11), (10.01,20.02,30.03)]]
        structure = self.build_structure([('A', coords[0]), ('B', coords[1]), ('C', coords[2]), ('D', coords[3])])
        chains_map, chains_group = identify_repeated_chains("test", structure, matching_mode='structure', params=self.params)
        self.assertEqual(len(chains_group), 1)
        self.assertEqual(chains_map['A'], chains_map['B'])
        self.assertEqual(chains_map['A'], chains_map['C'])
        self.assertEqual(chains_map['A'], chains_map['D'])

    def test_sequence_mode_identifies_redundant_chains(self):
        coords = [[(0,0,0), (1,0,0), (2,0,0)]]*3
        seq1 = ["ALA", "GLY", "VAL"]
        seq2 = ["ALA", "GLY", "VAL"]
        seq3 = ["ARG", "ASP", "GLU"]
        structure = self.build_structure([
            ('A', coords[0], seq1),
            ('B', coords[1], seq2),
            ('C', coords[2], seq3),
        ])
        chains_map, chains_group = identify_repeated_chains("test", structure, matching_mode='sequence', params=self.params)
        print(chains_map)
        print(chains_group)
        self.assertEqual(len(chains_group), 2)
        self.assertEqual(chains_map['A'], chains_map['B'])
        self.assertNotEqual(chains_map['A'], chains_map['C'])

    def test_default_mode_uses_mmcif_header(self):
        structure = self.build_structure([('A', [(0,0,0)]), ('B', [(1,0,0)]), ('C', [(2,0,0)])])
        structure.header = {
            'mmcif_dict': {
                '_entity_poly.entity_id': ['1', '2'],
                '_entity_poly.pdbx_strand_id': ['A,B', 'C']
            }
        }
        chains_map, chains_group = identify_repeated_chains("test", structure, matching_mode='default', params=self.params)
        self.assertIn(['A', 'B'], chains_group)
        self.assertIn(['C'], chains_group)
        self.assertEqual(chains_map['A'], chains_map['B'])

    def test_default_mode_warns_and_falls_back(self):
        coords = [[(0,0,0), (1,0,0), (7,9,11), (10,20,30)],
                  [(0,0,0), (1,0,0), (7,9,11), (10,20,30)]]
        seq = ["ALA", "GLY", "VAL", "ALA"]
        structure = self.build_structure([
            ('X', coords[0], seq),
            ('Y', coords[1], seq)
        ])
        structure.header = {}  # no mmCIF info
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            chains_map, chains_group = identify_repeated_chains("test", structure, matching_mode='default', params=self.params)
            self.assertTrue(any("Falling back" in str(wi.message) for wi in w))
        self.assertEqual(len(chains_group), 1)
        self.assertEqual(chains_map['X'], chains_map['Y'])

    def test_invalid_mode_raises_error(self):
        invalid_mode_params = PDBModelHyperparameters(options={"matching_mode" : "invalid_mode"})
        structure = self.build_structure([('A', [(0,0,0)])])
        with self.assertRaises(ValueError):
            identify_repeated_chains("test", structure,
                                     matching_mode='invalid_mode',
                                     params=invalid_mode_params)

if __name__ == "__main__":
    unittest.main()

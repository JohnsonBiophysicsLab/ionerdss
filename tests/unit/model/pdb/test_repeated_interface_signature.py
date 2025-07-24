"""
Unit tests for repeated_interface_signature.py

This module verifies correctness of geometric and sequence-based signature identification
and comparison utilities used in symmetric protein interface analysis.
"""

import unittest
import numpy as np

from ionerdss.math.coords import Coords
from ionerdss.math.angles import angles_from_points
from ionerdss.model.pdb.repeated_interface_signature import (
    identify_interface_signature,
    identify_interface_sequence_signature,
    signature_are_similar,
    signature_difference,
    build_signature,
    invert_signature
)


class DummyResidue:
    def __init__(self, name):
        self.name = name

    def get_resname(self):
        return self.name


class TestRepeatedInterfaceSignature(unittest.TestCase):

    def setUp(self):
        self.cg_model = {
            "chains": [DummyChain("A"), DummyChain("B")],
            "COMs": [Coords(0, 0, 0), Coords(2, 0, 0)],
            "interface_coords": [[Coords(1, 0, 0)], [Coords(1, 0, 0)]],
            "interfaces": [["B"], ["A"]],
            "interface_residues": [
                [[DummyResidue("GLY"), DummyResidue("ALA")]],
                [[DummyResidue("ALA"), DummyResidue("GLY")]]
            ]
        }

    def test_identify_interface_signature(self):
        sig, sig_inv, I_B, k = identify_interface_signature("A", "B", self.cg_model, 0)
        self.assertAlmostEqual(sig["dA"], 1.0)
        self.assertAlmostEqual(sig["dB"], 1.0)
        self.assertAlmostEqual(sig["dAB"], 0.0)
        self.assertEqual(k, 0)

    def test_sequence_signature(self):
        sig, sig_inv, _, k = identify_interface_sequence_signature("A", "B", 0, self.cg_model)
        self.assertEqual(sig["seqA"], "ALAGLY")
        self.assertEqual(sig["seqB"], "ALAGLY")
        self.assertEqual(sig_inv["seqA"], sig["seqB"])

    def test_signature_similarity(self):
        sig1 = build_signature(Coords(0, 0, 0), Coords(1, 0, 0), Coords(2, 0, 0), Coords(1, 0, 0))
        sig2 = build_signature(Coords(0.01, 0, 0), Coords(1.01, 0, 0), Coords(1.99, 0, 0), Coords(1.01, 0, 0))
        similar = signature_are_similar(sig1, sig2, 0.1, 0.1, 5.0)
        self.assertTrue(similar)

    def test_signature_difference(self):
        sig1 = build_signature(Coords(0, 0, 0), Coords(1, 0, 0), Coords(2, 0, 0), Coords(1, 0, 0))
        sig2 = invert_signature(sig1)
        diff = signature_difference(sig1, sig2)
        self.assertEqual(diff, 0.0)


class DummyChain:
    def __init__(self, id):
        self.id = id


if __name__ == '__main__':
    unittest.main()

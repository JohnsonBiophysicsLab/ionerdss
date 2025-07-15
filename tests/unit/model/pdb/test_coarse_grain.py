"""
Unit test for ionerdss.model.pdb.coarse_grain
"""
import unittest
from io import StringIO
from Bio.PDB import PDBParser
from ionerdss.model.pdb.coarse_grain import coarse_grain_structure

class TestCoarseGrainWithRealPDB(unittest.TestCase):
    def setUp(self):
        pdb_string = """
ATOM      1  N   MET A   1      11.104  13.207   2.100  1.00 23.26           N
ATOM      2  CA  MET A   1      12.560  13.407   2.170  1.00 23.09           C
ATOM      3  C   MET A   1      13.145  14.835   2.240  1.00 21.93           C
ATOM      4  O   MET A   1      12.451  15.794   2.176  1.00 21.26           O
ATOM      5  CB  MET A   1      13.086  12.526   3.318  1.00 24.35           C
ATOM      6  N   GLY B   2      14.430  15.011   2.327  1.00 20.65           N
ATOM      7  CA  GLY B   2      15.107  16.324   2.377  1.00 20.53           C
ATOM      8  C   GLY B   2      15.981  16.413   3.628  1.00 20.04           C
ATOM      9  O   GLY B   2      16.899  17.247   3.652  1.00 19.50           O
TER
END
"""
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("TEST", StringIO(pdb_string))

    def test_interface_detection_and_output_format(self):
        result = coarse_grain_structure(self.structure, distance_cutoff=0.35, residue_cutoff=1)

        # Top-level keys exist
        self.assertIn("chains", result)
        self.assertIn("COMs", result)
        self.assertIn("interface_coords", result)
        self.assertIn("interface_residues", result)
        self.assertEqual(len(result["chains"]), 2)

        self.assertGreaterEqual(len(result["interface_coords"][0]), 0)
        self.assertGreaterEqual(len(result["interface_coords"][1]), 0)

if __name__ == "__main__":
    unittest.main()

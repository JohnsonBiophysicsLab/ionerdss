"""
tests/test_hyperparameters.py

Unit tests for `PDBModelHyperparameters`.

This test suite verifies that:
  1) Default values are correctly assigned when no options are provided.
  2) User-supplied options override the defaults as expected.

The tests are self-contained and do not require external data files.
They serve as a guard against regressions when adding new hyperparameters
or changing default values in `model/pdb/hyperparameters.py`.

How to run:
    - With unittest:
        python -m unittest -v tests/test_hyperparameters.py
    - With pytest (also works):
        pytest -q tests/test_hyperparameters.py
"""

import unittest
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

class TestPDBModelHyperparameters(unittest.TestCase):

    def test_defaults(self):
        """Test that default hyperparameters are correctly assigned."""
        hp = PDBModelHyperparameters()
        self.assertEqual(hp.distance_cutoff, 0.6)
        self.assertEqual(hp.residue_cutoff, 3)
        # skip energy table, we are not using it any more
        # the code is there for backwards compatibility
        self.assertEqual(hp.rmsd_threshold, 2.0)
        self.assertEqual(hp.seq_threshold, 0.5)
        self.assertIsNone(hp.custom_aligner)
        self.assertEqual(hp.matching_mode, "default")
        self.assertEqual(hp.dist_threshold_intra, 3.5)
        self.assertEqual(hp.dist_threshold_inter, 3.5)
        self.assertEqual(hp.angle_threshold, 25.0)
        self.assertFalse(hp.standard_output)
        self.assertFalse(hp.verbose_mode)

    def test_overrides(self):
        """Test that options dictionary overrides defaults."""
        options = {
            "distance_cutoff": 1.2,
            "residue_cutoff": 10,
            "rmsd_threshold": 3.0,
            "seq_threshold": 0.9,
            "custom_aligner": "aligner",
            "matching_mode": "sequence",
            "dist_threshold_intra": 2.0,
            "dist_threshold_inter": 4.0,
            "angle_threshold": 30.0,
            "standard_output": True,
            "verbose_mode": True
        }
        hp = PDBModelHyperparameters(options)
        for key, value in options.items():
            self.assertEqual(getattr(hp, key), value)

if __name__ == "__main__":
    unittest.main()

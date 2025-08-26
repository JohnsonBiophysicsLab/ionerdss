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
    
    def setUp(self):
        self.options = {
            "distance_cutoff": 1.1,
            "residue_cutoff": 7,
            "energy_table": "etable",
            "rmsd_threshold": 3.3,
            "seq_threshold": 0.8,
            "custom_aligner": "aligner_fn",
            "matching_mode": "sequence",
            "dist_threshold_intra": 2.5,
            "dist_threshold_inter": 4.5,
            "angle_threshold": 15.0,
            "standard_output": True,
            "logger_level": "INFO",
        }
        self.hp = PDBModelHyperparameters(self.options)

    def test_repr_contains_all_fields(self):
        """__repr__ should include all hyperparameters in a single-line string."""
        rep = repr(self.hp)
        # Must start with class name
        self.assertTrue(rep.startswith("PDBModelHyperparameters("))
        # Each field=value pair should appear
        for key, value in self.options.items():
            self.assertIn(f"{key}=", rep)

    def test_repr_roundtrip_safe(self):
        """__repr__ should be evaluable if class is in scope (sanity check)."""
        rep = repr(self.hp)
        # Don't actually eval here (would need globals), but check it's a valid Python-ish form
        self.assertIn("PDBModelHyperparameters(", rep)
        self.assertIn("distance_cutoff=", rep)

    def test_str_human_readable(self):
        """__str__ should produce multi-line output with aligned fields."""
        s = str(self.hp)
        # Header line
        self.assertIn("PDBModelHyperparameters:", s)
        # Check a few specific lines
        self.assertIn("distance_cutoff      = 1.1", s)
        self.assertIn("residue_cutoff       = 7", s)
        self.assertIn("matching_mode        = sequence", s)
        self.assertIn("logger_level         = INFO", s)

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
        self.assertEqual(hp.logger_level, "INFO")

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
            "logger_level": "INFO"
        }
        hp = PDBModelHyperparameters(options)
        for key, value in options.items():
            self.assertEqual(getattr(hp, key), value)

if __name__ == "__main__":
    unittest.main()

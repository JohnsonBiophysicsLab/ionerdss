"""
Unit tests for ionerdss.model.pdb.api

Every hyperparameter field must appear both in the reference that
set_hyperparameters.__doc__ is generated from and in print_hyperparameters, so
neither can drift away from PDBModelHyperparameters again.
"""

import contextlib
import io
import re
import unittest
from dataclasses import fields
from unittest.mock import patch

from ionerdss.model.pdb import api
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.pdb.main import PDBModelBuilder

FIELD_NAMES = [f.name for f in fields(PDBModelHyperparameters)]
# units is internal; print_hyperparameters also leaves out the aligner object
DOCUMENTED = [name for name in FIELD_NAMES if name != "units"]
PRINTED = [name for name in DOCUMENTED if name != "chain_grouping_custom_aligner"]


class TestHyperparameterReference(unittest.TestCase):
    """The hyperparameter reference covers exactly the dataclass fields."""

    def test_sections_list_each_real_field_once(self):
        listed = [name for names in api._HYPERPARAMETER_SECTIONS.values() for name in names]
        self.assertEqual(sorted(set(listed) - set(FIELD_NAMES)), [])
        self.assertEqual(len(listed), len(set(listed)))

    def test_unlisted_field_falls_into_other(self):
        with patch.object(api, "_HYPERPARAMETER_SECTIONS",
                          {"Core Detection": ["interface_detect_distance_cutoff"]}):
            sections = api._hyperparameter_sections()
        self.assertIn("interface_detect_n_residue_cutoff", sections["Other"])
        self.assertNotIn("units", sections["Other"])

    def test_generated_docstring_documents_every_field(self):
        doc = api.set_hyperparameters.__doc__
        for name in DOCUMENTED:
            self.assertIn(f"- {name} (", doc)

    def test_generated_examples_use_real_fields(self):
        examples = api.set_hyperparameters.__doc__.split("Examples:", 1)[1]
        keywords = set(re.findall(r"(\w+)=", examples)) - {"workspace_path"}
        self.assertEqual(sorted(keywords - set(FIELD_NAMES)), [])

    def test_print_hyperparameters_lists_every_field(self):
        builder = PDBModelBuilder("1ABC")
        builder.set_hyperparameters()
        with contextlib.redirect_stdout(io.StringIO()):
            printed = builder.print_hyperparameters()
        for name in PRINTED:
            self.assertIn(f"  {name}: ", printed)


if __name__ == "__main__":
    unittest.main()

"""
tests/unit/model/pdb/test_io.py

Unit tests for model/pdb/io.py.

Purpose
-------
Verify that PDB/mmCIF downloading and parsing utilities behave correctly:
- download_pdb returns a valid path on success and raises RuntimeError on failure.
- parse_structure correctly uses MMCIFParser for `.cif` files and PDBParser for `.pdb`.

Testing strategy
----------------
- Mock network requests to avoid dependence on external RCSB servers.
- Create temporary fake PDB and CIF files to exercise parsing.
- Confirm correct Biopython parser is selected based on file extension.
- Validate exceptions are raised for failed downloads.

Assumptions
-----------
- Biopython is installed and provides PDBParser/MMCIFParser.
- requests is available and can be monkeypatched in tests.
- Test environment allows writing temporary files.

How to run
----------
pytest -q tests/unit/model/pdb/test_io.py
or
python -m unittest tests/unit/model/pdb/test_io.py
"""

import unittest
import tempfile
import os
import io
import gzip

from unittest import mock
from Bio.PDB import Structure
from ionerdss.model.pdb import io as pdb_io


class TestIO(unittest.TestCase):
    @mock.patch("requests.get")
    def test_download_pdb_success(self, mock_get):
        """Simulate successful download and decompression of CIF."""
        fake_content = gzip.compress(b"dummy cif content")
        mock_resp = mock.Mock()
        mock_resp.status_code = 200
        mock_resp.iter_content = lambda chunk_size: [fake_content]
        mock_get.return_value = mock_resp

        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = pdb_io.download_pdb("XXXX", tmpdir)
            self.assertTrue(os.path.exists(cif_path))
            with open(cif_path, "rb") as f:
                self.assertIn(b"dummy cif content", f.read())

    @mock.patch("requests.get")
    def test_download_pdb_failure_status(self, mock_get):
        """If RCSB returns error code, RuntimeError is raised."""
        mock_resp = mock.Mock()
        mock_resp.status_code = 404
        mock_get.return_value = mock_resp

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(RuntimeError):
                pdb_io.download_pdb("badid", tmpdir)

    @mock.patch("ionerdss.model.pdb.io.MMCIFParser")
    @mock.patch("ionerdss.model.pdb.io.PDBParser")
    def test_parse_structure_cif_and_pdb(self, mock_pdbparser, mock_cifparser):
        """Ensure correct parser is chosen based on file extension, no real parse needed."""
        # Arrange
        mock_cif = mock.Mock()
        mock_pdb = mock.Mock()
        mock_cifparser.return_value = mock_cif
        mock_pdbparser.return_value = mock_pdb
        mock_cif.get_structure.return_value = "CIF_STRUCT"
        mock_pdb.get_structure.return_value = "PDB_STRUCT"

        # Act
        result_cif = pdb_io.parse_structure("foo.cif")
        result_pdb = pdb_io.parse_structure("bar.pdb")

        # Assert: returns what the mock parsers returned
        self.assertEqual(result_cif, "CIF_STRUCT")
        self.assertEqual(result_pdb, "PDB_STRUCT")

        # Assert: correct parser was used
        mock_cif.get_structure.assert_called_once_with("foo", "foo.cif")
        mock_pdb.get_structure.assert_called_once_with("bar", "bar.pdb")


if __name__ == "__main__":
    unittest.main()

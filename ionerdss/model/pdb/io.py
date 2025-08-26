"""
model/pdb/io.py

PDB/CIF downloading and parsing utilities.

Purpose
-------
Provide a minimal I/O layer to fetch structures from the RCSB PDB archive
and parse them into Biopython `Structure` objects for downstream analysis
in the ionerdss pipeline.

Functions
---------
- download_pdb(pdb_id, save_dir)
    Download a gzipped PDB or mmCIF file directly from the RCSB PDB.
    Automatically decompresses `.gz` to `.cif`. Returns path to the
    decompressed file. Raises `RuntimeError` on failure.

- parse_structure(pdb_file)
    Detect file format (.cif or .pdb) and return a Biopython
    `Structure` object via `MMCIFParser` or `PDBParser`.

Design Notes
------------
- Downloads use HTTPS and a 30 second timeout.
- Currently only supports the default "assembly1" from RCSB.
- The caller must handle cleanup of files in `save_dir`.
- Parsing assumes input files are valid PDB or mmCIF.

Dependencies
------------
- Biopython (Bio.PDB)
- requests (for download)
- gzip, os (for file I/O)

Example
-------
>>> cif_path = download_pdb("1crn", "/tmp")
>>> structure = parse_structure(cif_path)
>>> len(list(structure.get_chains()))
1
"""

# model/pdb/io.py

import os
import gzip
import requests
from Bio.PDB import PDBList, MMCIFParser, PDBParser

def download_pdb(pdb_id, save_dir):
    """Download a PDB or CIF file given a PDB ID."""
    pdb_id = pdb_id.lower()
    decompressed_file = os.path.join(save_dir, f"{pdb_id}.cif")
    compressed_file = decompressed_file + ".gz"
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}-assembly1.cif.gz"

    try:
        response = requests.get(url, stream=True, timeout=30)
        if response.status_code == 200:
            with open(compressed_file, 'wb') as f_out:
                for chunk in response.iter_content(chunk_size=8192):
                    f_out.write(chunk)
            with gzip.open(compressed_file, 'rb') as f_in:
                with open(decompressed_file, 'wb') as f_out:
                    f_out.write(f_in.read())
            return decompressed_file
        else:
            raise RuntimeError(f"Failed to download PDB: {response.status_code}")
    except Exception as e:
        raise RuntimeError(f"Could not download {pdb_id}: {e}") from e

def parse_structure(pdb_file):
    """Parse .cif or .pdb file into Biopython Structure object."""
    parser = MMCIFParser(QUIET=True) if pdb_file.endswith('.cif') else PDBParser(QUIET=True)
    structure_id = os.path.basename(pdb_file).split('.')[0]
    return parser.get_structure(structure_id, pdb_file)

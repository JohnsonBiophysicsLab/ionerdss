"""
io.py

PDB/CIF downloading and parsing
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

#!/usr/bin/env python3
"""
Author: Sikao Guo
Modified by: Yue Moon Ying

This file is a .py version of the IPython notebooks for Coarse-grained model
generation from PDB code. 

This file includes a tutorial of how to automatically set up a coarse-grained model
of a protein complex given a pdb id as an input using `ionerdss`. The default we are
using is `8y7s`, which corresponds to benzaldehyde lyase mutant M6 from Herbiconiux
sp. SALV-R1. More information about the protein structure is available at:
https://www.rcsb.org/structure/8Y7S.


"""
import os
import urllib.request
import sys # for argv import
import ionerdss as ion
import subprocess # for running pymol
from os.path import expanduser, abspath # for expanding `~` in path

def parse_pdb_id(default_id="8y7s"):
    """
    Parses the command-line arguments to get a PDB ID, with a default id if
    the system argument not provided. We do not check whether the PDB actually
    exists in RCSB. If it does not exist, the program will throw an error when
    trying to fetch the PDB file.
    
    Returns:
        str (optional): The default lowercase PDB ID.
    
    Raises:
        SystemExit: If more than one argument is provided.
    """
    # accept command-line arguments
    # if provided more than 2 fields, send out a message to remind the user
    if len(sys.argv) > 2:
        print("Usage: ./coarse_grained_setup.py <PDB_code_string | OPTIONAL>")
        sys.exit(1)
    # if provided 2 fields, the second field is a PDB code
    elif len(sys.argv) == 2:
        return sys.argv[1].lower()
    else:
        return default_id.lower()

def setup_save_folder(pdb_id: str, folder_format: str = '~/Documents/{pdb_id}_dir') -> str:
    """
    Create and return an absolute save directory for the given PDB ID.

    Args:
        pdb_id (str): The 4-character PDB accession code.
        folder_format (str): A format string where '{pdb_id}' will be replaced by the PDB code.
                             Allows for custom directory structure (default: ~/Documents/{pdb_id}_dir).

    Returns:
        str: The absolute path to the created save folder.
    """
    folder_path = folder_format.format(pdb_id=pdb_id)
    abs_path = abspath(expanduser(folder_path))
    os.makedirs(abs_path, exist_ok=True)
    return abs_path

def download_pdb_thumbnail(pdb_id: str, save_folder: str) -> str:
    """
    Download the RCSB thumbnail image using only standard Python libraries.

    Args:
        pdb_id (str): The PDB code (4-character).
        save_folder (str): Directory where the image will be saved.

    Returns:
        str: Full path to the saved image.
    """
    image_url = f"https://cdn.rcsb.org/images/structures/{pdb_id}_assembly-1.jpeg"
    output_file = os.path.join(save_folder, f"{pdb_id}_thumbnail.jpeg")
    
    try:
        with urllib.request.urlopen(image_url, timeout=10) as response:
            with open(output_file, 'wb') as f:
                f.write(response.read())
        print(f"Thumbnail saved to {output_file}")
    except urllib.error.HTTPError as e:
        print(f"HTTP Error {e.code} while fetching thumbnail for {pdb_id}")
    except urllib.error.URLError as e:
        print(f"URL Error while fetching thumbnail: {e.reason}")
    except Exception as e:
        print(f"Unexpected error: {e}")
    
    return output_file


def render_pymol_image(save_folder: str, pdb_id: str) -> str:
    """
    Use PyMOL to render a static image of the regularized coarse-grained structure.

    Args:
        save_folder (str): The directory containing the PyMOL .pml script.
        pdb_id (str): The PDB accession code (used for naming/logging only).

    Returns:
        str: Full path to the generated image file, typically 'comparison_regularized.png'.

    Raises:
        FileNotFoundError: If PyMOL is not found in the system PATH.
        subprocess.CalledProcessError: If PyMOL exits with an error.
    """
    pml_path = os.path.join(save_folder, "visualize_regularized_coarse_grained.pml")
    try:
        subprocess.run(["pymol", "-cq", pml_path], check=True)
    except FileNotFoundError:
        print("PyMOL not found. Make sure 'pymol' is installed and in your PATH.")
    except subprocess.CalledProcessError as e:
        print(f"PyMOL execution failed: {e}")

    output_image = os.path.join(save_folder, f"comparison_regularized.png")
    if os.path.exists(output_image):
        print(f"PyMOL image saved to {output_image}")
    else:
        print("PyMOL did not produce the expected image.")
    return output_image

def process_pdb(pdb_id: str, folder_format: str = '~/Documents/{pdb_id}_dir'):
    """
    Main pipeline to:
      1. Create save folder
      2. Download RCSB thumbnail
      3. Initialize and process the PDBModel
      4. Run coarse-graining and regularization
      5. Render the resulting structure via PyMOL

    Args:
        pdb_id (str): The 4-character PDB accession code.
        folder_format (str): Optional format string for customizing output folder path.

    Notes:
        - Does not validate that the PDB actually exists in the RCSB.
        - PyMOL must be installed and accessible from the command line.
    """
    # get pdb id and working directory
    pdb_id = pdb_id.lower()
    save_folder = setup_save_folder(pdb_id, folder_format)
    download_pdb_thumbnail(pdb_id, save_folder)
    
    # initialize the pdb model given pdb id and save folder
    pdb_model = ion.PDBModel(pdb_id=pdb_id, save_dir=save_folder)

    # set up the coarse grain model from pdb model,
    # by detecting binding interfaces between chains based on atomic distances
    # for all tunable parameters see
    # https://ionerdss.readthedocs.io/en/latest/ionerdss.model.html#ionerdss.model.pdb_model.PDBModel.coarse_grain
    # set standard_output=True to see the determined interfaces
    pdb_model.coarse_grain(distance_cutoff=0.35, 
                        residue_cutoff=3,
                        show_coarse_grained_structure=False, 
                        save_pymol_script=False, 
                        standard_output=True)

    # identify, align, and regularize all molecular chains so that repeated chains share 
    # the same relative geometry. This method organizes molecule and interface objects 
    # accordingly and sets up reaction objects.
    # for all tunable parameters see
    # https://ionerdss.readthedocs.io/en/latest/ionerdss.model.html#ionerdss.model.pdb_model.PDBModel.regularize_repeated_chains
    pdb_model.regularize_repeated_chains(dist_threshold_intra=3.5, 
                                        dist_threshold_inter=3.5, 
                                        angle_threshold=25.0, 
                                        show_coarse_grained_structure=False, 
                                        save_pymol_script=True, 
                                        standard_output=False)
    
    render_pymol_image(save_folder, pdb_id)

"""
You can run from terminal:

```bash
python coarse_grain_pipeline.py
# or
python coarse_grain_pipeline.py 1hho
```

Let me know if you'd like to extend this with `argparse` (e.g., for optional parameters like `--folder_format`).

"""
if __name__ == "__main__":
    pdb_id = parse_pdb_id()
    print(f"Running coarse-grain pipeline for PDB ID: {pdb_id}")
    process_pdb(pdb_id)
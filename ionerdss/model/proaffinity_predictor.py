"""ProAffinity-GNN wrapper module for binding affinity prediction.

This module provides a self-contained interface to predict protein-protein binding
affinities using ProAffinity-GNN. It handles:
- PDB file download
- PDB file curation (filtering)
- PDBQT conversion using ADFR tools
- ProAffinity-GNN inference
- Unit conversion (k_BT to kJ/mol)
"""

import os
import sys
import subprocess
import urllib.request
import urllib.error
import numpy as np
import warnings

# Suppress known benign warnings from dependencies
# These warnings don't affect functionality and can be safely ignored
warnings.filterwarnings('ignore', category=FutureWarning, module='huggingface_hub.file_download', 
                       message='.*resume_download.*')
warnings.filterwarnings('ignore', message='.*Some weights of EsmModel were not initialized.*')
warnings.filterwarnings('ignore', message='.*You should probably TRAIN this model.*')


def download_pdb_direct(pdb_id, download_dir="pdbfiles", verbose=False):
    """
    Direct download from RCSB PDB using urllib.
    
    Args:
        pdb_id (str): PDB identifier (e.g., '1abc')
        download_dir (str): Directory to save downloaded files
        verbose (bool): Whether to print status messages
        
    Returns:
        str: Path to downloaded PDB file, or None if download failed
    """
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
    
    clean_id = pdb_id.strip().lower()
    url = f"https://files.rcsb.org/download/{clean_id}.pdb"
    filename = f"{clean_id}.pdb"
    filepath = os.path.join(download_dir, filename)
    
    try:
        if os.path.exists(filepath):
            if verbose: 
                print(f"Existed file: {clean_id.upper()} from {url}...")
            return filepath
        else:
            if verbose: 
                print(f"Downloading {clean_id.upper()} from {url}...")
            urllib.request.urlretrieve(url, filepath)
            if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
                if verbose: 
                    print(f"✓ Successfully downloaded: {filepath}")
                return filepath
            else:
                print(f"✗ Download failed or file is empty: {clean_id.upper()}")
                return None
            
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"✗ File not found (404): {clean_id.upper()}")
        else:
            print(f"✗ HTTP Error {e.code}: {clean_id.upper()}")
        return None
    except Exception as e:
        print(f"✗ Error downloading {clean_id.upper()}: {str(e)}")
        return None


def filter_pdb_file(input_pdb_path, output_pdb_path):
    """
    Filters a PDB file to include only 'ATOM' records with conventional residues.
    
    Args:
        input_pdb_path (str): Path to input PDB file
        output_pdb_path (str): Path to output filtered PDB file
    """
    conventional_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
        "DA", "DC", "DG", "DT",
        "A", "C", "G", "U", "T",
    }

    try:
        with open(input_pdb_path, 'r') as infile, \
             open(output_pdb_path, 'w') as outfile:
            for line in infile:
                if line.startswith("ATOM"):
                    residue_name = line[17:20].strip()
                    if residue_name in conventional_residues:
                        outfile.write(line)
                elif line.startswith("TER"):
                    outfile.write(line)
                    
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_pdb_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def pdb_to_pdbqt(pdbfile: str, adfr_path: str = None, ph: float = 7.4, verbose=False) -> str:
    """
    Converts a PDB file to PDBQT format using ADFR's prepare_receptor.
    
    Args:
        pdbfile (str): Path to input PDB file
        adfr_path (str): Path to ADFR installation (if None, searches in common locations)
        ph (float): pH value for adding hydrogens
        verbose (bool): Whether to print status messages
        
    Returns:
        str: Path to generated PDBQT file
        
    Raises:
        FileNotFoundError: If PDB file or prepare_receptor not found
        subprocess.CalledProcessError: If conversion fails
    """
    if not os.path.exists(pdbfile):
        raise FileNotFoundError(f"Input PDB file not found: {pdbfile}")
    
    # Locate prepare_receptor
    if adfr_path is None:
        # Search common locations
        possible_paths = [
            '/home/workspace/GitHub/ionerdss/ADFRsuite/bin/prepare_receptor',
            os.path.expanduser('~/ADFRsuite/bin/prepare_receptor'),
            '/usr/local/bin/prepare_receptor',
        ]
        for path in possible_paths:
            if os.path.exists(path):
                adfr_path = path
                break
    
    if adfr_path is None or not os.path.exists(adfr_path):
        raise FileNotFoundError(
            "prepare_receptor not found. Please install ADFR or specify adfr_path"
        )
    
    base_name, ext = os.path.splitext(pdbfile)
    pdbqtfile = f"{base_name}.pdbqt"
    tmp_pdb_file = f"{base_name}_tmp.pdb"

    try:
        if verbose:
            print(f"Filtering ATOM lines from {pdbfile} to {tmp_pdb_file}...")
        filter_pdb_file(pdbfile, tmp_pdb_file)
        if verbose:
            print("ATOM lines filtered successfully.")

        if verbose:
            print(f"Converting {tmp_pdb_file} to {pdbqtfile} using prepare_receptor...")
        
        command = [
            adfr_path,
            '-r', tmp_pdb_file,
            '-A', 'hydrogens',
            '-o', pdbqtfile
        ]
        
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        if verbose:
            print("Conversion successful.")
        if result.stdout and verbose:
            print("stdout:", result.stdout)
        if result.stderr and verbose:
            print("stderr:", result.stderr)

        if not os.path.exists(pdbqtfile):
            raise Exception(f"prepare_receptor ran but did not create {pdbqtfile}")

        return pdbqtfile

    except FileNotFoundError as e:
        print(f"Error: Command not found. Make sure prepare_receptor is available.")
        raise e
    except subprocess.CalledProcessError as e:
        print(f"Error during command execution:")
        print(f"Command: {' '.join(e.cmd)}")
        print(f"Return Code: {e.returncode}")
        print(f"Stdout: {e.stdout}")
        print(f"Stderr: {e.stderr}")
        raise e
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise e
    finally:
        if os.path.exists(tmp_pdb_file):
            os.remove(tmp_pdb_file)
            if verbose:
                print(f"Cleaned up temporary file: {tmp_pdb_file}")


def kbt_to_kj_mol(dG, T_kelvin=298.15):
    """
    Convert energy from k_B T units to kJ/mol.
    
    Args:
        dG (float or np.ndarray): Energy in k_B T units
        T_kelvin (float): Temperature in Kelvin
        
    Returns:
        float or np.ndarray: Energy in kJ/mol
    """
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    N_A = 6.02214076e23  # Avogadro's number in mol^-1
    
    kbt_joules = k_B * T_kelvin
    kbt_kj_mol = kbt_joules * N_A / 1000
    
    return kbt_kj_mol * np.array(dG)


def convert_pka_dG(pka, temperature=298.15):
    """
    Convert pKa to free energy (dG) in kJ/mol.
    
    Args:
        pka (float): pKa value
        temperature (float): Temperature in Kelvin
        
    Returns:
        float: Free energy in kJ/mol
    """
    K = 10**pka
    R = 8.314 / 1000  # kJ/(mol*K)
    dG = -R * temperature * np.log(K)
    return dG


def predict_proaffinity_binding_energy(pdb_id, chains, 
                                       model_weights_path=None,
                                       adfr_path=None,
                                       download_dir="pdbfiles",
                                       verbose=False):
    """
    Predict binding affinity using ProAffinity-GNN from a PDB ID.
    
    This is a self-contained function that handles the entire pipeline:
    1. Downloads PDB file
    2. Filters and curates the PDB
    3. Converts to PDBQT format
    4. Runs ProAffinity-GNN inference
    5. Returns result in kJ/mol
    
    Args:
        pdb_id (str): PDB identifier
        chains (str): Chain specification (e.g., 'A,B' or 'AB,CD')
        model_weights_path (str): Path to model.pkl file (default: same directory as this module)
        adfr_path (str): Path to ADFR prepare_receptor tool
        download_dir (str): Directory for PDB files
        verbose (bool): Whether to print status messages
        
    Returns:
        float: Predicted binding energy in kJ/mol (or np.nan if error)
    """
    try:
        # Step 1: Download PDB file
        if verbose:
            print(f"Processing {pdb_id.upper()}...")
            print("Step 1: Downloading PDB file...")
        pdbfile = download_pdb_direct(pdb_id, download_dir=download_dir, verbose=verbose)
        if pdbfile is None:
            raise Exception(f"Failed to download PDB: {pdb_id}")
        
        # Step 2: Convert to PDBQT
        if verbose:
            print("Step 2: Converting to PDBQT format...")
        pdbqtfile = pdb_to_pdbqt(pdbfile, adfr_path=adfr_path, verbose=verbose)
        
        # Step 3: Run ProAffinity-GNN inference
        if verbose:
            print("Step 3: Running ProAffinity-GNN inference...")
        
        # Import ProAffinity module from same directory
        from .ProAffinity_GNN_inference import run_proaffinity_inference
        
        # Set default model weights path to same directory as this module
        if model_weights_path is None:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            model_weights_path = os.path.join(current_dir, 'model.pkl')
        
        if not os.path.exists(model_weights_path):
            raise Exception(f"Model weights not found at: {model_weights_path}")
        
        # Run inference (returns energy in kJ/mol directly)
        dG_kj_mol = run_proaffinity_inference(pdbqtfile, chains, weights_path=model_weights_path, verbose=verbose)
        
        if verbose:
            print(f"Predicted binding energy: {dG_kj_mol:.4f} kJ/mol")
        
        return dG_kj_mol
        
    except Exception as e:
        print(f"Error in ProAffinity prediction for {pdb_id}: {str(e)}")
        return np.nan


# Convenience wrapper matching the Test.ipynb interface
def run_proaffinity_from_pdbid(pdb_id, chains, **kwargs):
    """
    Simplified interface matching Test.ipynb usage pattern.
    
    Args:
        pdb_id (str): PDB identifier
        chains (str): Chain specification
        **kwargs: Additional arguments passed to predict_proaffinity_binding_energy
        
    Returns:
        float: Binding energy in kJ/mol
    """
    return predict_proaffinity_binding_energy(pdb_id, chains, **kwargs)

"""PDBModel module for generating NERDSS molecule types and reactions from a PDB structure.

This module defines the `PDBModel` class, which extends the `Model` class to generate NERDSS molecule types,
reactions, and corresponding files from a PDB structure.
"""

import os
import gzip
import requests
import numpy as np
import math
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from Bio.PDB import PDBList, MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1
from scipy.spatial import KDTree
from sklearn.cluster import KMeans
from .components import MoleculeType, MoleculeInterface, Reaction, Model, MoleculeTemplate, CoarseGrainedMolecule, BindingInterfaceTemplate, BindingInterface, ReactionTemplate
from .pdb.geometry import rigid_transform_chains
from ..utils.coords import Coords


def apply_rigid_transform(R, t, point):
    """Apply rigid transformation (rotation R + translation t) to a point."""
    return R @ point + t


def extract_ca_coords(chain):
    """Extract CA coordinates from a BioPython Chain object.
    
    Parameters
    ----------
    chain : Bio.PDB.Chain
        BioPython Chain object
        
    Returns
    -------
    np.ndarray
        Array of CA coordinates with shape (N, 3)
    """
    coords = []
    for residue in chain:
        if is_aa(residue) and 'CA' in residue:
            coords.append(residue['CA'].coord)
    return np.array(coords)


def extract_aligned_ca_coords(chain1, chain2):
    """Extract CA coordinates from two chains, ensuring they have the same length.
    
    Parameters
    ----------
    chain1, chain2 : Bio.PDB.Chain
        BioPython Chain objects
        
    Returns
    -------
    tuple of np.ndarray
        Tuple of (chain1_coords, chain2_coords) with matching lengths
    """
    coords1 = extract_ca_coords(chain1)
    coords2 = extract_ca_coords(chain2)
    
    # Take the minimum length to ensure matching shapes
    min_len = min(len(coords1), len(coords2))
    if min_len == 0:
        raise ValueError(f"One of the chains has no CA atoms: chain1 has {len(coords1)}, chain2 has {len(coords2)}")
    
    return coords1[:min_len], coords2[:min_len]

class PDBModel(Model):
    """Handles the generation of NERDSS molecule types and reactions from a PDB structure.

    Attributes:
        pdb_file (str): Path to the PDB structure file.
        pdb_id (str): PDB ID of the structure.
        save_dir (str): Directory to save the output files.
    """

    def __init__(self, pdb_file: str = None, pdb_id: str = None, save_dir: str = None):
        """Initializes a PDBModel object.

        Args:
            pdb_file (str, optional): Path to the PDB structure file. Defaults to None.
            pdb_id (str, optional): PDB ID of the structure. Defaults to None.
            save_dir (str, optional): Directory to save output files. Defaults to None.

        Raises:
            ValueError: If neither `pdb_file` nor `pdb_id` is provided.
        """
        if save_dir.startswith("~"):
            save_dir = os.path.expanduser(save_dir)
        super().__init__(save_dir)
        self.pdb_file = pdb_file
        self.pdb_id = pdb_id
        self.save_dir = os.path.abspath(save_dir) if save_dir else os.getcwd()
        if not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)

        if not self.pdb_file and not self.pdb_id:
            raise ValueError("Either 'pdb_file' or 'pdb_id' must be provided.")

        if not self.pdb_file:
            self.pdb_file = self.download_pdb()

        if self.pdb_file and not pdb_id:
            self.pdb_id = os.path.basename(self.pdb_file).split('.')[0].lower()

        self.all_atoms_structure = self.pdb_parser()

        self.all_chains = []
        self.all_COM_chains_coords = []
        self.all_interfaces = []
        self.all_interfaces_coords = []
        self.all_interfaces_residues = []

        self.chains_map = {}  # Records the mapping of original chain IDs to molecular types used in NERDSS
        self.chains_group = []  # Groups chains with the same MOL_ID or entity_id or similar stucture as homologous

        # used to store the information of the molecules and interfaces for NERDSS model
        self.molecule_list: list[CoarseGrainedMolecule] = []
        self.molecules_template_list = []
        self.interface_list: list[BindingInterface] = []
        self.interface_template_list = []
        self.binding_chains_pairs = []
        self.binding_energies = []
        self.reaction_list: list[ReactionType] = []
        self.reaction_template_list = []

    def download_pdb(self) -> str:
        """Downloads the PDB structure file.

        Returns:
            str: Path to the downloaded PDB file.

        Raises:
            ValueError: If the PDB ID is invalid or the file cannot be retrieved.
        """
        if not self.pdb_id or len(self.pdb_id) != 4:
            raise ValueError("Invalid PDB ID. PDB IDs must be four characters long.")

        pdb_id_upper = self.pdb_id.upper()
        pdb_id_lower = self.pdb_id.lower()
        assembly_url = f"https://files.rcsb.org/download/{pdb_id_upper}-assembly1.cif.gz"
        compressed_file = os.path.join(self.save_dir, f"{pdb_id_lower}-assembly1.cif.gz")
        decompressed_file = os.path.join(self.save_dir, f"{pdb_id_lower}.cif")

        try:
            response = requests.get(assembly_url, stream=True)
            if response.status_code == 200:
                with open(compressed_file, 'wb') as file:
                    for chunk in response.iter_content(chunk_size=8192):
                        file.write(chunk)
                print(f"Successfully downloaded assembly file: {compressed_file}")

                with gzip.open(compressed_file, 'rb') as f_in:
                    with open(decompressed_file, 'wb') as f_out:
                        f_out.write(f_in.read())
                return decompressed_file
            else:
                print(f"Assembly file not available for {pdb_id_upper} (status code: {response.status_code})")
        except requests.RequestException as e:
            print(f"Failed to download assembly file for {pdb_id_upper}: {e}")

        try:
            print(f"Downloading the CIF file for {pdb_id_upper}...")
            pdbl = PDBList()
            pdbl.retrieve_pdb_file(pdb_id_upper, pdir=self.save_dir, file_format="mmCif")
            return decompressed_file
        except Exception as e:
            raise ValueError(f"Failed to download PDB file for {pdb_id_upper}: {e}")

    def pdb_parser(self):
        """Parses the .cif or .pdb file into a Biopython Structure object.

        Returns:
            Bio.PDB.Structure.Structure: The parsed structure containing all atoms.

        Raises:
            ValueError: If the file format is not .cif or .pdb.
        """
        if self.pdb_file.endswith('.cif'):
            parser = MMCIFParser(QUIET=True)
        elif self.pdb_file.endswith('.pdb'):
            parser = PDBParser(QUIET=True)
        else:
            raise ValueError("Unsupported file format. Only .cif and .pdb files are supported.")

        structure_id = os.path.basename(self.pdb_file).split('.')[0]
        structure = parser.get_structure(structure_id, self.pdb_file)
        return structure

    def coarse_grain(
        self,
        distance_cutoff=0.35,
        residue_cutoff=3,
        predict_affinity=False,
        adfr_path=None,
        show_coarse_grained_structure=False,
        save_pymol_script=False,
        standard_output=False
    ) -> None:
        """Coarse grains the PDB structure by detecting binding interfaces between chains based on atomic distances.

        Args:
            distance_cutoff (float, optional): Max distance (nm) for atoms to be considered in contact. Defaults to 0.35.
            residue_cutoff (int, optional): Minimum residue pair count to be considered a valid interface. Defaults to 3.
            predict_affinity (bool, optional): Whether to predict binding affinity using ProAffinity-GNN. Defaults to False.
            adfr_path (str, optional): Path to ADFR prepare_receptor tool. Required if predict_affinity is True. Defaults to None.
            show_coarse_grained_structure (bool, optional): Whether to visualize the coarse-grained structure. Defaults to False.
            save_pymol_script (bool, optional): Whether to save a PyMOL script for visualization. Defaults to False.
            standard_output (bool, optional): Whether to print detected interfaces. Defaults to False.
        """
        # self.all_chains = list(self.all_atoms_structure.get_chains())
        self.all_chains = sorted(self.all_atoms_structure.get_chains(), key=lambda chain: chain.id)
        self.all_COM_chains_coords = []
        self.all_interfaces = []
        self.all_interfaces_coords = []
        self.all_interfaces_residues = []
        self.all_chains_radius = []

        self.interface_energies: list[float] = []
        self.all_interface_energies: list[list[float]] = []

        # Initialize interface lists
        num_chains = len(self.all_chains)
        for _ in range(num_chains):
            self.all_interfaces.append([])
            self.all_interfaces_coords.append([])
            self.all_interfaces_residues.append([])
            self.all_interface_energies.append([])

        # Calculate the center of mass (COM) for each chain
        for chain in self.all_chains:
            atom_coords = [atom.coord for residue in chain for atom in residue if is_aa(residue)]
            if not atom_coords:
                self.all_COM_chains_coords.append(None)
                continue

            # Calculate the COM
            avg_coords = np.mean(atom_coords, axis=0)
            self.all_COM_chains_coords.append(Coords(*avg_coords))

            # Calculate the radius of gyration of each chain based on its COM
            if len(atom_coords) > 1:
                distances = np.linalg.norm(atom_coords - avg_coords, axis=1)
                radius = np.sqrt(np.mean(distances ** 2))
            else:
                radius = 0.0
            self.all_chains_radius.append(radius)

        # Helper function to compute bounding box for a chain
        def compute_bounding_box(chain):
            atom_coords = np.array([atom.coord for residue in chain for atom in residue if is_aa(residue)])
            if atom_coords.size == 0:
                return None, None
            min_coords = np.min(atom_coords, axis=0)
            max_coords = np.max(atom_coords, axis=0)
            return min_coords, max_coords
        
        # Precompute bounding boxes for all chains
        bounding_boxes = [compute_bounding_box(chain) for chain in self.all_chains]

        # Store chain pairs that need affinity prediction
        affinity_prediction_pairs = []
        pair_to_indices = {}  # Maps (i,j) to index in affinity_prediction_pairs
        
        # Helper function to process a pair of chains
        def process_chain_pair(i, j):
            if self.all_COM_chains_coords[i] is None or self.all_COM_chains_coords[j] is None:
                return

            min_box1, max_box1 = bounding_boxes[i]
            min_box2, max_box2 = bounding_boxes[j]

            # Skip if bounding boxes are farther apart than the cutoff distance
            if np.any(min_box2 > max_box1 + distance_cutoff * 10) or np.any(max_box2 < min_box1 - distance_cutoff * 10):
                return
            
            chain1 = self.all_chains[i]
            chain2 = self.all_chains[j]

            atom_coords_chain1 = []
            ca_coords_chain1 = []
            residue_ids_chain1 = []
            residue_types_chain1 = []
            atom_coords_chain2 = []
            ca_coords_chain2 = []
            residue_ids_chain2 = []
            residue_types_chain2 = []

            for residue1 in chain1:
                if not is_aa(residue1) or 'CA' not in residue1:
                    continue
                res_type1 = residue1.get_resname().upper()
                for atom1 in residue1:
                    atom_coords_chain1.append(atom1.coord)
                    ca_coords_chain1.append(residue1['CA'].coord)
                    residue_ids_chain1.append(residue1.id[1])
                    residue_types_chain1.append(res_type1)

            for residue2 in chain2:
                if not is_aa(residue2) or 'CA' not in residue2:
                    continue
                res_type2 = residue2.get_resname().upper()
                for atom2 in residue2:
                    atom_coords_chain2.append(atom2.coord)
                    ca_coords_chain2.append(residue2['CA'].coord)
                    residue_ids_chain2.append(residue2.id[1])
                    residue_types_chain2.append(res_type2)

            if len(ca_coords_chain1) == 0 or len(ca_coords_chain2) == 0:
                return

            # Build KDTree for chain2
            tree = KDTree(atom_coords_chain2)
            indices = tree.query_ball_point(atom_coords_chain1, r=distance_cutoff * 10)

            interface1 = []
            interface1_coords = []
            interface1_types = []
            interface2 = []
            interface2_coords = []
            interface2_types = []

            # Collect interface residues based on KDTree results
            for idx1, neighbors in enumerate(indices):
                if neighbors:
                    if residue_ids_chain1[idx1] not in interface1:
                        interface1.append(residue_ids_chain1[idx1])
                        interface1_coords.append(ca_coords_chain1[idx1])
                        interface1_types.append(residue_types_chain1[idx1])

                    for idx2 in neighbors:
                        if residue_ids_chain2[idx2] not in interface2:
                            interface2.append(residue_ids_chain2[idx2])
                            interface2_coords.append(ca_coords_chain2[idx2])
                            interface2_types.append(residue_types_chain2[idx2])
            
            # Store results if any interfaces were found
            if len(interface1) >= residue_cutoff and len(interface2) >= residue_cutoff:
                # Store interface information for later processing
                avg_coords1 = np.mean(interface1_coords, axis=0)
                avg_coords2 = np.mean(interface2_coords, axis=0)
                
                # Store interface data temporarily (energy will be added later)
                result_data = {
                    'i': i,
                    'j': j,
                    'interface1': sorted(interface1),
                    'interface2': sorted(interface2),
                    'avg_coords1': avg_coords1,
                    'avg_coords2': avg_coords2,
                }
                
                # If affinity prediction is enabled, store chain pair info for batch processing
                if predict_affinity:
                    chain1_id = self.all_chains[i].id
                    chain2_id = self.all_chains[j].id
                    pair_idx = len(affinity_prediction_pairs)
                    affinity_prediction_pairs.append({
                        'pdb_file': self.pdb_file,
                        'chains': f"{chain1_id},{chain2_id}",
                        'result_data': result_data
                    })
                    pair_to_indices[(i, j)] = pair_idx
                else:
                    # Use default energy immediately
                    result_data['total_energy'] = -16 * 8.314/1000 * 298  # -16RT in kJ/mol
                    
                return result_data
            return None

        # Parallelize chain pair processing to detect interfaces
        with ThreadPoolExecutor() as executor:
            futures = [executor.submit(process_chain_pair, i, j) for i in range(num_chains - 1) for j in range(i + 1, num_chains)]
            interface_results = []
            for future in futures:
                result = future.result()
                if result is not None:
                    interface_results.append(result)
        
        # Batch predict affinities if needed
        if predict_affinity and affinity_prediction_pairs:
            # Print informational message about ProAffinity-GNN usage
            print("\n" + "="*80)
            print("NOTE: Using simplified ProAffinity-GNN predictor")
            print("="*80)
            print("This is an easy-to-use version that skips sequence alignment with")
            print("canonical FASTA sequences. For better accuracy and advanced options,")
            print("please visit: https://github.com/legendzzy/ProAffinity-GNN")
            print("="*80 + "\n")
            

            if standard_output:
                print(f"Running batch affinity prediction for {len(affinity_prediction_pairs)} chain pairs...")
                print("Loading ESM2 model once for all predictions...")
            
            try:
                from .proaffinity_predictor import predict_proaffinity_binding_energy_batch
                
                # Run batch predictions
                binding_energies = predict_proaffinity_binding_energy_batch(
                    predictions_list=affinity_prediction_pairs,
                    adfr_path=adfr_path,
                    verbose=standard_output
                )
                
                # Assign predicted energies back to results
                for pair_info, binding_energy in zip(affinity_prediction_pairs, binding_energies):
                    result_data = pair_info['result_data']
                    if np.isnan(binding_energy):
                        # Fall back to default if prediction fails
                        result_data['total_energy'] = -16 * 8.314/1000 * 298  # -16RT in kJ/mol
                        if standard_output:
                            print(f"Warning: ProAffinity prediction failed for chains {pair_info['chains']}, using default energy")
                    else:
                        result_data['total_energy'] = binding_energy
                        if standard_output:
                            print(f"Predicted binding energy for chains {pair_info['chains']}: {binding_energy:.2f} kJ/mol")
                            
            except Exception as e:
                # Fall back to default energies if batch prediction fails
                if standard_output:
                    print(f"Warning: Batch affinity prediction failed: {e}")
                    print("Using default energies for all chain pairs")
                for pair_info in affinity_prediction_pairs:
                    pair_info['result_data']['total_energy'] = -16 * 8.314/1000 * 298
        
        # Now populate the interface lists with the results
        for result_data in interface_results:
            i = result_data['i']
            j = result_data['j']
            total_energy = result_data['total_energy']
            
            self.all_interfaces[i].append(self.all_chains[j].id)
            self.all_interfaces_coords[i].append(Coords(*result_data['avg_coords1']))
            self.all_interfaces_residues[i].append(result_data['interface1'])
            self.all_interface_energies[i].append(total_energy)
            
            self.all_interfaces[j].append(self.all_chains[i].id)
            self.all_interfaces_coords[j].append(Coords(*result_data['avg_coords2']))
            self.all_interfaces_residues[j].append(result_data['interface2'])
            self.all_interface_energies[j].append(total_energy)

        for i in range(num_chains):
            sorted_indices = sorted(range(len(self.all_interfaces[i])), key=lambda k: self.all_interfaces[i][k])
            self.all_interfaces[i] = [self.all_interfaces[i][k] for k in sorted_indices]
            self.all_interfaces_coords[i] = [self.all_interfaces_coords[i][k] for k in sorted_indices]
            self.all_interfaces_residues[i] = [self.all_interfaces_residues[i][k] for k in sorted_indices]
            self.all_interface_energies[i] = [self.all_interface_energies[i][k] for k in sorted_indices]

        # Print detected interfaces
        if standard_output:
            print("Binding interfaces detected:")
            for i, chain in enumerate(self.all_chains):
                print(f"Chain {chain.id}:")
                print(f"  Center of Mass (COM): {self.all_COM_chains_coords[i]}")
                print(f"  Interfaces: {self.all_interfaces[i]}")
                print("  Interface Coordinates: ")
                for j, interface_coord in enumerate(self.all_interfaces_coords[i]):
                    print(f"    {interface_coord}")
                    print(f"    Interface Energy: {self.all_interface_energies[i][j]:.2f} kJ/mol")

        # Save PyMOL script
        if save_pymol_script:
            self.save_original_coarse_grained_structure()

        # Plot the original coarse-grained structure
        if show_coarse_grained_structure:
            self.plot_original_coarse_grained_structure()

    def plot_original_coarse_grained_structure(self):
        """Visualizes the original coarse-grained structure, showing each chain’s COM and interface coordinates before regularization."""
        all_points = []
        chain_ids = []
        for chain in self.all_chains:
            chain_id = chain.id
            chain_ids.append(chain_id)
            com_coord = self.all_COM_chains_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])]
            interface_coords = self.all_interfaces_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])]
            points = []
            points.append([com_coord.x, com_coord.y, com_coord.z])
            for interface_coord in interface_coords:
                points.append([interface_coord.x, interface_coord.y, interface_coord.z])
            all_points.append(points)
        self.plot_points_3d(all_points, chain_ids)

    def plot_points_3d(self, points, chain_ids=None):
        """Plots sets of 3D points for multiple chains in a single 3D Matplotlib figure.

        Args:
            points (list): A list of arrays, each of shape (N, 3), representing a chain’s COM + interface sites.
            chain_ids (list, optional): A list of labels for each chain. Defaults to None.
        """

        # Prepare a 3D figure
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection="3d")

        # Generate a color cycle for different chains
        colors = plt.cm.get_cmap("tab10", len(points))

        for i, chain in enumerate(points):
            # chain[0] is the center of mass (COM)
            com = chain[0]
            # The remaining points in the chain are interface points
            interfaces = chain[1:]

            # Pick a color for this chain
            color = colors(i)

            # Plot the COM
            ax.scatter(com[0], com[1], com[2],
                    color=color,
                    s=70,  # size of marker
                    marker="o",
                    label=f"Chain {chain_ids[i]} COM" if chain_ids != None else None)

            # Plot interfaces and lines to the COM
            for j, interface in enumerate(interfaces):
                # Plot the interface point
                ax.scatter(interface[0], interface[1], interface[2],
                        color=color,
                        s=50,
                        marker="^")  # or any shape you like

                # Draw a line from the COM to this interface
                xs = [com[0], interface[0]]
                ys = [com[1], interface[1]]
                zs = [com[2], interface[2]]
                ax.plot(xs, ys, zs, color=color, linewidth=1)

        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.set_title("Original Coarse-Grained Structure")

        ax.legend(loc="best")

        plt.tight_layout()
        plt.show()

    def save_original_coarse_grained_structure(self, output_cif: str = "original_coarse_grained_structure.cif", pymol_script: str = "original_visualize_coarse_grained.pml"):
        """Saves the original coarse-grained structure (COM and interface coordinates for each chain)
        to a CIF file and generates a PyMOL script for quick visualization.

        Args:
            output_cif (str, optional): Output .cif filename. Defaults to "original_coarse_grained_structure.cif".
            pymol_script (str, optional): Output .pml filename for PyMOL. Defaults to "original_visualize_coarse_grained.pml".
        """
        output_cif = os.path.join(self.save_dir, output_cif)
        pymol_script = os.path.join(self.save_dir, pymol_script)

        with open(output_cif, 'w') as cif_file:
            atom_id = 1

            # Write CIF header
            cif_file.write("# Coarse-grained structure CIF file\n")
            cif_file.write("data_coarse_grained\n")
            cif_file.write("_audit_conform_dict.text 'Original coarse-grained model generated by ionerdss'\n")
            cif_file.write("loop_\n")
            cif_file.write("_atom_site.group_PDB\n")
            cif_file.write("_atom_site.id\n")
            cif_file.write("_atom_site.label_atom_id\n")
            cif_file.write("_atom_site.label_comp_id\n")
            cif_file.write("_atom_site.label_asym_id\n")
            cif_file.write("_atom_site.Cartn_x\n")
            cif_file.write("_atom_site.Cartn_y\n")
            cif_file.write("_atom_site.Cartn_z\n")
            cif_file.write("_atom_site.occupancy\n")
            cif_file.write("_atom_site.B_iso_or_equiv\n")
            cif_file.write("_atom_site.type_symbol\n")

            # Write COM atoms for each chain
            for i, chain in enumerate(self.all_chains):
                if not self.all_COM_chains_coords[i]:
                    continue
                com = self.all_COM_chains_coords[i]
                cif_file.write(
                    f"ATOM  {atom_id:5d}  COM  MOL {chain.id}  "
                    f"{com.x:8.3f} {com.y:8.3f} {com.z:8.3f}  1.00  0.00  C\n"
                )
                atom_id += 1

                # Write interface atoms for the current chain
                for j, interface_coord in enumerate(self.all_interfaces_coords[i]):
                    cif_file.write(
                        f"ATOM  {atom_id:5d}  INT  MOL {chain.id}  "
                        f"{interface_coord.x:8.3f} {interface_coord.y:8.3f} {interface_coord.z:8.3f}  1.00  0.00  O\n"
                    )
                    atom_id += 1

        print(f"Coarse-grained structure saved to {output_cif}.")

        # Generate PyMOL script for visualization
        with open(pymol_script, 'w') as pml_file:
            pml_file.write("# PyMOL script to visualize coarse-grained structure\n")
            pml_file.write(f"load {output_cif}, coarse_grained\n")
            pml_file.write("hide everything\n")
            pml_file.write("show spheres, name COM\n")
            pml_file.write("show spheres, name INT\n")
            pml_file.write("set sphere_scale, 1.0\n")
            pml_file.write("color red, name COM\n")
            pml_file.write("color blue, name INT\n")
            
            # Create pseudo-atoms for COM and interfaces and draw lines
            atom_index = 1
            for i, chain in enumerate(self.all_chains):
                com = self.all_COM_chains_coords[i]
                if not com:
                    continue
                # Make a pseudoatom for the chain's COM
                pml_file.write(
                    f"pseudoatom com_{chain.id}, pos=[{com.x:.3f}, {com.y:.3f}, {com.z:.3f}], color=red\n"
                )
                
                # For each interface, create a pseudoatom and connect it to the COM
                for j, interface_coord in enumerate(self.all_interfaces_coords[i], start=1):
                    pml_file.write(
                        f"pseudoatom int_{chain.id}_{j}, pos=[{interface_coord.x:.3f}, "
                        f"{interface_coord.y:.3f}, {interface_coord.z:.3f}], color=blue\n"
                    )
                    # Use f-strings so {atom_index} is replaced numerically
                    pml_file.write(f"distance line{atom_index}, com_{chain.id}, int_{chain.id}_{j}\n")
                    pml_file.write(f"set dash_width, 4, line{atom_index}\n")
                    pml_file.write(f"set dash_gap, 0.5, line{atom_index}\n")
                    atom_index += 1

            pml_file.write("set sphere_transparency, 0.2\n")
            pml_file.write("bg_color white\n")
            pml_file.write("zoom all\n")

            pml_file.write(f"load {self.pdb_file}\n")
            pml_file.write("spectrum chain\n")
            pdb_file = os.path.basename(self.pdb_file).split('.')[0]
            pml_file.write(f"hide everything, {pdb_file}\n")
            pml_file.write(f"show cartoon, {pdb_file}\n")
            pml_file.write(f"set cartoon_transparency, 0.7, {pdb_file}\n")
            save_figure = os.path.join(self.save_dir, "comparison_initial.png")
            pml_file.write(f"png {save_figure}, 800, 800, 150, 1\n")

        print(f"PyMOL script saved to {pymol_script}. Run 'pymol {pymol_script}' to visualize the coarse-grained structure.")

    def regularize_homologous_chains(self, dist_threshold_intra=3.5, dist_threshold_inter=3.5, angle_threshold=25.0, show_coarse_grained_structure=False, save_pymol_script=False, standard_output=False):
        """
        Aligns and regularizes all molecular chains so that homologous chains share 
        the same relative geometry. This method organizes molecule and interface objects 
        accordingly and sets up reaction objects.

        Args:
            dist_threshold_intra (float): Distance threshold for intra-chain similarity. Defaults to 3.5 angstrom.
            dist_threshold_inter (float): Distance threshold for inter-chain similarity. Defaults to 3.5 angstrom.
            angle_threshold (float): Angle threshold for similarity. Defaults to 25.0 degree.
            show_coarse_grained_structure (bool): Whether to visualize the regularized coarse-grained structure. Defaults to False.
            save_pymol_script (bool): Whether to save a PyMOL script for visualization. Defaults to False.
            standard_output (bool): Whether to print detected interfaces. Defaults to False.
        """
        self.identify_homologous_chains()
        if not self.chains_group:
            self._assign_original_chain_ids()
        for group in self.chains_group:
            group.sort()
        self.chains_group.sort()
        print("Homologous chain groups identified:")
        print(self.chains_group)

        self.molecule_list: list[CoarseGrainedMolecule] = []
        self.molecules_template_list = []
        self.interface_list: list[BindingInterface] = []
        self.interface_template_list = []
        self.interface_signatures = []

        for group in self.chains_group:
            # print(f"Start parsing chain group / molecule template {group}")
            mol_temp_name = self.chains_map[group[0]]
            is_existing_mol_temp, idx = self._is_existing_mol_temp(mol_temp_name)
            if is_existing_mol_temp:
                # print(f"This is an existed mol template {mol_temp_name}")
                molecule_template = self.molecules_template_list[idx]
            else:
                molecule_template = MoleculeTemplate(mol_temp_name)
                # print(f"New mol template {mol_temp_name} is created.")
                self.molecules_template_list.append(molecule_template)

            for j, chain_id in enumerate(group):
                # print(f"Start parsing chain / molecule {chain_id}")
                mol_name = chain_id
                is_existing_mol, mol_index = self._is_existing_mol(mol_name)
                if is_existing_mol:
                    # print(f"This is an existing molecule {mol_name}")
                    molecule: CoarseGrainedMolecule = self.molecule_list[mol_index]
                    molecule.radius = self.all_chains_radius[self.all_chains.index([chain for chain in self.all_chains if chain.id == mol_name][0])]
                    molecule.diffusion_translation, molecule.diffusion_rotation = self._compute_diffusion_constants_nm_us(molecule.radius / 10.0)
                    molecule.my_template.diffusion_translation, molecule.my_template.diffusion_rotation = molecule.diffusion_translation, molecule.diffusion_rotation
                else:
                    molecule = CoarseGrainedMolecule(mol_name)
                    # print(f"New molecule {mol_name} is created.")
                    molecule.my_template = molecule_template
                    molecule.coord = self.all_COM_chains_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == mol_name][0])]
                    molecule.radius = self.all_chains_radius[self.all_chains.index([chain for chain in self.all_chains if chain.id == mol_name][0])]
                    molecule.diffusion_translation, molecule.diffusion_rotation = self._compute_diffusion_constants_nm_us(molecule.radius / 10.0)
                    self.molecule_list.append(molecule)
                    molecule_template.radius = molecule.radius
                    molecule_template.diffusion_translation, molecule_template.diffusion_rotation = molecule.diffusion_translation, molecule.diffusion_rotation
                
                # loop the interface of this chain (molecule)
                for i, interface_id in enumerate(self.all_interfaces[self.all_chains.index([chain for chain in self.all_chains if chain.id == mol_name][0])]):
                    A = mol_name
                    B = interface_id # this is the chain name of the partner
                    partner_mol_template_name = self.chains_map[B]
                    # print(f"Parsing the interface {interface_id} for molecule {mol_name}; its binding partner is molecule {B} via its interface {A}")
                    is_existing_mol_temp, idx = self._is_existing_mol_temp(partner_mol_template_name)
                    if is_existing_mol_temp:
                        # print(f"molecule {B} already has its template created.")
                        partner_molecule_template = self.molecules_template_list[idx]
                    else:
                        partner_molecule_template = MoleculeTemplate(partner_mol_template_name)
                        # print(f"new mol template {partner_mol_template_name} created for molecule {B}.")
                        self.molecules_template_list.append(partner_molecule_template)

                    is_existing_mol, partner_mol_index = self._is_existing_mol(B)
                    if is_existing_mol:
                        # print(f"molecule {B} is already created.")
                        partner_molecule = self.molecule_list[partner_mol_index]
                    else:
                        partner_molecule = CoarseGrainedMolecule(B)
                        # print(f"New molecule {B} is created.")
                        partner_molecule.my_template = partner_molecule_template
                        partner_molecule.coord = self.all_COM_chains_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == B][0])]
                        self.molecule_list.append(partner_molecule)

                    COM_A = self.all_COM_chains_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == A][0])]
                    I_A = self.all_interfaces_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == A][0])][i]
                    COM_B = self.all_COM_chains_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == B][0])]
                    for k, partner_interface_id in enumerate(self.all_interfaces[self.all_chains.index([chain for chain in self.all_chains if chain.id == B][0])]):
                        if partner_interface_id == A:
                            I_B = self.all_interfaces_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == B][0])][k]
                            R_B = self.all_interfaces_residues[self.all_chains.index([chain for chain in self.all_chains if chain.id == B][0])][k]
                            E_B = self.all_interface_energies[self.all_chains.index([chain for chain in self.all_chains if chain.id == B][0])][k]
                            break

                    signature = {
                        "dA": np.linalg.norm([(COM_A - I_A).x, (COM_A - I_A).y, (COM_A - I_A).z]),
                        "dB": np.linalg.norm([(COM_B - I_B).x, (COM_B - I_B).y, (COM_B - I_B).z]),
                        "dAB": np.linalg.norm([(I_A - I_B).x, (I_A - I_B).y, (I_A - I_B).z]),
                        "thetaA": self._calc_angle(COM_A, I_A, I_B),
                        "thetaB": self._calc_angle(COM_B, I_B, I_A)
                    }

                    # print the signature
                    # print(f"Parsing signature: {signature}")

                    is_existing_sig = False

                    for existing_sig in self.interface_signatures:
                        if self._sig_are_similar(signature, existing_sig, dist_threshold_intra, dist_threshold_inter, angle_threshold):
                            is_existing_sig = True
                            break

                    if not is_existing_sig:
                        # print("this is a new signature. added to list.")
                        self.interface_signatures.append(signature)
                        signature_conjugated = {
                            "dA": signature["dB"],
                            "dB": signature["dA"],
                            "dAB": signature["dAB"],
                            "thetaA": signature["thetaB"],
                            "thetaB": signature["thetaA"]
                        }
                        self.interface_signatures.append(signature_conjugated)
                        # print(f"the conjugated signature: {signature_conjugated} is also added to the list.")

                        # build the interface template pairs for both molecule templates, need to check if this is homo dimerization or hetero
                        is_homo = False
                        if self.chains_map[A] != self.chains_map[B]:
                            pass
                        else:
                            if abs(signature["dA"] - signature["dB"]) > dist_threshold_intra or abs(signature["thetaA"] - signature["thetaB"]) > angle_threshold:
                                pass
                            else:
                                is_homo = True

                        if is_homo:
                            # only need to build the interface template once
                            interface_template_id_prefix = self.chains_map[A]

                            # determine the sufffix of this interface_template
                            tmp_count = 1
                            for interface_temp in molecule_template.interface_template_list:
                                interface_temp_id = interface_temp.name
                                if interface_temp_id.startswith(interface_template_id_prefix):
                                    tmp_count += 1

                            interface_template_id_suffix = str(tmp_count)
                            interface_template_id = interface_template_id_prefix + interface_template_id_suffix
                            interface_template = BindingInterfaceTemplate(interface_template_id)
                            interface_template.signature = signature
                            if j == 0:
                                interface_template.coord = self.all_interfaces_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])][i] - molecule.coord
                            else:
                                # align the current chain to the first chain in the group, then get the relative position of interface to COM
                                chain1 = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == group[0]][0])]
                                chain2 = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])]
                                R, t = rigid_transform_chains(chain2, chain1)
                                Q = []
                                Q_COM_coord = self.all_COM_chains_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])]
                                Q.append([Q_COM_coord.x, Q_COM_coord.y, Q_COM_coord.z])
                                temp_coord = self.all_interfaces_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])][i]
                                Q.append([temp_coord.x, temp_coord.y, temp_coord.z])
                                Q2 = []
                                for point in Q:
                                    transformed_point = apply_rigid_transform(R, t, np.array(point))
                                    Q2.append(transformed_point)
                                interface_template.coord = Coords(Q2[1][0] - Q2[0][0], Q2[1][1] - Q2[0][1], Q2[1][2] - Q2[0][2])
                            molecule_template.interface_template_list.append(interface_template)
                            self.interface_template_list.append(interface_template)
                            partner_interface_template = interface_template
                            partner_molecule_template = molecule_template
                        else:
                            # add interface template 1
                            interface_template_id_prefix = self.chains_map[B]

                            # determine the sufffix of this interface_template
                            tmp_count = 1
                            for interface_temp in molecule_template.interface_template_list:
                                interface_temp_id = interface_temp.name
                                if interface_temp_id.startswith(interface_template_id_prefix):
                                    tmp_count += 1

                            interface_template_id_suffix = str(tmp_count)
                            interface_template_id = interface_template_id_prefix + interface_template_id_suffix
                            interface_template = BindingInterfaceTemplate(interface_template_id)
                            interface_template.signature = signature
                            if j == 0:
                                interface_template.coord = self.all_interfaces_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])][i] - molecule.coord
                            else:
                                # align the current chain to the first chain in the group, then get the relative position of interface to COM
                                chain1 = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == group[0]][0])]
                                chain2 = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])]
                                R, t = rigid_transform_chains(chain2, chain1)
                                Q = []
                                Q_COM_coord = self.all_COM_chains_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])]
                                Q.append([Q_COM_coord.x, Q_COM_coord.y, Q_COM_coord.z])
                                temp_coord = self.all_interfaces_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])][i]
                                Q.append([temp_coord.x, temp_coord.y, temp_coord.z])
                                Q2 = []
                                for point in Q:
                                    transformed_point = apply_rigid_transform(R, t, np.array(point))
                                    Q2.append(transformed_point)
                                interface_template.coord = Coords(Q2[1][0] - Q2[0][0], Q2[1][1] - Q2[0][1], Q2[1][2] - Q2[0][2])
                            molecule_template.interface_template_list.append(interface_template)
                            self.interface_template_list.append(interface_template)

                            # add interface template 2
                            interface_template_id_prefix = self.chains_map[A]

                            # determine the sufffix of this interface_template
                            tmp_count = 1
                            for interface_temp in molecule_template.interface_template_list:
                                interface_temp_id = interface_temp.name
                                if interface_temp_id.startswith(interface_template_id_prefix):
                                    tmp_count += 1

                            interface_template_id_suffix = str(tmp_count)
                            interface_template_id = interface_template_id_prefix + interface_template_id_suffix
                            partner_interface_template = BindingInterfaceTemplate(interface_template_id)
                            partner_interface_template.signature = signature_conjugated
                            B_group = None
                            for g in self.chains_group:
                                if B in g:
                                    B_group = g

                            if B == B_group[0]:
                                partner_interface_template.coord = I_B - partner_molecule.coord
                            else:
                                # align the current chain to the first chain in the group, then get the relative position of interface to COM
                                chain1 = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == B_group[0]][0])]
                                chain2 = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == B][0])]
                                R, t = rigid_transform_chains(chain2, chain1)
                                Q = []
                                Q_COM_coord = self.all_COM_chains_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == B][0])]
                                Q.append([Q_COM_coord.x, Q_COM_coord.y, Q_COM_coord.z])
                                temp_coord = I_B
                                Q.append([temp_coord.x, temp_coord.y, temp_coord.z])
                                Q2 = []
                                for point in Q:
                                    transformed_point = apply_rigid_transform(R, t, np.array(point))
                                    Q2.append(transformed_point)
                                partner_interface_template.coord = Coords(Q2[1][0] - Q2[0][0], Q2[1][1] - Q2[0][1], Q2[1][2] - Q2[0][2])
                            partner_molecule_template.interface_template_list.append(partner_interface_template)
                            self.interface_template_list.append(partner_interface_template)

                    else:
                        # print("this is an existing signature. using the existing interface template.")
                        # find the interface_template and partner_interface_template
                        interface_template = None
                        partner_interface_template = None
                        signature_conjugated = {
                            "dA": signature["dB"],
                            "dB": signature["dA"],
                            "dAB": signature["dAB"],
                            "thetaA": signature["thetaB"],
                            "thetaB": signature["thetaA"]
                        }
                        for mol_temp in self.molecules_template_list:
                            for interface_temp in mol_temp.interface_template_list:
                                if self._sig_are_similar(signature, interface_temp.signature, dist_threshold_intra, dist_threshold_inter, angle_threshold):
                                    interface_template = interface_temp
                                    molecule_template = mol_temp
                                    # print(f"using {mol_temp.name} - {interface_temp.name}")
                                    break
                        for mol_temp in self.molecules_template_list:
                            for interface_temp in mol_temp.interface_template_list:
                                if self._sig_are_similar(signature_conjugated, interface_temp.signature, dist_threshold_intra, dist_threshold_inter, angle_threshold):
                                    partner_interface_template = interface_temp
                                    partner_molecule_template = mol_temp
                                    # print(f"using {mol_temp.name} - {interface_temp.name}")
                                    break

                    # build the interfaces for molecules, link the interface template to interface

                    is_existing_interface, _ = self._is_existing_interface(interface_id, molecule)

                    if not is_existing_interface:
                        # print(f"Creating new interface {interface_id} for molecule {mol_name}")
                        # create the interface
                        interface = BindingInterface(B)
                        interface.my_template = interface_template
                        interface.coord = self.all_interfaces_coords[self.all_chains.index([chain for chain in self.all_chains if chain.id == A][0])][i]
                        interface.my_residues = self.all_interfaces_residues[self.all_chains.index([chain for chain in self.all_chains if chain.id == A][0])][i]
                        interface.energy = self.all_interface_energies[self.all_chains.index([chain for chain in self.all_chains if chain.id == A][0])][i]
                        interface.my_template.energy = interface.energy
                        self.interface_list.append(interface)
                        molecule.interface_list.append(interface)

                        # print(f"Creating new interface {A} for partner molecule {B}")
                        # create the interface for the partner molecule
                        partner_interface = BindingInterface(A)
                        partner_interface.my_template = partner_interface_template
                        partner_interface.coord = I_B
                        partner_interface.my_residues = R_B
                        partner_interface.energy = E_B
                        partner_interface.my_template.energy = E_B
                        self.interface_list.append(partner_interface)
                        partner_molecule.interface_list.append(partner_interface)

                        # add the chains pair to self.binding_chains_pairs
                        if chain_id < interface_id:
                            binding_chains_pair = (chain_id, interface_id)
                        else:
                            binding_chains_pair = (interface_id, chain_id)
                        if binding_chains_pair not in self.binding_chains_pairs:
                            self.binding_chains_pairs.append(binding_chains_pair)
                    else:
                        # print(f"Interface {interface_id} already exists for molecule {mol_name}")
                        # print(f"Interface {A} already exists for molecule {B}")
                        pass

        # update the interfaces list of each molecule based on the molecule template
        for group in self.chains_group:
            for i, chain_id in enumerate(group):
                # determin the COM and interfaces of the corresponding molecule template
                molecule_template = [mol_template for mol_template in self.molecules_template_list if mol_template.name == self.chains_map[chain_id]][0]
                molecule_0 = [mol for mol in self.molecule_list if mol.name == group[0]][0]
                com_coord = molecule_0.coord
                interface_coords = [interface_template.coord + com_coord for interface_template in molecule_template.interface_template_list]
                interface_template_ids = [interface_template.name for interface_template in molecule_template.interface_template_list]

                # calculate the R and t for the rigid transformation
                if i == 0:
                    # calculate the normal_point for this molecule
                    molecule = [mol for mol in self.molecule_list if mol.name == chain_id][0]
                    molecule.normal_point = [com_coord.x, com_coord.y, com_coord.z + 1] # normal_point - COM is [0,0,1]
                    # no need to transform the first chain
                    continue
                else:
                    chain1 = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == group[0]][0])]
                    chain2 = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])]
                    R, t = rigid_transform_chains(chain1, chain2)
                    com_coord_transformed = apply_rigid_transform(R, t, np.array([com_coord.x, com_coord.y, com_coord.z]))
                    interface_coords_transformed = []
                    for interface_coord in interface_coords:
                        interface_coord_transformed = apply_rigid_transform(R, t, np.array([interface_coord.x, interface_coord.y, interface_coord.z]))
                        interface_coords_transformed.append(interface_coord_transformed)
                    normal_point_transformed = apply_rigid_transform(R, t, np.array([com_coord.x, com_coord.y, com_coord.z + 1]))
                    # update the COM and interfaces of the molecule
                    molecule = [mol for mol in self.molecule_list if mol.name == chain_id][0]
                    molecule.coord = Coords(com_coord_transformed[0], com_coord_transformed[1], com_coord_transformed[2])
                    for j, interface in enumerate(molecule.interface_list):
                        # find the corresponding interface template
                        interface_template_id = interface.my_template.name
                        for k, intf_template in enumerate(interface_template_ids):
                            if interface_template_id == intf_template:
                                interface.coord = Coords(interface_coords_transformed[k][0], interface_coords_transformed[k][1], interface_coords_transformed[k][2])
                                break
                    molecule.normal_point = [normal_point_transformed[0], normal_point_transformed[1], normal_point_transformed[2]]

        self._update_interface_templates_free_required_list()

        self.binding_chains_pairs.sort()
        self.molecule_list.sort(key=lambda m: m.name)
        self.molecules_template_list.sort(key=lambda mt: mt.name)
        self.interface_list.sort(key=lambda i: i.name)
        self.interface_template_list.sort(key=lambda it: it.name)

        # print("binding chains pairs:")
        # for pair in self.binding_chains_pairs:
        #     print(pair)
        # print("molecule list:")
        # for molecule in self.molecule_list:
        #     print(molecule)
        # print("molecule template list:")
        # for molecule_template in self.molecules_template_list:
        #     print(molecule_template)
        # print("interface list:")
        # for interface in self.interface_list:
        #     print(interface)
        # print("interface template list:")
        # for interface_template in self.interface_template_list:
        #     print(interface_template)

        self._build_reactions()

        if standard_output:
            print("Molecules Template and Reactions Template After Regularization:")
            for molecule_template in self.molecules_template_list:
                print(molecule_template)
            for reaction_template in self.reaction_template_list:
                print(reaction_template)

            print("Molecules and Reactions:")
            for molecule in self.molecule_list:
                print(molecule)
            for reaction in self.reaction_list:
                print(reaction)

        if show_coarse_grained_structure:
            self.plot_regularized_structure()

        if save_pymol_script:
            self.save_regularized_coarse_grained_structure()

        self._generate_model_data()

    def _generate_model_data(self) -> None:
        """Generates molecule types and reactions and saves the model."""
        
        # Step 1: Generate molecule types
        molecule_types = []
        for mol_template in self.molecules_template_list:
            mol_name = mol_template.name
            mol_interfaces = []
            for intf_template in mol_template.interface_template_list:
                iface = MoleculeInterface(name=intf_template.name, coord=intf_template.coord)
                mol_interfaces.append(iface)
            molecule = MoleculeType(name=mol_name, interfaces=mol_interfaces, diffusion_translation=mol_template.diffusion_translation, diffusion_rotation=mol_template.diffusion_rotation)
            molecule_types.append(molecule)

        # Step 2: Generate reactions
        reactions = []
        for reaction_template in self.reaction_template_list:
            brad = getattr(reaction_template, 'binding_radius', 1.0)
            brad = brad / 10.0 # convert to nm
            bind_anlges = getattr(reaction_template, 'binding_angles', None)
            if bind_anlges:
                bind_anlges = tuple(angle for angle in bind_anlges)
            norm1 = getattr(reaction_template, 'norm1', [])
            norm2 = getattr(reaction_template, 'norm2', [])
            norm1 = tuple(n for n in norm1)
            norm2 = tuple(n for n in norm2)
            ka: float = getattr(reaction_template, 'ka', 0.0)
            kb: float = getattr(reaction_template, 'kb', 0.0)
            reaction = ReactionType(name=reaction_template.expression, binding_radius=brad, binding_angles=bind_anlges, norm1=norm1, norm2=norm2, ka=ka, kb=kb)
            reactions.append(reaction)

        # Step 3: Save model data
        self.molecule_types = molecule_types
        self.reactions = reactions
        self.name = self.pdb_id
        model_file = f"{self.save_dir}/{self.pdb_id}_model.json"
        self.save_model(model_file)

        print(f"Model saved successfully to {model_file}")

    def plot_regularized_structure(self):
        """
        Visualizes the molecular structure after regularizing molecules.

        This method plots the center of mass (COM) and interface coordinates for 
        each molecule in the system.
        """
        all_points = []
        chain_ids = []
        for _, mol in enumerate(self.molecule_list):
            chain_ids.append(mol.name)
            com_coord = mol.coord
            interface_coords = [interface.coord for interface in mol.interface_list]
            points = []
            points.append([com_coord.x, com_coord.y, com_coord.z])
            for interface_coord in interface_coords:
                points.append([interface_coord.x, interface_coord.y, interface_coord.z])
            all_points.append(points)
        self.plot_points_3d(all_points, chain_ids)

    def save_regularized_coarse_grained_structure(self, output_cif: str = "regularized_coarse_grained_structure.cif", pymol_script: str = "visualize_regularized_coarse_grained.pml"):
        """
        Saves the regularized coarse-grained molecular structure, including centers of mass (COM) 
        and interface coordinates, to a CIF file and generates a PyMOL script for visualization.

        This method performs the following tasks:
        1. Writes a CIF file containing coarse-grained molecular structure data.
        2. Generates a PyMOL script to visualize the saved structure with color-coded spheres 
        representing COMs and interface points, along with dashed lines connecting them.

        Args:
            output_cif : str, optional
                The filename for the output CIF file. Defaults to "regularized_coarse_grained_structure.cif".
            pymol_script : str, optional
                The filename for the output PyMOL visualization script. Defaults to "visualize_regularized_coarse_grained.pml".
        """
        output_cif = os.path.join(self.save_dir, output_cif)
        pymol_script = os.path.join(self.save_dir, pymol_script)

        with open(output_cif, 'w') as cif_file:
            atom_id = 1

            # Write CIF header
            cif_file.write("# Coarse-grained structure CIF file\n")
            cif_file.write("data_coarse_grained\n")
            cif_file.write("_audit_conform_dict.text 'Regularized coarse-grained model generated by ionerdss'\n")
            cif_file.write("loop_\n")
            cif_file.write("_atom_site.group_PDB\n")
            cif_file.write("_atom_site.id\n")
            cif_file.write("_atom_site.label_atom_id\n")
            cif_file.write("_atom_site.label_comp_id\n")
            cif_file.write("_atom_site.label_asym_id\n")
            cif_file.write("_atom_site.Cartn_x\n")
            cif_file.write("_atom_site.Cartn_y\n")
            cif_file.write("_atom_site.Cartn_z\n")
            cif_file.write("_atom_site.occupancy\n")
            cif_file.write("_atom_site.B_iso_or_equiv\n")
            cif_file.write("_atom_site.type_symbol\n")

            # Loop over regularized molecules
            for mol in self.molecule_list:
                if not mol.coord:
                    continue
                # Write COM (center-of-mass) as 'COM' atom
                cif_file.write(
                    f"ATOM  {atom_id:5d}  COM  MOL {mol.name}  "
                    f"{mol.coord.x:8.3f} {mol.coord.y:8.3f} {mol.coord.z:8.3f}  1.00  0.00  C\n"
                )
                atom_id += 1

                # Write each interface atom
                for intf in mol.interface_list:
                    cif_file.write(
                        f"ATOM  {atom_id:5d}  INT  MOL {mol.name}  "
                        f"{intf.coord.x:8.3f} {intf.coord.y:8.3f} {intf.coord.z:8.3f}  1.00  0.00  O\n"
                    )
                    atom_id += 1

        print(f"Regularized coarse-grained structure saved to {output_cif}.")

        # Generate PyMOL script
        with open(pymol_script, 'w') as pml_file:
            pml_file.write("# PyMOL script to visualize regularized coarse-grained structure\n")
            pml_file.write(f"load {output_cif}, coarse_grained\n")
            pml_file.write("hide everything\n")
            pml_file.write("show spheres, name COM\n")
            pml_file.write("show spheres, name INT\n")
            pml_file.write("set sphere_scale, 1.0\n")
            pml_file.write("color red, name COM\n")
            pml_file.write("color blue, name INT\n")

            atom_index = 1
            for mol in self.molecule_list:
                if not mol.coord:
                    continue
                # Create a pseudoatom in PyMOL for the COM
                pml_file.write(
                    f"pseudoatom com_{mol.name}, pos=[{mol.coord.x:.3f}, {mol.coord.y:.3f}, {mol.coord.z:.3f}], color=red\n"
                )

                # For each interface, create a pseudoatom and connect it to COM
                for j, intf in enumerate(mol.interface_list, start=1):
                    pml_file.write(
                        f"pseudoatom int_{mol.name}_{j}, pos=[{intf.coord.x:.3f}, {intf.coord.y:.3f}, {intf.coord.z:.3f}], color=blue\n"
                    )
                    pml_file.write(
                        f"distance line{atom_index}, com_{mol.name}, int_{mol.name}_{j}\n"
                    )
                    # Use f-strings so {atom_index} is replaced numerically
                    pml_file.write(f"set dash_width, 4, line{atom_index}\n")
                    pml_file.write(f"set dash_gap, 0.5, line{atom_index}\n")
                    atom_index += 1

            pml_file.write("set sphere_transparency, 0.2\n")
            pml_file.write("bg_color white\n")
            pml_file.write("zoom all\n")

            pml_file.write(f"load {self.pdb_file}\n")
            pml_file.write("spectrum chain\n")
            pdb_file = os.path.basename(self.pdb_file).split('.')[0]
            pml_file.write(f"hide everything, {pdb_file}\n")
            pml_file.write(f"show cartoon, {pdb_file}\n")
            pml_file.write(f"set cartoon_transparency, 0.7, {pdb_file}\n")
            save_figure = os.path.join(self.save_dir, "comparison_regularized.png")
            pml_file.write(f"png {save_figure}, 800, 800, 150, 1\n")

        print(f"PyMOL script saved to {pymol_script}.")

    def _update_interface_templates_free_required_list(self):
        """
        Updates the `required_free_list` attribute for each interface template by checking 
        potential steric clashes among binding partners within the same molecule template.

        This method iterates over groups of molecular chains and identifies steric clashes 
        between interface partners. If a steric clash is detected, the corresponding interface 
        templates must remain unbound.

        The procedure follows these steps:
        1. Iterate over molecular chain groups.
        2. Check each interface template's first occurrence within a group.
        3. Compare binding partners for steric clashes with interfaces from previous chains.
        4. If a steric clash is detected, update the `required_free_list` to enforce 
        the constraint.
        """
        for group in self.chains_group:
            for i, chain_id in enumerate(group):
                if i == 0:
                    continue
                # find the molecule in the list
                molecule = [mol for mol in self.molecule_list if mol.name == chain_id][0]
                # loop the interfaces list of the molecule
                for interface in molecule.interface_list:
                    # determine if this interface appears first time
                    interface_id = interface.name
                    interface_template_id = interface.my_template.name
                    first_appearance = True
                    for j in range(i):
                        chain_id_2 = group[j]
                        molecule_2 = [mol for mol in self.molecule_list if mol.name == chain_id_2][0]
                        for interface_2 in molecule_2.interface_list:
                            interface_id_2 = interface_2.name
                            interface_template_id_2 = interface_2.my_template.name
                            if interface_template_id == interface_template_id_2:
                                first_appearance = False
                                break
                    if first_appearance:
                        # check the steric clashes between the partner to this interface and partner to the partners to other interfaces of previous chains; two interfaces belong to different interface tempalte
                        my_partner_chain_id = interface_id
                        my_partner_chain = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == my_partner_chain_id][0])]
                        my_chain = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id][0])]
                        for j in range(i):
                            chain_id_2 = group[j]
                            molecule_2 = [mol for mol in self.molecule_list if mol.name == chain_id_2][0]
                            for interface_2 in molecule_2.interface_list:
                                interface_id_2 = interface_2.name
                                interface_template_id_2 = interface_2.my_template.name
                                if interface_template_id != interface_template_id_2:
                                    another_partner_chain_id = interface_id_2
                                    another_partner_chain = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == another_partner_chain_id][0])]
                                    another_chain = self.all_chains[self.all_chains.index([chain for chain in self.all_chains if chain.id == chain_id_2][0])]
                                    R, t = rigid_transform_chains(my_chain, another_chain)
                                    # rotate the CA atoms of my_partner_chain and check the steric clashes with CA atoms of another_partner_chain
                                    my_partner_chain_CA_coords = []
                                    for residue in my_partner_chain:
                                        if is_aa(residue) and 'CA' in residue:
                                            my_partner_chain_CA_coords.append(residue['CA'].coord)
                                    my_partner_chain_CA_coords_transformed = []
                                    for coord in my_partner_chain_CA_coords:
                                        coord_transformed = apply_rigid_transform(R, t, coord)
                                        my_partner_chain_CA_coords_transformed.append(coord_transformed)
                                    another_partner_chain_CA_coords = []
                                    for residue in another_partner_chain:
                                        if is_aa(residue) and 'CA' in residue:
                                            another_partner_chain_CA_coords.append(residue['CA'].coord)
                                    if check_steric_clashes(np.array(my_partner_chain_CA_coords_transformed), np.array(another_partner_chain_CA_coords)):
                                        molecule_template_id = self.chains_map[chain_id]
                                        molecule_template = [mol_template for mol_template in self.molecules_template_list if mol_template.name == molecule_template_id][0]
                                        interface_template_1 = [interface_template for interface_template in molecule_template.interface_template_list if interface_template.name == interface_template_id][0]
                                        interface_template_2 = [interface_template for interface_template in molecule_template.interface_template_list if interface_template.name == interface_template_id_2][0]
                                        if interface_template_id not in interface_template_2.required_free_list:
                                            interface_template_2.required_free_list.append(interface_template_id)
                                        if interface_template_id_2 not in interface_template_1.required_free_list:
                                            interface_template_1.required_free_list.append(interface_template_id_2)

    def _build_reactions(self):
        """
        Constructs `Reaction` objects for each binding pair in `binding_chains_pairs`, 
        including angle calculations. Also creates `ReactionTemplate` objects if they 
        do not already exist for the corresponding reactants.

        This method:
        1. Iterates over molecular binding pairs.
        2. Retrieves the molecules and interfaces involved.
        3. Computes reaction binding angles using molecular coordinates.
        4. Constructs `Reaction` objects and updates the `reaction_list`.
        5. Ensures `ReactionTemplate` objects exist, creating them if necessary.
        """
        for binding_pair in self.binding_chains_pairs:
            molecule_1: CoarseGrainedMolecule = [mol for mol in self.molecule_list if mol.name == binding_pair[0]][0]
            molecule_2: CoarseGrainedMolecule = [mol for mol in self.molecule_list if mol.name == binding_pair[1]][0]
            interface_1: BindingInterface = [interface for interface in molecule_1.interface_list if interface.name == binding_pair[1]][0]
            interface_2: BindingInterface = [interface for interface in molecule_2.interface_list if interface.name == binding_pair[0]][0]

            # build the reaction
            reaction = Reaction()
            reaction.reactants = []
            reaction.products = []
            reaction.binding_angles = []
            reaction.expression = ""
            reaction.reactants.append((molecule_1, interface_1))
            reaction.reactants.append((molecule_2, interface_2))
            reaction.products.append(f"{molecule_1.name}({interface_1.name}!1).{molecule_2.name}({interface_2.name}!1)")
            reaction.expression = f"{molecule_1.name}({interface_1.name}) + {molecule_2.name}({interface_2.name}) <-> {molecule_1.name}({interface_1.name}!1).{molecule_2.name}({interface_2.name}!1)"
            c1 = np.array([molecule_1.coord.x, molecule_1.coord.y, molecule_1.coord.z])
            c2 = np.array([molecule_2.coord.x, molecule_2.coord.y, molecule_2.coord.z])
            i1 = np.array([interface_1.coord.x, interface_1.coord.y, interface_1.coord.z])
            i2 = np.array([interface_2.coord.x, interface_2.coord.y, interface_2.coord.z])
            n1 = np.array(molecule_1.normal_point)
            n2 = np.array(molecule_2.normal_point)
            theta1, theta2, phi1, phi2, omega, sigma_magnitude = angles(c1, c2, i1, i2, n1, n2)
            if len(molecule_1.my_template.interface_template_list) == 1:
                phi1 = 'nan'
            if len(molecule_2.my_template.interface_template_list) == 1:
                phi2 = 'nan'
            reaction.binding_angles = [theta1, theta2, phi1, phi2, omega]
            reaction.norm1 = [0,0,1]
            reaction.norm2 = [0,0,1]
            reaction.binding_radius = sigma_magnitude

            # calculate the rates
            RT: float = 8.314 * 298 / 1000 # kJ/mol at 298 K
            energy: float = interface_1.energy   # kJ/mol
            C0: float = 0.6022 # unit nm^-3 / M

            reaction.kd = np.exp(energy / RT) # unit M
            reaction.ka = 1200 # unit nm^3/us
            reaction.kb = reaction.kd * reaction.ka * C0 * 1e6 # unit /s
            reaction.energy = energy

            self.reaction_list.append(reaction)
            # print("Reaction:")
            # print(reaction.expression)
            # print("Angles:")
            # print(reaction.binding_angles)
            # print("Sigma:")
            # print(reaction.binding_radius)
            # print("c1:")
            # print(np.array([molecule_1.coord.x, molecule_1.coord.y, molecule_1.coord.z]))
            # print("p1:")
            # print(np.array([interface_1.coord.x, interface_1.coord.y, interface_1.coord.z]))
            # print("c2:")
            # print(np.array([molecule_2.coord.x, molecule_2.coord.y, molecule_2.coord.z]))
            # print("p2:")
            # print(np.array([interface_2.coord.x, interface_2.coord.y, interface_2.coord.z]))

            # build the reaction template if it does not exist
            molecule_1_template_id = self.chains_map[molecule_1.name]
            molecule_2_template_id = self.chains_map[molecule_2.name]
            interface_1_template_id = interface_1.my_template.name
            interface_2_template_id = interface_2.my_template.name
            reactants = []
            if molecule_1_template_id < molecule_2_template_id:
                reactants.append(f"{molecule_1_template_id}({interface_1_template_id})")
                reactants.append(f"{molecule_2_template_id}({interface_2_template_id})")
            elif molecule_1_template_id == molecule_2_template_id:
                if interface_1_template_id < interface_2_template_id:
                    reactants.append(f"{molecule_1_template_id}({interface_1_template_id})")
                    reactants.append(f"{molecule_2_template_id}({interface_2_template_id})")
                else:
                    reactants.append(f"{molecule_2_template_id}({interface_2_template_id})")
                    reactants.append(f"{molecule_1_template_id}({interface_1_template_id})")
            else:
                reactants.append(f"{molecule_2_template_id}({interface_2_template_id})")
                reactants.append(f"{molecule_1_template_id}({interface_1_template_id})")
            existed = False
            for reaction_template in self.reaction_template_list:
                if reaction_template.reactants == reactants:
                    existed = True
                    self.reaction_list[-1].my_template = reaction_template
                    # print("My Reaction Template:")
                    # print(reaction_template.expression)
                    # print("Template Angles:")
                    # print(reaction_template.binding_angles)
                    # print("Template Sigma:")
                    # print(reaction_template.binding_radius)
                    break
            if not existed:
                reaction_template = ReactionTemplate()
                reaction_template.reactants = reactants
                reaction_template.products = []
                reaction_template.products.append(f"{molecule_1_template_id}({interface_1_template_id}!1).{molecule_2_template_id}({interface_2_template_id}!1)")
                # reactants and products do not include the interfaces that need to be free, but the expression does
                free_list_1 = ""
                molecule_template_1 = [mol_template for mol_template in self.molecules_template_list if mol_template.name == molecule_1_template_id][0]
                interface_template_1 = [interface_template for interface_template in molecule_template_1.interface_template_list if interface_template.name == interface_1_template_id][0]
                free_interface_template_list_1 = interface_template_1.required_free_list
                for free_interface in free_interface_template_list_1:
                    free_list_1 += f", {free_interface}"
                tmp_reactant_1 = f"{molecule_1_template_id}({interface_1_template_id}{free_list_1})"

                free_list_2 = ""
                molecule_template_2 = [mol_template for mol_template in self.molecules_template_list if mol_template.name == molecule_2_template_id][0]
                interface_template_2 = [interface_template for interface_template in molecule_template_2.interface_template_list if interface_template.name == interface_2_template_id][0]
                free_interface_template_list_2 = interface_template_2.required_free_list
                for free_interface in free_interface_template_list_2:
                    free_list_2 += f", {free_interface}"
                tmp_reactant_2 = f"{molecule_2_template_id}({interface_2_template_id}{free_list_2})"

                tmp_product = f"{molecule_1_template_id}({interface_1_template_id}!1{free_list_1}).{molecule_2_template_id}({interface_2_template_id}!1{free_list_2})"

                reaction_template.expression = f"{tmp_reactant_1} + {tmp_reactant_2} <-> {tmp_product}"

                reaction_template.binding_angles = reaction.binding_angles
                reaction_template.binding_radius = reaction.binding_radius
                reaction_template.norm1 = reaction.norm1
                reaction_template.norm2 = reaction.norm2

                reaction_template.kd = reaction.kd
                reaction_template.kb = reaction.kb
                reaction_template.ka = reaction.ka
                reaction_template.energy = reaction.energy

                self.reaction_template_list.append(reaction_template)
                self.reaction_list[-1].my_template = reaction_template
                # print("My Reaction Template:")
                # print(reaction_template.expression)
                # print("Template Angles:")
                # print(reaction_template.binding_angles)
                # print("Template Sigma:")
                # print(reaction_template.binding_radius)

    def _calc_angle(self, P, Q, R):
        """
        Calculates the angle at point Q formed by vectors P->Q and R->Q.

        Args:
            P (Coords): The first point.
            Q (Coords): The vertex point where the angle is calculated.
            R (Coords): The third point.

        Returns:
            float: The angle in degrees.
        """
        v1 = [(Q - P).x, (Q - P).y, (Q - P).z]
        v2 = [(R - Q).x, (R - Q).y, (R - Q).z]
        theta = np.degrees(math.acos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
        return theta
    
    def _sig_are_similar(self, sig1, sig2, dist_threshold_intra, dist_threshold_inter, angle_threshold):
        """
        Compares two groups of interface interaction geometry signatures.

        Args:
            sig1 (dict): The first interface signature.
            sig2 (dict): The second interface signature.
            dist_threshold_intra (float): Distance threshold for intra-molecular comparisons.
            dist_threshold_inter (float): Distance threshold for inter-molecular comparisons.
            angle_threshold (float): Angle threshold for comparisons.

        Returns:
            bool: True if the signatures are similar within the given thresholds, False otherwise.
        """
        for key in ("dA", "dB"):
            if abs(sig1[key] - sig2[key]) > dist_threshold_intra:
                return False
        for key in ("dAB",):
            if abs(sig1[key] - sig2[key]) > dist_threshold_inter:
                return False
        for key in ("thetaA", "thetaB"):
            if abs(sig1[key] - sig2[key]) > angle_threshold:
                return False
        return True
    
    def _is_existing_mol_temp(self, mol_temp_name):
        """
        Checks if a molecule template with the given name exists in the molecule template list.

        Args:
            mol_temp_name (str): The name of the molecule template to check.

        Returns:
            tuple: (bool, int or None)
                - True and index if the template exists.
                - False and None otherwise.
        """
        for i, mol_temp in enumerate(self.molecules_template_list):
            if mol_temp.name == mol_temp_name:
                return True, i
        return False, None
    
    def _is_existing_mol(self, mol_name):
        """
        Checks if a molecule with the given name exists in the molecule list.

        Args:
            mol_name (str): The name of the molecule to check.

        Returns:
            tuple: (bool, int or None)
                - True and index if the molecule exists.
                - False and None otherwise.
        """
        for i, mol in enumerate(self.molecule_list):
            if mol.name == mol_name:
                return True, i
        return False, None
    
    def _is_existing_interface(self, interface_name, molecule):
        """
        Checks if an interface with the given name exists in the molecule's interface list.

        Args:
            interface_name (str): The name of the interface to check.
            molecule (CoarseGrainedMolecule): The molecule to check within.

        Returns:
            tuple: (bool, int or None)
                - True and index if the interface exists.
                - False and None otherwise.
        """
        for i, interface in enumerate(molecule.interface_list):
            if interface.name == interface_name:
                return True, i
        return False, None
    
    def _is_existing_sig(self, sig, dist_threshold_intra=2.5, dist_threshold_inter=2.5, angle_threshold=25.0):
        """
        Checks if a given interface signature already exists in the interface_signatures list.

        Args:
            sig (dict): The interface signature to check.
            dist_threshold_intra (float, optional): Distance threshold for intra-molecular comparisons. Defaults to 2.5.
            dist_threshold_inter (float, optional): Distance threshold for inter-molecular comparisons. Defaults to 2.5.
            angle_threshold (float, optional): Angle threshold for comparisons. Defaults to 25.0.

        Returns:
            bool: True if the signature exists, False otherwise.
        """
        for existing_sig in self.interface_signatures:
            if self._sig_are_similar(sig, existing_sig, dist_threshold_intra=dist_threshold_intra, dist_threshold_inter=dist_threshold_inter, angle_threshold=angle_threshold):
                return True
        return False
    
    def _compute_diffusion_constants_nm_us(self, R_nm, T=298.0, eta=1e-3):
        """
        Compute translational and rotational diffusion constants.
        
        Args:
            R_nm (float): radius in nanometers
            T (float): temperature in Kelvin
            eta (float): viscosity in Pa·s (default: water)
        
        Returns:
            tuple: (D_t in nm^2/μs, D_r in rad^2/μs)
        """
        kB = 1.380649e-23  # J/K
        R_m = R_nm * 1e-9  # convert nm to meters

        # Diffusion constants
        D_t_m2_per_s = kB * T / (6 * np.pi * eta * R_m)
        D_r_rad2_per_s = kB * T / (8 * np.pi * eta * R_m**3)

        # Convert units
        D_t_nm2_per_us = D_t_m2_per_s * 1e12  # m²/s → nm²/μs
        D_r_rad2_per_us = D_r_rad2_per_s * 1e-6  # rad²/s → rad²/μs

        return D_t_nm2_per_us, D_r_rad2_per_us
    
    def _get_default_energy_table(self):
        """Returns energy table for residue-residue interactions.

        Reference:
            Miyazawa, S., & Jernigan, R. L. (1996). Residue-residue potentials 
            with a favorable contact pair term and an unfavorable high packing density term,
            for simulation and threading. J Mol Biol, 256(3), 623–644.

        Returns:
            dict: A symmetric dictionary with residue pair tuples as keys and contact energies (in RT units) as values.
        """
        residues = [
            'CYS', 'MET', 'PHE', 'ILE', 'LEU', 'VAL', 'TRP', 'TYR', 'ALA', 'GLY',
            'THR', 'SER', 'ASN', 'GLN', 'ASP', 'GLU', 'HIS', 'ARG', 'LYS', 'PRO'
        ]

        # Extracted from the upper triangle of the table (manually transcribed)
        energy_matrix = [
            [-5.44],
            [-4.99, -5.46],
            [-5.80, -6.56, -7.26],
            [-5.50, -6.02, -6.84, -6.54],
            [-5.83, -6.41, -7.28, -7.04, -7.37],
            [-4.96, -5.32, -6.29, -6.05, -6.48, -5.52],
            [-4.95, -5.55, -6.16, -5.78, -6.14, -5.18, -5.06],
            [-4.16, -4.91, -5.66, -5.25, -5.67, -4.62, -4.66, -4.17],
            [-3.57, -3.94, -4.81, -4.58, -4.91, -4.04, -3.82, -3.36, -2.72],
            [-3.16, -3.39, -4.13, -3.78, -4.16, -3.38, -3.42, -3.01, -2.31, -2.24],
            [-3.11, -3.51, -4.28, -4.03, -4.34, -3.46, -3.22, -3.01, -2.32, -2.08, -2.12],
            [-2.86, -3.03, -4.02, -3.52, -3.92, -3.05, -2.99, -2.78, -2.01, -1.82, -1.96, -1.67],
            [-2.59, -2.95, -3.75, -3.24, -3.74, -2.83, -3.07, -2.76, -1.84, -1.74, -1.88, -1.58, -1.68],
            [-2.85, -3.30, -4.10, -3.67, -4.04, -3.07, -3.11, -2.97, -1.89, -1.66, -1.90, -1.49, -1.71, -1.54],
            [-2.41, -2.57, -3.48, -3.17, -3.40, -2.48, -2.84, -2.76, -1.70, -1.59, -1.80, -1.63, -1.68, -1.46, -1.21],
            [-2.27, -2.89, -3.56, -3.27, -3.59, -2.67, -2.99, -2.79, -1.51, -1.22, -1.74, -1.48, -1.51, -1.42, -1.02, -0.91],
            [-3.60, -3.98, -4.77, -4.14, -4.54, -3.58, -3.98, -3.52, -2.41, -2.15, -2.42, -2.11, -2.08, -1.98, -2.32, -2.15, -3.05],
            [-2.57, -3.12, -3.98, -3.63, -4.03, -3.07, -3.41, -3.16, -1.83, -1.72, -1.90, -1.62, -1.64, -1.80, -2.29, -2.27, -2.16, -1.55],
            [-1.95, -2.48, -3.36, -3.01, -3.37, -2.49, -2.69, -2.60, -1.31, -1.15, -1.31, -1.05, -1.21, -1.29, -1.68, -1.80, -1.35, -0.59, -0.12],
            [-3.07, -3.45, -4.25, -3.76, -4.20, -3.32, -3.73, -3.19, -2.03, -1.87, -1.90, -1.57, -1.53, -1.73, -1.33, -1.26, -2.25, -1.70, -0.97, -1.75]
        ]

        energy_table = {}

        for i, res_i in enumerate(residues):
            for j, res_j in enumerate(residues[:i+1]):
                energy = energy_matrix[i][j] + 2.27  # Adjusted energy value
                energy_table[(res_i, res_j)] = energy
                energy_table[(res_j, res_i)] = energy  # symmetry

        return energy_table

    def identify_homologous_chains(self):
        """
        Identifies homologous chains in the molecular structure and populates `self.chain_map` 
        and `self.chain_groups`. Attempts to parse the header from PDB/CIF files first; 
        if unsuccessful, falls back to sequence alignment.
        """
        if self.pdb_file.endswith('.pdb'):
            self._parse_pdb_header()
        elif self.pdb_file.endswith('.cif'):
            self._parse_cif_header()
        if not self.chains_map:
            self._find_homologous_chains_by_alignment()

    def _parse_pdb_header(self):
        """
        Parses the PDB file header to extract homologous chain information 
        (MOL_ID, CHAIN). Populates `self.chain_map` and `self.chain_groups` 
        based on the identified molecular groups.
        """
        try:
            with open(self.pdb_file, 'r') as file:
                current_mol_id = None
                chains_group = []

                for line in file:
                    if line.startswith("COMPND"):
                        if "MOL_ID:" in line:
                            current_mol_id = line.split(":")[1].strip().split(";")[0]
                        elif "CHAIN:" in line and current_mol_id:
                            chains = line.split(":")[1].strip().split(";")[0].split(",")
                            chains_group.append(chains)

                # Group chains with the same MOL_ID as homologous
                available_NERDSS_mol_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
                for i, chains in enumerate(chains_group):
                    for chain in chains:
                        self.chains_map[chain] = available_NERDSS_mol_ids[i]
                self.chains_group = chains_group

                if self.chains_map:
                    print("Homologous chains identified using PDB header:")
                    print(self.chains_map)

        except Exception as e:
            print(f"Failed to parse PDB header for homologous chains: {str(e)}")
            # Attempt to find homologous chains using sequence alignment
            self._find_homologous_chains_by_alignment()

    def _parse_cif_header(self):
        """
        Parses the CIF file header to extract homologous chain information 
        (entity_id). Populates `self.chain_map` and `self.chain_groups`.
        """
        try:
            with open(self.pdb_file, 'r') as file:
                section_found = False
                section_contents = []
                entity_ids = []
                chains_group = []

                for line in file:
                    if line.startswith("loop_"):
                        next_line = next(file).strip()
                        if next_line.startswith("_entity_poly.entity_id"):
                            section_found = True
                            continue
                    if section_found:
                        # record all the contents between the loop_ and the next loop_
                        if line.startswith("loop_"):
                            break
                        section_contents.append(line.strip())
                # loop through the contents to find the chain and entity_id
                for line in section_contents:
                    # split the line by whitespace
                    line = line.split()
                    # if the first element is a number, it is the entity_id
                    if line[0].isdigit():
                        entity_ids.append(line[0])
                    # if the first element includes , split by , to get the chains
                    elif "," in line[0]:
                        chains_group.append(line[0].split(","))

                # Group chains with the same entity_id as homologous
                available_NERDSS_mol_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
                for i, chains in enumerate(chains_group):
                    for chain in chains:
                        self.chains_map[chain] = available_NERDSS_mol_ids[i]
                self.chains_group = chains_group

                if self.chains_map:
                    print("Homologous chains identified using CIF header:")
                    print(self.chains_map)

        except Exception as e:
            print(f"Failed to parse CIF file for homologous chains: {str(e)}")
            # Attempt to find homologous chains using sequence alignment
            self._find_homologous_chains_by_alignment()

    def _find_homologous_chains_by_alignment(self, seq_identity_threshold: float = 90.0):
        """
        Identifies homologous chains by performing global sequence alignment on 
        amino acid sequences. Chains with identity above the specified threshold 
        are grouped together.

        Args:
            seq_identity_threshold (float, optional): Minimum sequence identity 
                percentage to classify chains as homologous. Defaults to 90.0.
        """
        try:
            similar_chains = []
            chains = list(self.all_atoms_structure.get_chains())
            chain_sequences = {}

            for chain in chains:
                sequence = "".join(seq1(residue.resname) for residue in chain.get_residues() if is_aa(residue))
                chain_sequences[chain.id] = sequence

            # Set up the aligner with the same settings as pairwise2.align.globalxx
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 1.0
            aligner.mismatch_score = 0.0
            aligner.open_gap_score = -1.0
            aligner.extend_gap_score = -0.5

            for i, chain1 in enumerate(chains):
                for j, chain2 in enumerate(chains):
                    if i >= j:
                        continue

                    seq1_chain = chain_sequences[chain1.id]
                    seq2_chain = chain_sequences[chain2.id]

                    # Calculate sequence length for identity calculation
                    max_length = max(len(seq1_chain), len(seq2_chain))
                    if max_length == 0:
                        continue

                    # Get alignment
                    alignment = aligner.align(seq1_chain, seq2_chain)[0]
                    
                    # Calculate identity percentage
                    matches = sum(a == b for a, b in zip(alignment[0], alignment[1]) 
                                if a != '-' and b != '-')
                    identity = (matches / max_length) * 100

                    if identity < seq_identity_threshold:
                        continue
                    similar_chains.append((chain1.id, chain2.id))

            graph = defaultdict(set)
            for chain1, chain2 in similar_chains:
                graph[chain1].add(chain2)
                graph[chain2].add(chain1)

            visited = set()
            groups = []

            def dfs(chain, group):
                visited.add(chain)
                group.add(chain)
                for neighbor in graph[chain]:
                    if neighbor not in visited:
                        dfs(neighbor, group)

            for chain in graph:
                if chain not in visited:
                    group = set()
                    dfs(chain, group)
                    groups.append(list(group))
            self.chains_group = groups
            available_NERDSS_mol_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
            for i, chains in enumerate(groups):
                for chain in chains:
                    self.chains_map[chain] = available_NERDSS_mol_ids[i]

            if self.chains_map:
                print("Homologous chains identified using sequence alignment:")
                print(self.chains_map)

        except Exception as e:
            print(f"Failed to find homologous chains using sequence alignment: {str(e)}")

    def _assign_original_chain_ids(self):
        """
        Assigns original chain IDs as molecular types if no homologous chains are detected. 
        Each chain receives a unique letter from A-Z.
        """
        chains = list(self.all_atoms_structure.get_chains())
        available_NERDSS_mol_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
        for i, chain in enumerate(chains):
            self.chains_map[chain.id] = available_NERDSS_mol_ids[i]
        self.chains_group = [[chain.id] for chain in chains]
        print("Using original chain IDs as molecular types:")
        print(self.chains_map)

# -------------------------------------------------------------------------
# helper Classes
# -------------------------------------------------------------------------

class MoleculeTemplate:
    """
    Represents a molecule type in NERDSS, including the molecule's center of mass (COM) 
    and a list of binding interfaces.

    Attributes:
        name (str): Identifier for the molecule type.
        interface_template_list (list): A list of BindingInterfaceTemplate objects that 
            describe the molecule’s binding sites.
        normal_point (list): Default normal vector direction.
    """

    def __init__(self, name: str):
        """
        Initializes a MoleculeTemplate.

        Args:
            name (str): Name/identifier of the molecule type.
        """
        self.name = name
        self.interface_template_list = []
        self.normal_point = [0,0,1]
        self.diffusion_translation = None
        self.diffusion_rotation = None
        self.radius = None

    def __str__(self):
        interfaces = "\n  ".join(str(it) for it in self.interface_template_list)
        return f"Molecule Template: {self.name}\n  Interfaces:\n  {interfaces}"
    
    def __eq__(self, other):
        if not isinstance(other, MoleculeTemplate):
            return False
        return self.name == other.name

class BindingInterfaceTemplate:
    """
    Represents a binding interface template between molecules.

    Attributes:
        name (str): Identifier of the interface template.
        coord (Coords): Relative coordinates of the interface.
        my_residues (list): Residues forming this interface.
        required_free_list (list): Other interface templates that must remain unbound 
            for this interface to bind.
        signature (dict): Stores interface geometry information.
    """

    def __init__(self, name: str):
        """
        Initializes a BindingInterfaceTemplate.

        Args:
            name (str): Identifier for the interface template.
        """
        self.name = name
        self.coord = None
        self.my_residues = []
        self.required_free_list = [] # The list of interface templates that need to be free to bind to this interface template
        self.signature = {}
        self.energy = None

    def __str__(self):
        residues = ", ".join(self.my_residues)
        required_free = ", ".join(self.required_free_list)
        return (f"Interface Template: {self.name}\n"
                f"  Coordinates: {self.coord}\n"
                f"  Residues: {residues}\n"
                f"  Required Free: {required_free}")
    
    def __eq__(self, other):
        if not isinstance(other, BindingInterfaceTemplate):
            return False
        # TODO: check this
        return self.name == other.name

class CoarseGrainedMolecule:
    """
    Represents a coarse-grained molecule in NERDSS, potentially derived from a PDB chain.

    Attributes:
        name (str): Identifier of the molecule.
        my_template (MoleculeTemplate): Reference to the associated molecule template.
        coord (Coords): Center-of-mass coordinates.
        interface_list (list): List of binding interfaces.
        normal_point (list): Normal vector direction.
    """

    def __init__(self, name: str):
        """
        Initializes a CoarseGrainedMolecule.

        Args:
            name (str): Name/identifier of the molecule.
        """
        self.name = name
        self.my_template = None
        self.coord = None
        self.interface_list: list[BindingInterface] = []
        self.normal_point = None
        self.diffusion_translation = None
        self.diffusion_rotation = None
        self.radius = None

    def __str__(self):
        interfaces = "\n  ".join(str(interface) for interface in self.interface_list)
        return (f"CoarseGrainedMolecule: {self.name}\n"
                f"  Template: {self.my_template}\n"
                f"  Coordinates: {self.coord}\n"
                f"  Interfaces:\n  {interfaces}")
    
    def __repr__(self):
        # Similar to __str__ but more formal for debugging
        return self.name
    
    def __eq__(self, other):
        if not isinstance(other, CoarseGrainedMolecule):
            return False
        return self.name == other.name
    
    def __hash__(self):
        return hash(self.name)

class BindingInterface:
    """
    Represents a binding interface between molecules.

    Attributes:
        name (str): Identifier of the binding interface.
        coord (Coords): Position of the interface.
        my_template (BindingInterfaceTemplate): Reference to the associated interface template.
        my_residues (list): Residues included in the interface.
        signature (dict): Stores interface geometry information.
    """

    def __init__(self, name: str):
        """
        Initializes a BindingInterface.

        Args:
            name (str): Identifier for the binding interface.
        """
        self.name = name
        self.coord = None
        self.my_template = None
        self.my_residues = []
        self.signature = {}
        self.energy: float = None

    def __str__(self):
        return (f"BindingInterface: {self.name}\n"
                f"  Template: {self.my_template}\n"
                f"  Coordinates: {self.coord}\n"
                f"  Residue Count: {len(self.my_residues)}\n"
                f"  Residues: {self.my_residues}")
    
    def __eq__(self, other):
        if not isinstance(other, BindingInterface):
            return False
        return self.my_template == other.my_template

class ReactionTemplate:
    """
    Defines a reaction template between two MoleculeTemplates.

    Attributes:
        expression (str): Textual representation of the reaction.
        reactants (list): List of reactant molecule/interface templates.
        products (list): List of product molecule/interface templates.
        binding_angles (tuple): Tuple describing binding angles (theta1, theta2, phi1, phi2, omega).
        binding_radius (float): Distance between binding interfaces.
        norm1 (list): Normal vector of the first reactant.
        norm2 (list): Normal vector of the second reactant.
    """

    def __init__(self):
        """
        Initializes a ReactionTemplate with default values.
        """
        self.expression = None
        self.reactants = None
        self.products = None
        self.binding_angles = None
        self.binding_radius = None
        self.norm1 = None
        self.norm2 = None
        self.kd = None
        self.ka = None
        self.kb = None
        self.energy = None

    def __str__(self):
        return (f"Reaction Template: {self.expression}\n"
                f"  Reactants: {self.reactants}\n"
                f"  Products: {self.products}\n"
                f"  Binding Angles: {self.binding_angles}\n"
                f"  Binding Radius: {self.binding_radius / 10:.6f} nm\n"
                f"  norm1: {self.norm1}\n"
                f"  norm2: {self.norm2}")
    
    def __eq__(self, other):
        if not isinstance(other, ReactionTemplate):
            return False
        return self.expression == other.expression

class Reaction:
    """
    Represents an actual reaction between two Molecule objects.

    Attributes:
        expression (str): Textual representation of the reaction.
        reactants (list): List of reactant molecules/interfaces.
        products (list): List of product molecules/interfaces.
        binding_angles (tuple): Tuple describing binding angles (theta1, theta2, phi1, phi2, omega).
        binding_radius (float): Distance between binding interfaces.
        norm1 (list): Normal vector of the first reactant.
        norm2 (list): Normal vector of the second reactant.
    """

    def __init__(self):
        """
        Initializes a Reaction with default values.
        """
        self.expression = None
        self.reactants = None
        self.products = None
        self.binding_angles = None
        self.binding_radius = None
        self.norm1 = None
        self.norm2 = None
        self.my_template = None
        self.kd: float = None
        self.ka: float = None
        self.kb: float = None
        self.energy: float = None

    def __str__(self):
        return (f"Reaction: {self.expression}\n"
                f"  Reactants: {self.reactants}\n"
                f"  Products: {self.products}\n"
                f"  Binding Angles: {self.binding_angles}\n"
                f"  Binding Radius: {self.binding_radius / 10:.6f} nm")

    def __repr__(self):
        return f"Reaction({self.expression})"

    def __eq__(self, other):
        if not isinstance(other, Reaction):
            return False
        return self.expression == other.expression
    
    def __hash__(self):
        return hash(self.expression)

# -------------------------------------------------------------------------
# helper functions - geometry transformation
# -------------------------------------------------------------------------

def rigid_transform_3d(points_a: np.ndarray, points_b: np.ndarray):
    """
    Computes a rigid transformation (rotation + translation) that aligns 
    `points_a` to `points_b` using Singular Value Decomposition (SVD).

    Args:
        points_a (np.ndarray): Shape (N, 3), first set of 3D points.
        points_b (np.ndarray): Shape (N, 3), second set of 3D points.

    Returns:
        tuple:
            - np.ndarray: 3x3 rotation matrix `R`
            - np.ndarray: 3-element translation vector `t`
    """
    assert len(points_a) == len(points_b), "Point sets must be same length."
    centroid_a = points_a[0]
    centroid_b = points_b[0]
    pa = points_a[1:] - centroid_a
    pb = points_b[1:] - centroid_b
    h = pa.T @ pb
    u, s, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    t = centroid_b - r @ centroid_a
    return r, t


def apply_rigid_transform(r: np.ndarray, t: np.ndarray, point: np.ndarray):
    """
    Applies a rigid transformation (rotation + translation) to a point.

    Args:
        r (np.ndarray): A 3x3 rotation matrix.
        t (np.ndarray): A 3-element translation vector.
        point (np.ndarray): A shape (3,) array representing a point.

    Returns:
        np.ndarray: Transformed point.
    """
    return (r @ point.T).T + t


def rigid_transform_chains(chain1, chain2):
    """
    Aligns chain1 to chain2 by:
    1. Extracting amino acid sequences.
    2. Performing sequence alignment.
    3. Identifying matching residues.
    4. Computing a coarse-grained set of representative points.
    5. Computing a rigid transformation.

    Args:
        chain1 (Bio.PDB.Chain.Chain): First molecular chain.
        chain2 (Bio.PDB.Chain.Chain): Second molecular chain.

    Returns:
        tuple:
            - np.ndarray: 3x3 rotation matrix `R`
            - np.ndarray: 3-element translation vector `t`
    """

    # Step 1: Extract sequences from both chains
    def extract_sequence(chain):
        """Extracts the amino acid sequence from a chain."""
        return "".join(seq1(residue.resname) for residue in chain.get_residues() if is_aa(residue))

    sequence1 = extract_sequence(chain1)
    sequence2 = extract_sequence(chain2)

    # Step 2: Find the best overlap between the two sequences using PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5
    
    alignments = aligner.align(sequence1, sequence2)
    alignment = alignments[0]  # Get the best alignment

    aligned_seq1 = alignment[0]
    aligned_seq2 = alignment[1]

    # Step 3: Identify matching residue pairs in the aligned sequences
    residue_pairs = []
    idx1, idx2 = 0, 0
    residues1 = [res for res in chain1 if is_aa(res)]
    residues2 = [res for res in chain2 if is_aa(res)]

    for a1, a2 in zip(aligned_seq1, aligned_seq2):
        if a1 == '-' or a2 == '-':
            if a1 != '-':
                idx1 += 1
            if a2 != '-':
                idx2 += 1
            continue
        residue_pairs.append((residues1[idx1]['CA'].coord, residues2[idx2]['CA'].coord))
        idx1 += 1
        idx2 += 1

    # Step 4: Group residues into four spatially groups
    def group_residues(residues, n_groups=4):
        """Groups residues into n_groups based on their spatial proximity."""
        coords = np.array([res for res, _ in residues])
        kmeans = KMeans(n_clusters=n_groups).fit(coords)
        groups = [[] for _ in range(n_groups)]
        for i, label in enumerate(kmeans.labels_):
            groups[label].append(residues[i])
        return groups

    groups = group_residues(residue_pairs)

    # Step 5: Compute the average position of each group and COM
    P = [np.mean([res[0] for res in group], axis=0) for group in groups]
    Q = [np.mean([res[1] for res in group], axis=0) for group in groups]
    P.insert(0, np.mean([res[0] for res in residue_pairs], axis=0))
    Q.insert(0, np.mean([res[1] for res in residue_pairs], axis=0))

    P = np.array(P)
    Q = np.array(Q)

    # Step 6: Apply rigid transformation
    R, t = rigid_transform_3d(P, Q)

    return R, t


def check_steric_clashes(points_1, points_2, cutoff: float = 3.5, number_threshold: int = 2):
    """
    Detects steric clashes between two sets of molecular points.

    Args:
        points_1 (np.ndarray): N x 3 coordinates for the first molecule.
        points_2 (np.ndarray): M x 3 coordinates for the second molecule.
        cutoff (float, optional): Distance threshold (default: 3.5 Å).
        threshold (int, optional): Minimum number of close contacts to flag a clash (default: 2).

    Returns:
        bool: True if a steric clash is detected, False otherwise.
    """
    tree = KDTree(points_2)
    clashes = tree.query_ball_point(points_1, r=cutoff)
    return any(len(clash) >= number_threshold for clash in clashes)


# -------------------------------------------------------------------------
# binding angles calculation functions
# -------------------------------------------------------------------------

def unit(x:np.ndarray, eps=10**-6) -> np.ndarray:
    """
    Normalizes a vector to unit length, handling numerical precision errors.

    Args:
        x (np.ndarray): Input vector.
        eps (float, optional): Small threshold for numerical stability.

    Returns:
        np.ndarray: Unit vector."
    """
    x_norm = np.linalg.norm(x)
    if abs(x_norm-1) < eps:
        return x
    elif x_norm < eps:
        return np.zeros(3)
    else:
        return x/x_norm

def _clip_cosine_value(x: float, eps=10**-6) -> float:
    """
    Ensures cosine values remain in the range [-1, 1] for numerical stability.

    Args:
        x (float): Input value.
        eps (float, optional): Small numerical threshold.

    Returns:
        float: Corrected value within [-1, 1].
    """
    if x < -1 and abs(x+1) < eps:
        return -1
    elif x > 1 and abs(x-1) < eps:
        return 1
    elif -1 <= x <= 1:
        return x
    else:
        raise ValueError(f'{x} is out of the range of sin/cos')

def calculate_phi(v:np.ndarray, n:np.ndarray, sigma:np.ndarray, eps=10**-6) -> float:
    """
    Computes the phi angle given three vectors.

    Args:
        v (np.ndarray): Direction vector.
        n (np.ndarray): Normal vector.
        sigma (np.ndarray): Sigma direction.

    Returns:
        float: Computed phi angle.
    """

    # calculate phi
    t1 = unit(np.cross(v, sigma))
    t2 = unit(np.cross(v, n))
    phi = math.acos(_clip_cosine_value(np.dot(t1, t2)))

    # determine the sign of phi (+/-)
    v_uni = unit(v)
    n_proj = n - v_uni * np.dot(v_uni, n)
    sigma_proj = sigma - v_uni * np.dot(v_uni, sigma)
    phi_dir = unit(np.cross(sigma_proj, n_proj))

    if np.dot(v_uni, phi_dir) > 0:
        phi = -phi
    else:
        phi = phi
    
    return phi

def angles(com1, com2, int_site1, int_site2, normal_point1, normal_point2, eps=10**-6):
    """
    Computes binding angles for two molecules based on their center-of-mass (COM),
    interface sites, and normal vectors.

    This function determines five angles (theta1, theta2, phi1, phi2, omega) that describe
    the relative orientation between two interacting molecules.

    Args:
        com1 (array-like): Coordinates of the center-of-mass of molecule 1.
        com2 (array-like): Coordinates of the center-of-mass of molecule 2.
        int_site1 (array-like): Coordinates of the interface site on molecule 1.
        int_site2 (array-like): Coordinates of the interface site on molecule 2.
        normal_point1 (array-like): A point defining the orientation of molecule 1.
        normal_point2 (array-like): A point defining the orientation of molecule 2.
        eps (float, optional): Small numerical threshold to prevent division errors. Default is 1e-6.

    Returns:
        tuple:
            - str: Theta1 (binding angle in radians, formatted as string).
            - str: Theta2 (binding angle in radians, formatted as string).
            - str: Phi1 (torsion angle in radians, formatted as string, or 'nan' if undefined).
            - str: Phi2 (torsion angle in radians, formatted as string, or 'nan' if undefined).
            - str: Omega (twist angle in radians, formatted as string).
            - float: Distance between interface sites (sigma magnitude).
    """

    # Convert sequences into arrays for convinience
    com1 = np.array(com1)
    com2 = np.array(com2)
    int_site1 = np.array(int_site1)
    int_site2 = np.array(int_site2)
    normal_point1 = np.array(normal_point1)
    normal_point2 = np.array(normal_point2)

    # Get Vectors
    v1 = int_site1 - com1 # from COM to interface (particle 1)
    v2 = int_site2 - com2  # from COM to interface (particle 2)
    sigma1 = int_site1 - int_site2 # sigma, from p2 to p1
    sigma2 = int_site2 - int_site1  # sigma, from p1 to p2
    n1 = unit(normal_point1 - com1) # normal vector for p1
    n2 = unit(normal_point2 - com2) # normal vector for p2

    # Calculate the magnititude of sigma
    sigma_magnitude = np.linalg.norm(sigma1)

    # Calculate theta1 and theta2
    costheta1 = np.dot(v1, sigma1) / np.linalg.norm(v1) / np.linalg.norm(sigma1)
    costheta2 = np.dot(v2, sigma2) / np.linalg.norm(v2) / np.linalg.norm(sigma2)
    theta1 = math.acos(_clip_cosine_value(costheta1))
    theta2 = math.acos(_clip_cosine_value(costheta2))

    # check geometry
    errormsg = ''
    iferror = False # determine if v // n
    if np.linalg.norm(np.cross(n1, v1)) < eps:
        iferror = True
        errormsg += f'\n\tn1 ({n1}) and v1 ({v1}) parallel, phi1 not available'
    if np.linalg.norm(np.cross(n2, v2)) < eps:
        iferror = True
        errormsg += f'\n\tn2 ({n2}) and v2 ({v2}) parallel, phi2 not available'
    if iferror:
        raise ValueError(errormsg)

    # determine if phi1 exists (v1 // sigma1 ?)
    if np.linalg.norm(np.cross(sigma1, v1)) < eps:
        phi1 = float('nan')
        # omega_parallel = True
        omega_t1 = unit(np.cross(sigma1, n1))
    else:
        phi1 = calculate_phi(v1, n1, sigma1)
        omega_t1 = unit(np.cross(sigma1, v1))

    # determine if phi2 exists (v2 // sigma2 ?)
    if np.linalg.norm(np.cross(sigma2, v2)) < eps:
        phi2 = float('nan')
        # omega_parallel = True
        omega_t2 = unit(np.cross(sigma1, n2))
    else:
        phi2 = calculate_phi(v2, n2, sigma2)
        omega_t2 = unit(np.cross(sigma1, v2))

    # calculate omega (both cases are same)
    omega = math.acos(_clip_cosine_value(np.dot(omega_t1, omega_t2)))
    # determine the sign of omega (+/-)
    sigma1_uni = unit(sigma1)
    sigma1xomega_t1 = np.cross(sigma1, omega_t1)
    sigma1xomega_t2 = np.cross(sigma1, omega_t2)
    omega_dir = unit(np.cross(sigma1xomega_t1, sigma1xomega_t2))
    if np.dot(sigma1_uni, omega_dir) > 0:
        omega = -omega
    else:
        omega = omega

    if abs(theta1 - np.pi) < eps:
        theta1 = 'M_PI'
    else:
        theta1 = "%.6f" % theta1
    if abs(theta2 - np.pi) < eps:
        theta2 = 'M_PI'
    else:
        theta2 = "%.6f" % theta2
    if abs(phi1 - np.pi) < eps:
        phi1 = 'M_PI'
    else:
        phi1 = "%.6f" % phi1
    if abs(phi2 - np.pi) < eps:
        phi2 = 'M_PI'
    else:
        phi2 = "%.6f" % phi2
    if abs(omega - np.pi) < eps:
        omega = 'M_PI'
    else:
        omega = "%.6f" % omega

    return theta1, theta2, phi1, phi2, omega, sigma_magnitude

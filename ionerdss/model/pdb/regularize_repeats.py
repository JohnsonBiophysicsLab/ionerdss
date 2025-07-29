"""
regularize_repeats.py

Provides utilities to reindex and relabel coarse-grained chain/interface data into a 
canonical, compact format for further modeling and export.

In the NERDSS pipeline, the coarse-graining step returns raw data:
- Chain IDs are strings (e.g., 'A', 'B', 'E')
- Interfaces are lists of neighbors using original IDs
- Interface order may vary between repeated chains

This module:
- Reindexes chains to `0, 1, 2, ..., N-1`
- Reindexes interfaces per chain to `[0, 1, ..., M-1]` with consistent ordering
- Ensures symmetry-equivalent chains have consistent interface maps
- Optionally merges chain/interface data by repeated_chainy group

Functions
---------
- regularize_model(raw_data, chains_map=None)
    Canonicalize and simplify the structure of the coarse-grained model.
"""

import numpy as np
from Bio.PDB.Polypeptide import is_aa

from ionerdss.model.components import (
    MoleculeTemplate,
    MoleculeType,
    MoleculeInterface,
    CoarseGrainedMolecule,
    BindingInterface,
    BindingInterfaceTemplate,
)
from ionerdss.model.pdb.repeated_interface_signature import (
    signature_are_similar,
    signature_difference,
    build_signature,
    invert_signature,
    signature_hash,
)
from ionerdss.model.pdb.geometry import rigid_transform_chains
from ionerdss.utils.coords import Coords
from ionerdss.utils.diffusion_constant import compute_diffusion_constants_nm_us
from ionerdss.utils.rigid_transform import apply_rigid_transform
from ionerdss.utils.steric_clash import check_clashes_between_two_sets

def process_chain_group(group, chains_map, molecule_template_list, molecule_list, cg_model,
                        dist_threshold_intra,
                        dist_threshold_inter,
                        angle_threshold,
                        signature_to_template_map=None):
    """
    Processes a group of molecular chains and assigns molecule and interface templates.

    For each chain in the group, it creates or reuses molecule templates and molecules,
    then constructs and registers their interfaces. Interface signatures are compared
    and stored to detect repeated patterns.

    Args:
        group (list of str): List of chain IDs in the group.
        chains_map (dict): Mapping of chain ID to molecule template name.
        molecule_template_list (list): Existing list of MoleculeTemplate objects.
        molecule_list (list): Existing list of CoarseGrainedMolecule objects.
        cg_model (dict): Coarse-grained structure data.
        dist_threshold_intra (float): Intra-molecular distance threshold.
        dist_threshold_inter (float): Inter-molecular distance threshold.
        angle_threshold (float): Angular deviation threshold.

    Returns:
        tuple: (interface_signatures, interface_template_list, interface_list, binding_chains_pairs)
    """
    # Initialize local collections and thresholds for this group
    if signature_to_template_map is None:
        signature_to_template_map = {}
    interface_signatures = []
    interface_template_list = []
    interface_list = []
    binding_chains_pairs = []

    molecule_template = get_or_create_molecule_template(chains_map[group[0]],
                                                        molecule_template_list)

    for j, chain_id in enumerate(group):
        molecule = get_or_create_molecule(
            chain_id, molecule_template, molecule_list, cg_model
        )

        process_interfaces_for_chain(
            chain_id, group, j, molecule, cg_model, chains_map,
            molecule_list, molecule_template_list, interface_signatures,
            dist_threshold_intra, dist_threshold_inter, angle_threshold,
            interface_template_list, interface_list, binding_chains_pairs,
            signature_to_template_map
        )

    return interface_signatures, interface_template_list, interface_list, binding_chains_pairs

def get_or_create_molecule_template(name, template_list):
    """
    Retrieves a molecule template by name or creates a new one if it doesn't exist.

    Args:
        name (str): Name of the molecule template.
        template_list (list): List of MoleculeTemplate objects.

    Returns:
        MoleculeTemplate: The found or newly created molecule template.
    """
    exists, idx = _is_existing_mol_temp(name, template_list)
    if exists:
        return template_list[idx]
    new_template = MoleculeTemplate(name)
    template_list.append(new_template)
    return new_template

def get_or_create_molecule(name, template, molecule_list, cg_model):
    """
    Retrieves a molecule by name or creates a new CoarseGrainedMolecule if it doesn't exist.

    Initializes its geometry, radius, and diffusion properties based on the CG model.

    Args:
        name (str): Molecule name (usually chain ID).
        template (MoleculeTemplate): Template associated with the molecule.
        molecule_list (list): List of CoarseGrainedMolecule instances.
        cg_model (dict): Structure data including COMs and radii.

    Returns:
        CoarseGrainedMolecule: The found or newly created molecule object.
    """
    exists, idx = _is_existing_mol(name, molecule_list)
    if exists:
        molecule = molecule_list[idx]
    else:
        molecule = CoarseGrainedMolecule(name)
        molecule.coord = get_chain_data(name, cg_model, "COMs")
        molecule.radius = get_chain_data(name, cg_model, "radii")
        molecule_list.append(molecule)
        template.radius = molecule.radius
    molecule.my_template = template
    molecule.diffusion_translation, molecule.diffusion_rotation = compute_diffusion_constants_nm_us(molecule.radius / 10.0)
    molecule.my_template.diffusion_translation = molecule.diffusion_translation
    molecule.my_template.diffusion_rotation = molecule.diffusion_rotation
    return molecule

def get_chain_data(chain_id, cg_model, key):
    """
    Retrieves chain-specific data from the coarse-grained model using a key.

    Args:
        chain_id (str): Chain identifier.
        cg_model (dict): Coarse-grained structural data.
        key (str): Key to extract data (e.g., "COMs", "interface_coords").

    Returns:
        Any: The value corresponding to the key and chain ID.
    """
    return cg_model[key][
        cg_model["chains"].index(
            [chain for chain in cg_model["chains"] if chain.id == chain_id][0]
        )
    ]

def process_interfaces_for_chain(
    chain_id, group, j, molecule, cg_model, chains_map,
    molecule_list, molecule_template_list, interface_signatures,
    dist_threshold_intra, dist_threshold_inter, angle_threshold,
    interface_template_list, interface_list, binding_chains_pairs,
    signature_to_template_map):
    """
    Processes all interfaces for a given chain within its group.

    Creates or reuses interface templates, computes signature geometry,
    and registers interfaces between molecule pairs. Detects repeated
    interfaces using distance and angle thresholds.

    Args:
        chain_id (str): Chain ID being processed.
        group (list): List of chain IDs in the group.
        j (int): Index of the current chain within the group.
        molecule (CoarseGrainedMolecule): The molecule corresponding to this chain.
        cg_model (dict): Coarse-grained model data.
        chains_map (dict): Chain-to-molecule mapping.
        molecule_list (list): All molecules created so far.
        molecule_template_list (list): All molecule templates.
        interface_signatures (list): List of known interface signatures.
        dist_threshold_intra (float): Intra-chain distance threshold.
        dist_threshold_inter (float): Inter-chain distance threshold.
        angle_threshold (float): Angular tolerance for matching.
        interface_template_list (list): Interface template objects.
        interface_list (list): Interface instances to populate.
        binding_chains_pairs (list): List of all interface partner chain pairs.

    Returns:
        None
    """
    molecule_name = chain_id
    for i, interface_id in enumerate(get_chain_data(molecule_name, cg_model, "interfaces")):
        partner_name = interface_id
        partner_template_name = chains_map[partner_name]
        partner_template = get_or_create_molecule_template(partner_template_name, molecule_template_list)
        partner_molecule = get_or_create_molecule(partner_name, partner_template, molecule_list, cg_model)

        COM_A = get_chain_data(molecule_name, cg_model, "COMs")
        I_A = get_chain_data(molecule_name, cg_model, "interface_coords")[i]
        COM_B = get_chain_data(partner_name, cg_model, "COMs")

        I_B, R_B, E_B = None, None, None
        for k, pid in enumerate(get_chain_data(partner_name, cg_model, "interfaces")):
            if pid == molecule_name:
                I_B = get_chain_data(partner_name, cg_model, "interface_coords")[k]
                R_B = get_chain_data(partner_name, cg_model, "interface_residues")[k]
                E_B = get_chain_data(partner_name, cg_model, "interface_energies")[k]
                break

        signature = build_signature(COM_A, I_A, COM_B, I_B)
        inverted_signature = invert_signature(signature)
        sig_key = signature_hash(signature)
        inv_key = signature_hash(inverted_signature)

        if sig_key in signature_to_template_map:
            interface_template = signature_to_template_map[sig_key]
            partner_template_template = signature_to_template_map[inv_key]
            register_interfaces(
                chain_id, interface_id, molecule, partner_molecule,
                interface_template, partner_template_template,
                cg_model, I_B, R_B, E_B, i,
                interface_list, binding_chains_pairs
            )
        else:
            is_homodimerize = determine_homo_dimerization(signature, molecule_name, partner_name, chains_map, dist_threshold_intra, angle_threshold)
            interface_templates = build_new_interface_templates(
                is_homodimerize, j, i, group, chain_id, molecule_name, partner_name, cg_model,
                signature, chains_map, molecule_template_list, interface_template_list
            )

            # Register templates in the map
            if is_homodimerize:
                signature_to_template_map[sig_key] = interface_templates[0]
                signature_to_template_map[inv_key] = interface_templates[0]
            else:
                signature_to_template_map[sig_key] = interface_templates[0]
                signature_to_template_map[inv_key] = interface_templates[1]

            for interface_template in interface_templates:
                register_interfaces(
                    chain_id, interface_id, molecule, partner_molecule,
                    interface_template,
                    partner_template if interface_template != partner_template else interface_template,
                    cg_model, I_B, R_B, E_B, i,
                    interface_list, binding_chains_pairs
                )

def determine_homo_dimerization(signature, mol_a, mol_b,
                                chains_map, dist_thresh, angle_thresh):
    """
    Determines whether the interaction between two molecules is a homo-dimerization
    based on their mapped template identity and symmetry of interaction geometry.

    Args:
        signature (dict): Geometric signature with 'dA', 'dB', 'thetaA', 'thetaB'.
        mol_A (str): Chain ID of molecule A.
        mol_B (str): Chain ID of molecule B.
        chains_map (dict): Mapping from chain ID to molecule template name.
        dist_thresh (float): Distance threshold for symmetry.
        angle_thresh (float): Angular threshold for symmetry.

    Returns:
        bool: True if the interaction is symmetric and chains share the same template.
    """
    if chains_map[mol_a] != chains_map[mol_b]:
        return False
    return (
        abs(signature["dA"] - signature["dB"]) < dist_thresh and
        abs(signature["thetaA"] - signature["thetaB"]) < angle_thresh
    )

def build_new_interface_templates(
    is_homo, j, i, group, chain_id, mol_A, mol_B, cg_model,
    signature, chains_map, mol_template_list, interface_template_list
):
    """
    Constructs one or two BindingInterfaceTemplate(s) depending on homo- or hetero-dimer.

    Args:
        is_homo (bool): Whether the interaction is symmetric (homo-dimer).
        j (int): Index of current chain.
        i (int): Index of interface.
        group (list): Group of chain IDs.
        chain_id (str): Current chain ID.
        mol_A (str): Chain ID of molecule A.
        mol_B (str): Chain ID of molecule B.
        cg_model (dict): CG model with COMs and interface coords.
        signature (dict): Geometric signature.
        chains_map (dict): Chain-to-template map.
        mol_template_list (list): Existing molecule templates.
        interface_template_list (list): Output list to append new templates.

    Returns:
        list of BindingInterfaceTemplate: One or two templates
    """
    def get_interface_template_id(partner_prefix):
        suffix = sum(1 for mt in mol_template_list for it in mt.interface_template_list if it.name.startswith(partner_prefix)) + 1
        return f"{partner_prefix}_{suffix}"

    if is_homo:
        # Interface is on mol_A, binding to mol_B (same template)
        partner_prefix = chains_map[mol_B]  # or mol_A
        interface_template_id = get_interface_template_id(partner_prefix)
        interface_template = BindingInterfaceTemplate(interface_template_id)
        interface_template.signature = signature
        offset = compute_interface_offset(j, i, group, chain_id, cg_model)
        interface_template.coord = offset
        print(f"[DEBUG] {interface_template_id} offset = {offset}")

        # Assign to correct molecule template (chain_id is mol_A)
        my_template = get_or_create_molecule_template(chains_map[chain_id], mol_template_list)
        my_template.interface_template_list.append(interface_template)
        interface_template_list.append(interface_template)
        return [interface_template]

    else:
        templates = []
        for mol, sig, side_chain_id in [
            (mol_A, signature, chain_id),
            (mol_B, invert_signature(signature), mol_B if chain_id == mol_A else mol_A)
        ]:
            partner_prefix = chains_map[mol_B if mol == mol_A else mol_A]
            interface_template_id = get_interface_template_id(partner_prefix)
            interface_template = BindingInterfaceTemplate(interface_template_id)
            interface_template.signature = sig
            offset = compute_interface_offset(j, i, group, side_chain_id, cg_model)
            interface_template.coord = offset
            print(f"[DEBUG] {interface_template_id} offset = {offset}")

            # Attach to molecule template that interface lives on
            my_template = get_or_create_molecule_template(chains_map[side_chain_id], mol_template_list)
            my_template.interface_template_list.append(interface_template)
            interface_template_list.append(interface_template)
            templates.append(interface_template)

        return templates


def compute_interface_offset(j, i, group, chain_id, cg_model):
    """
    Computes the relative coordinate of the interface site with respect to COM.
    Applies rigid transform if not the first chain in the group.

    Args:
        j (int): Chain index.
        i (int): Interface index.
        group (list): Group of chains.
        chain_id (str): Chain ID of current molecule.
        cg_model (dict): CG data containing coords.

    Returns:
        Coords: Relative offset from COM to interface site.
    """
    if j == 0:
        pt = get_chain_data(chain_id, cg_model, "interface_coords")[i]
        com = get_chain_data(chain_id, cg_model, "COMs")
        offset = np.array([pt.x, pt.y, pt.z]) - np.array([com.x, com.y, com.z])
        return Coords(*offset)
    else:
        chain1 = [c for c in cg_model["chains"] if c.id == group[0]][0]
        chain2 = [c for c in cg_model["chains"] if c.id == chain_id][0]
        R, t = rigid_transform_chains(chain2, chain1)
        Q_COM = get_chain_data(chain_id, cg_model, "COMs")
        Q_pt = get_chain_data(chain_id, cg_model, "interface_coords")[i]
        Q2 = [apply_rigid_transform(R, t, pt.to_numpy()) for pt in [Q_COM, Q_pt]]
        assert isinstance(Q_COM, Coords), f"Expected Coords, got {type(Q_COM)}"
        return Coords(*(Q2[1] - Q2[0]))

def match_existing_templates(signature, signature_conj, mol_template_list, dist_intra, dist_inter, angle):
    """
    Finds the closest matching existing interface templates for a given pair of signatures.

    Args:
        signature (dict): Signature from A to B.
        signature_conj (dict): Conjugated signature from B to A.
        mol_template_list (list): List of molecule templates.
        dist_intra (float): Intra distance threshold.
        dist_inter (float): Inter distance threshold.
        angle (float): Angle threshold.

    Returns:
        tuple: (interface_template_A, interface_template_B)
    """
    def find_match(sig):
        matches = [
            (it, mt) for mt in mol_template_list for it in mt.interface_template_list
            if signature_are_similar(sig, it.signature, dist_intra, dist_inter, angle)
        ]
        if not matches:
            raise ValueError(f"No matching interface template for signature: {sig}")
        return sorted(matches, key=lambda pair: signature_difference(sig, pair[0].signature))[0]

    return find_match(signature), find_match(signature_conj)

def register_interfaces(
    chain_id, interface_id, molecule, partner_molecule, interface_template,
    partner_template, cg_model, I_B, R_B, E_B, i,
    interface_list, binding_chains_pairs
):
    """
    Creates and registers BindingInterface objects and links them to molecules.

    Args:
        chain_id (str): Chain ID of current molecule.
        interface_id (str): Partner chain ID.
        molecule (CoarseGrainedMolecule): Molecule object.
        partner_molecule (CoarseGrainedMolecule): Partner molecule.
        interface_template (BindingInterfaceTemplate): Template for this interface.
        partner_template (BindingInterfaceTemplate): Template for partner.
        cg_model (dict): Structure data.
        I_B, R_B, E_B: Coordinates, residues, and energy for partner.
        i (int): Interface index.
        interface_list (list): Global list to register new interfaces.
        binding_chains_pairs (list): List of all chain-chain binding pairs.
    """
    if not _is_existing_interface(interface_id, interface_list)[0]:
        interface = BindingInterface(interface_id)
        interface.my_template = interface_template
        interface.coord = get_chain_data(chain_id, cg_model, "interface_coords")[i]
        interface.my_residues = get_chain_data(chain_id, cg_model, "interface_residues")[i]
        interface.energy = get_chain_data(chain_id, cg_model, "interface_energies")[i]
        interface.my_template.energy = interface.energy
        interface_list.append(interface)
        molecule.interface_list.append(interface)

        partner_interface = BindingInterface(chain_id)
        partner_interface.my_template = partner_template
        partner_interface.coord = I_B
        partner_interface.my_residues = R_B
        partner_interface.energy = E_B
        partner_interface.my_template.energy = E_B
        interface_list.append(partner_interface)
        partner_molecule.interface_list.append(partner_interface)

        pair = tuple(sorted((chain_id, interface_id)))
        if pair not in binding_chains_pairs:
            binding_chains_pairs.append(pair)

def update_molecule_interfaces(
    chains_groups, chains_map,
    mol_template_list, molecule_list,
    cg_model
):
    """
    Aligns and updates the interface coordinates for all chains in a group so that
    their geometry is consistent with the first chain in the group.

    For homologous chains, applies a rigid body transformation to align interface
    coordinates and update COM and normal vectors.

    Args:
        chains_groups (list): List of chain ID groups.
        chains_map (dict): Mapping from chain IDs to molecule template names.
        mol_template_list (list): List of molecule templates.
        molecule_list (list): List of CoarseGrainedMolecule objects.
        cg_model (dict): Coarse-grained structural data.
    """
    for group in chains_groups:
        for i, chain_id in enumerate(group):
            mol_template = next(mt for mt in mol_template_list if mt.name == chains_map[chain_id])
            mol0 = next(m for m in molecule_list if m.name == group[0])
            com_coord = mol0.coord
            interface_coords = [it.coord + com_coord for it in mol_template.interface_template_list]
            interface_ids = [it.name for it in mol_template.interface_template_list]

            if i == 0:
                mol = next(m for m in molecule_list if m.name == chain_id)
                mol.normal_point = [com_coord.x, com_coord.y, com_coord.z + 1]
                continue

            chain1 = [c for c in cg_model["chains"] if c.id == group[0]][0]
            chain2 = [c for c in cg_model["chains"] if c.id == chain_id][0]
            R, t = rigid_transform_chains(chain1, chain2)
            com_trans = apply_rigid_transform(R, t, np.array([com_coord.x, com_coord.y, com_coord.z]))
            intf_coords_trans = [apply_rigid_transform(R, t, np.array([pt.x, pt.y, pt.z])) for pt in interface_coords]
            norm_pt = apply_rigid_transform(R, t, np.array([com_coord.x, com_coord.y, com_coord.z + 1]))

            mol = next(m for m in molecule_list if m.name == chain_id)
            mol.coord = Coords(*com_trans)
            mol.normal_point = list(norm_pt)
            for j, interface in enumerate(mol.interface_list):
                for k, tid in enumerate(interface_ids):
                    if interface.my_template.name == tid:
                        interface.coord = Coords(*intf_coords_trans[k])
                        break


def regularize_repeated_chains(cg_model, chains_map, chains_groups,
                               dist_threshold_intra=3.5,
                               dist_threshold_inter=3.5,
                               angle_threshold=25.0,
                               standard_output=False):
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
    if not chains_groups:
        raise ValueError("Invalid chains group!")
    for group in chains_groups:
        group.sort()
    chains_groups.sort()

    # ... [rest of the regularization logic remains mostly unchanged] ...    

    # check if the structure has homologous chains
    # if any element in chains_group has more than one chain, then it has homologous chains
    has_homologous_chains = any(len(group) > 1 for group in chains_groups)
    if not has_homologous_chains:
        dist_threshold_intra = 0.0
        dist_threshold_inter = 0.0
        angle_threshold = 0.0

    molecule_list = []
    molecule_template_list = []
    interface_template_list = []       # This is the one you return
    interface_list = []                # Ditto

    # Collect results from all groups
    all_interface_signatures = []
    all_binding_pairs = []
    signature_to_template_map = {}

    for group in chains_groups:
        sigs, templates, interfaces, pairs = process_chain_group(
            group, chains_map, molecule_template_list, molecule_list, cg_model,
            dist_threshold_intra, dist_threshold_inter,
            angle_threshold, signature_to_template_map
        )
        all_interface_signatures.extend(sigs)
        interface_template_list.extend(templates)
        interface_list.extend(interfaces)
        all_binding_pairs.extend(pairs)

    # Post-processing:
    update_molecule_interfaces(
        chains_groups, chains_map,
        molecule_template_list, molecule_list,
        cg_model
    )

    _update_interface_templates_free_required_list(chains_groups,
                                                   molecule_list,
                                                   cg_model["chains"],
                                                   chains_map, molecule_template_list)

    all_binding_pairs.sort()
    molecule_list.sort(key=lambda m: m.name)
    molecule_template_list.sort(key=lambda mt: mt.name)
    interface_list.sort(key=lambda i: i.name)
    interface_template_list.sort(key=lambda it: it.name)

    if standard_output:
        print("Molecules Template and Reactions Template After Regularization:")
        for molecule_template in molecule_template_list:
            print(molecule_template)

        print("Molecules:")
        for molecule in molecule_list:
            print(molecule)

    molecule_types = _generate_molecule_types(molecule_template_list)

   # Return structured model data
    return {
        "molecule_templates": molecule_template_list,
        "molecules": molecule_list,
        "interface_templates": interface_template_list,
        "interfaces": interface_list,
        "binding_pairs": all_binding_pairs,
        "molecule_types": molecule_types,
    }

# This is a modularized scaffold for the original parsing loop
# The functions are broken down by responsibility

# ... (existing code) ...

def _update_interface_templates_free_required_list(
    chains_group, molecule_list, all_chains, chains_map, molecule_template_list
):
    """
    Updates the `required_free_list` attribute for each interface template by checking 
    potential steric clashes among binding partners within the same molecule template.
    """
    print("\n[DEBUG] Entering _update_interface_templates_free_required_list")
    print("---------------------------------------------------------------")
    print(f"Number of chain groups: {len(chains_group)}")
    for idx, group in enumerate(chains_group):
        print(f"  Group {idx}: {group}")
    print(f"Total molecules: {len(molecule_list)}")
    print("  Molecule names: ", [m.name for m in molecule_list])
    print(f"Total chains in cg_model: {len(all_chains)}")
    print("  All chain IDs: ", [chain.id for chain in all_chains])
    print(f"Total molecule templates: {len(molecule_template_list)}")
    for mt in molecule_template_list:
        print(f"  MoleculeTemplate '{mt.name}' has interface templates: {[it.name for it in mt.interface_template_list]}")
    print("Chains map (chain_id -> molecule_template):")
    for k, v in chains_map.items():
        print(f"  {k} -> {v}")
    print("---------------------------------------------------------------\n")
    
    def extract_prefix(name):
        """Extracts the part before the last underscore."""
        return name.rsplit("_", 1)[0] if "_" in name else name

    for group in chains_group:
        for i, chain_id in enumerate(group):
            if i == 0:
                continue
            molecule = [mol for mol in molecule_list if mol.name == chain_id][0]
            for interface in molecule.interface_list:
                interface_id = interface.name
                interface_template_id = interface.my_template.name

                first_appearance = True
                for j in range(i):
                    chain_id_2 = group[j]
                    molecule_2 = [mol for mol in molecule_list if mol.name == chain_id_2][0]
                    for interface_2 in molecule_2.interface_list:
                        if interface_template_id == interface_2.my_template.name:
                            first_appearance = False
                            break
                    if not first_appearance:
                        break

                if first_appearance:
                    my_partner_chain = [c for c in all_chains if c.id == interface_id][0]
                    my_chain = [c for c in all_chains if c.id == chain_id][0]

                    for j in range(i):
                        chain_id_2 = group[j]
                        molecule_2 = [mol for mol in molecule_list if mol.name == chain_id_2][0]
                        for interface_2 in molecule_2.interface_list:
                            interface_template_id_2 = interface_2.my_template.name
                            if interface_template_id != interface_template_id_2:
                                another_partner_chain = [c for c in all_chains if c.id == interface_2.name][0]
                                another_chain = [c for c in all_chains if c.id == chain_id_2][0]

                                R, t = rigid_transform_chains(my_chain, another_chain)

                                my_coords = [res['CA'].coord for res in my_partner_chain if is_aa(res) and 'CA' in res]
                                my_coords_trans = [apply_rigid_transform(R, t, coord) for coord in my_coords]
                                other_coords = [res['CA'].coord for res in another_partner_chain if is_aa(res) and 'CA' in res]

                                if check_clashes_between_two_sets(np.array(my_coords_trans), np.array(other_coords)):
                                    mol_template_id = chains_map[chain_id]
                                    mol_template = [mt for mt in molecule_template_list if mt.name == mol_template_id][0]
                                    
                                    # Match intf1
                                    prefix1 = extract_prefix(interface_template_id)
                                    matches1 = [it for it in mol_template.interface_template_list if extract_prefix(it.name) == prefix1]
                                    if not matches1:
                                        raise ValueError(
                                            f"Interface template '{interface_template_id}' not found in molecule '{mol_template.name}'. "
                                            f"Available: {[it.name for it in mol_template.interface_template_list]}"
                                        )
                                    intf1 = matches1[0]

                                    # Match intf2
                                    prefix2 = extract_prefix(interface_template_id_2)
                                    matches2 = [it for it in mol_template.interface_template_list if extract_prefix(it.name) == prefix2]
                                    if not matches2:
                                        raise ValueError(
                                            f"Interface template '{interface_template_id_2}' not found in molecule '{mol_template.name}'. "
                                            f"Available: {[it.name for it in mol_template.interface_template_list]}"
                                        )
                                    intf2 = matches2[0]

                                    if interface_template_id not in intf2.required_free_list:
                                        intf2.required_free_list.append(interface_template_id)
                                    if interface_template_id_2 not in intf1.required_free_list:
                                        intf1.required_free_list.append(interface_template_id_2)
    return

def _generate_molecule_types(molecule_template_list):
    """Generates molecule types."""

    # Step 1: Generate molecule types
    molecule_types = []
    for mol_template in molecule_template_list:
        mol_name = mol_template.name
        mol_interfaces = []
        for intf_template in mol_template.interface_template_list:
            iface = MoleculeInterface(name=intf_template.name, coord=intf_template.coord)
            mol_interfaces.append(iface)
            print(f"[DEBUG] Template: {intf_template.name}, Coord: {intf_template.coord}")
        molecule = MoleculeType(name=mol_name, interfaces=mol_interfaces, translational_diffusion_constant=mol_template.diffusion_translation, rotational_diffusion_constant=mol_template.diffusion_rotation)
        molecule_types.append(molecule)
    return molecule_types

def _is_existing_mol_temp(molecule_template_name, molecule_template_list):
    """
    Checks if a molecule template with the given name exists in the molecule template list.

    Args:
        molecule_template_name (str): The name of the molecule template to check.
        molecule_template_list (list): The list of molecule template

    Returns:
        tuple: (bool, int or None)
            - True and index if the template exists.
            - False and None otherwise.
    """
    for i, mol_temp in enumerate(molecule_template_list):
        if mol_temp.name == molecule_template_name:
            return True, i
    return False, None

def _is_existing_mol(molecule_name, molecule_list):
    """
    Checks if a molecule with the given name exists in the molecule list.

    Args:
        molecule_name (str): The name of the molecule to check.

    Returns:
        tuple: (bool, int or None)
            - True and index if the molecule exists.
            - False and None otherwise.
    """
    for i, mol in enumerate(molecule_list):
        if mol.name == molecule_name:
            return True, i
    return False, None

def _is_existing_interface(interface_name, interface_list):
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
    for i, interface in enumerate(interface_list):
        if interface.name == interface_name:
            return True, i
    return False, None

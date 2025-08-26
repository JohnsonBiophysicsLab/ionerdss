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
import logging
from copy import deepcopy

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
    find_matching_signature_hash
)
from ionerdss.model.pdb.geometry import rigid_transform_chains
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.utils.coords import Coords
from ionerdss.utils.diffusion_constant import compute_diffusion_constants_nm_us
from ionerdss.utils.rigid_transform import apply_rigid_transform
from ionerdss.utils.steric_clash import check_clashes_between_two_sets

# set up logger
logger = logging.getLogger("ionerdss.model.pdb")       # module-level logger


def process_chain_group(group, chains_map, molecule_template_list,
                        molecule_list, cg_model,
                        params: PDBModelHyperparameters,
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
            params,
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
        molecule.coords = get_chain_data(name, cg_model, "COMs")
        molecule.radius = get_chain_data(name, cg_model, "radii")
        molecule_list.append(molecule)
        template.radius = molecule.radius
    molecule.my_template = template
    molecule.diffusion_translation, molecule.diffusion_rotation = compute_diffusion_constants_nm_us(
        molecule.radius / 10.0)
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
        params: PDBModelHyperparameters,
        interface_template_list, interface_list, binding_chains_pairs,
        signature_to_template_map):
    """
    Processes all interfaces for a given chain within its group.

    Debug prints have been added to trace flow and parameter values.
    
    Args:
        binding_partner_map (map): (chain_idx, iface_idx) -> (partner_chain_idx, partner_iface_idx)
    """

    # debug log the current chain
    logger.debug(
        "Processing chain=%s, group=%s, index=%d", chain_id, group, j)
    logger.debug(
        "Current molecule=%s", molecule.name if hasattr(molecule, 'name') else molecule)

    molecule_name = chain_id

    # get interfaces of this molecule
    # the interfaces are named after partner molecule name
    interfaces = get_chain_data(molecule_name, cg_model, "interfaces")

    # debug log the number of interfaces in the current chain
    logger.debug(
        "Found %d interfaces for chain %s", len(interfaces), molecule_name)

    for i, partner_name in enumerate(interfaces):

        # debug log the current interface
        logger.debug(
            "Start processing new interface...")
        logger.debug(
            "Interface %d: %s -> %s", i, molecule_name, partner_name)

        partner_template_name = chains_map[partner_name]
        logger.debug(
            "Partner=%s, template=%s", partner_name, partner_template_name)

        # create template and mlecule for partner
        partner_template = get_or_create_molecule_template(
            partner_template_name, molecule_template_list)
        partner_molecule = get_or_create_molecule(
            partner_name, partner_template, molecule_list, cg_model)

        # get coordinates of center of mass and interfaces
        com_a = get_chain_data(molecule_name, cg_model, "COMs")
        intf_a = get_chain_data(molecule_name, cg_model, "interface_coords")[i]
        com_b = get_chain_data(partner_name, cg_model, "COMs")

        logger.debug(
            "COM_A=%s, I_A=%s, COM_B=%s", com_a, intf_a, com_b)

        intf_b, res_b, energy_b = None, None, None
        for k, pid in enumerate(get_chain_data(partner_name, cg_model, "interfaces")):
            if pid == molecule_name:
                intf_b = get_chain_data(
                    partner_name, cg_model, "interface_coords")[k]
                res_b = get_chain_data(partner_name, cg_model,
                                       "interface_residues")[k]
                energy_b = get_chain_data(partner_name, cg_model,
                                          "interface_energies")[k]
                logger.debug(
                    "Found reciprocal interface: I_B=%s,\n\tresidues=%s,\n\tenergies=%s", intf_b, res_b, energy_b)
                break

        # build signature hash of the forward and reverse interaction
        # forward interaction = this -> partner
        # reverse interaction = partner -> this
        signature = build_signature(com_a, intf_a, com_b, intf_b)
        inverted_signature = invert_signature(signature)
        sig_key = signature_hash(signature)
        inv_key = signature_hash(inverted_signature)

        logger.debug(
            "Signature hash=%s, inverted=%s", sig_key, inv_key)

        # check if sig_key is in the signature_to_template_map
        # get a subset map of signature_to_template_map that corresponds
        # to the signature hash of the interface on this molecule
        subset_map = {
            k: v
            for k, v in signature_to_template_map.items()
            if hasattr(v, "name") and isinstance(v.name, str)
            and v.name.startswith(chains_map[molecule_name])
        }
        logger.debug(
            "signature_to_template_map subset : \n %s", subset_map)
        matching_key, interface_template = find_matching_signature_hash(
            sig_key, subset_map, params)
        # partner is searched from all interfaces
        matching_inv_key, interface_template = find_matching_signature_hash(
            inv_key, signature_to_template_map, params)
        logger.debug(
            "signature_to_template_map matchkey : \n %s; %s",
            matching_key, matching_inv_key
        )
        logger.debug(
            "signature_to_template_map : \n %s",
            chains_map[molecule_name]
        )
        logger.debug(
            "signature_to_template_map : \n %s",
            signature_to_template_map
        )

        if matching_key is not None and matching_inv_key is not None:
            logger.debug("Reusing existing template")
            interface_template = signature_to_template_map[matching_key]
            partner_template_template = signature_to_template_map[matching_inv_key]
            register_interfaces(
                chain_id, partner_name, molecule, partner_molecule,
                interface_template, partner_template_template,
                cg_model, intf_b, res_b, energy_b, i,
                interface_list, binding_chains_pairs
            )
        else:
            logger.debug("Building new interface templates")
            is_homodimerize = determine_homo_dimerization(
                signature, molecule_name, partner_name, chains_map,
                params
            )
            logger.debug("[Homodimerization? : %s", is_homodimerize)

            interface_templates = build_new_interface_templates(
                is_homodimerize, j, i, group, chain_id, molecule_name,
                partner_name, cg_model, signature, chains_map,
                molecule_template_list, interface_template_list
            )
            logger.debug(
                "Created %d new templates", len(interface_templates))

            # Register templates in the map
            if is_homodimerize:
                signature_to_template_map[sig_key] = interface_templates[0]
                #signature_to_template_map[inv_key] = interface_templates[0]
                logger.debug(
                    "Stored single homodimer template in map")
            else:
                signature_to_template_map[sig_key] = interface_templates[0]
                signature_to_template_map[inv_key] = interface_templates[1]
                logger.debug(
                    "Stored heterodimer templates in map")

            for t in interface_templates:
                logger.debug(
                    "Registering interface template= %s", t)
                register_interfaces(
                    chain_id, partner_name, molecule, partner_molecule,
                    t,
                    partner_template if t != partner_template else t,
                    cg_model, intf_b, res_b, energy_b, i,
                    interface_list, binding_chains_pairs
                )

    # debug log end of processing and separator
    logger.debug(
        "Finished processing chain=%s", chain_id)
    logger.debug(
        "......................................")


def determine_homo_dimerization(signature, mol_a, mol_b,
                                chains_map, params: PDBModelHyperparameters):
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
        abs(signature["dA"] - signature["dB"]) < params.dist_threshold_intra and
        abs(signature["thetaA"] - signature["thetaB"]) < params.angle_threshold
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
        binding_partner_map (map): (chain_idx, iface_idx) -> (partner_chain_idx, partner_iface_idx)

    Returns:
        list of BindingInterfaceTemplate: One or two templates
    """
    def get_interface_template_id(molecule_prefix):
        suffix = sum(1 for mt in mol_template_list for it in mt.interface_template_list if it.name.startswith(
            molecule_prefix)) + 1
        return f"{molecule_prefix}_{suffix}"

    # prepare chain index enumeration
    chain_index = {getattr(c, "id", c): idx for idx, c in enumerate(cg_model["chains"])}


    if is_homo:
        # Case 1 : homodimerize
        # Interface is on mol_A, binding to mol_B (same template)
        molecule_prefix = chains_map[mol_A]
        interface_template_id = get_interface_template_id(molecule_prefix)
        interface_template = BindingInterfaceTemplate(interface_template_id)
        interface_template.signature = signature
        offset = compute_interface_offset(j, i, group, chain_id, cg_model)
        interface_template.coords = offset

        # debug log the offset of the current interface
        logger.debug("interface %s offset = %s", interface_template_id, offset)

        # Assign to correct molecule template (chain_id is mol_A)
        my_template = get_or_create_molecule_template(
            chains_map[chain_id], mol_template_list)
        my_template.interface_template_list.append(interface_template)
        interface_template_list.append(interface_template)
        return [interface_template]

    else:
        # Case 2 : two species binding
        templates = []
        for mol, sig, side_chain_id in [
            (mol_A, signature, chain_id),
            (mol_B, invert_signature(signature),
             mol_B if chain_id == mol_A else mol_A)
        ]:
            gidx = chain_index[chain_id]     # chain_id like "L"
                                             # note gidx is global index
                                             # whereas j is local index
                                             # within chain group
            # determine current group
            # determine current interface in terms of (chain id, interface id)
            if mol == mol_A:
                molecule_prefix = chains_map[mol_A]
                current_chain_id = gidx
                current_interface_id = i
            else:
                molecule_prefix = chains_map[mol_B]
                current_chain_id, current_interface_id = cg_model["binding_partner_map"][(gidx,i)]

            # get all keys with value == current group prefix
            current_group = [k for k, v in chains_map.items() if v == molecule_prefix]
            interface_template_id = get_interface_template_id(molecule_prefix)
            interface_template = BindingInterfaceTemplate(
                interface_template_id)
            interface_template.signature = sig
            offset = compute_interface_offset(
                current_chain_id, current_interface_id,
                current_group, side_chain_id, cg_model)
            interface_template.coords = offset

            # debug log the offset of the current interface
            logger.debug("interface %s offset = %s", interface_template_id, offset)

            # Attach to molecule template that interface lives on
            my_template = get_or_create_molecule_template(
                chains_map[side_chain_id], mol_template_list)
            my_template.interface_template_list.append(interface_template)
            interface_template_list.append(interface_template)
            templates.append(interface_template)

        return templates

def compute_interface_offset(j, i, group, chain_id, cg_model):
    """
    Computes the relative coordinate of the interface site with respect to COM.
    Applies rigid transform if not the first chain in the group.

    Args:
        j (int): Chain index within the group.
        i (int): Interface index within the chain.
        group (list): Group of chains (list of chain IDs).
        chain_id (str): Chain ID of current molecule.
        cg_model (dict): CG data containing coords.

    Returns:
        Coords: Relative offset from COM to interface site.
    """

    logger.debug(
        "[compute_interface_offset] START chain_id=%s, j=%d, i=%d, group=%s",
        chain_id, j, i, group
    )

    if j == 0:
        pt = get_chain_data(chain_id, cg_model, "interface_coords")[i]
        com = get_chain_data(chain_id, cg_model, "COMs")
        offset = np.array([pt.x, pt.y, pt.z]) - np.array([com.x, com.y, com.z])

        logger.debug(
            "[compute_interface_offset] j==0 (reference chain)\n"
            "  COM = (%f, %f, %f)\n"
            "  Pt  = (%f, %f, %f)\n"
            "  Offset = (%f, %f, %f)",
            com.x, com.y, com.z,
            pt.x, pt.y, pt.z,
            offset[0], offset[1], offset[2]
        )
        return Coords(*offset)

    else:
        chain1 = [c for c in cg_model["chains"] if c.id == group[0]][0]
        chain2 = [c for c in cg_model["chains"] if c.id == chain_id][0]

        logger.debug(
            "[compute_interface_offset] j>0 (aligned chain)\n"
            "  Reference chain = %s\n"
            "  Current chain   = %s",
            chain1.id, chain2.id
        )

        # Compute rigid transform from current chain to reference chain
        R, t = rigid_transform_chains(chain2, chain1)
        logger.debug(
            "[compute_interface_offset] Rigid transform (R,t):\n"
            "  R = %s\n"
            "  t = %s",
            np.array2string(R, precision=3), np.array2string(t, precision=3)
        )

        Q_COM = get_chain_data(chain_id, cg_model, "COMs")
        Q_pt = get_chain_data(chain_id, cg_model, "interface_coords")[i]
        logger.debug(
            "[compute_interface_offset] Raw COM = %s, Raw Pt = %s",
            Q_COM, Q_pt
        )

        # Apply rigid transform to both COM and interface point
        Q2 = [apply_rigid_transform(R, t, pt.to_numpy())
              for pt in [Q_COM, Q_pt]]

        logger.debug(
            "[compute_interface_offset] Transformed COM = %s\n"
            "  Transformed Pt  = %s",
            Q2[0], Q2[1]
        )

        assert isinstance(Q_COM, Coords), f"Expected Coords, got {type(Q_COM)}"

        offset = Q2[1] - Q2[0]
        logger.debug(
            "[compute_interface_offset] Final Offset = (%f, %f, %f)",
            offset[0], offset[1], offset[2]
        )
        return Coords(*offset)



def match_existing_templates(signature, signature_conj, mol_template_list, params):
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
            if signature_are_similar(sig, it.signature, params)
        ]
        if not matches:
            raise ValueError(
                f"No matching interface template for signature: {sig}")
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
        interface.coords = get_chain_data(
            chain_id, cg_model, "interface_coords")[i]
        interface.my_residues = get_chain_data(
            chain_id, cg_model, "interface_residues")[i]
        interface.energy = get_chain_data(
            chain_id, cg_model, "interface_energies")[i]
        interface.my_template.energy = interface.energy
        interface_list.append(interface)
        molecule.interface_list.append(interface)

        partner_interface = BindingInterface(chain_id)
        partner_interface.my_template = partner_template
        partner_interface.coords = I_B
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
            mol_template = next(
                mt for mt in mol_template_list if mt.name == chains_map[chain_id])
            mol0 = next(m for m in molecule_list if m.name == group[0])
            com_coords = mol0.coords
            interface_coords = [
                it.coords + com_coords for it in mol_template.interface_template_list]
            interface_ids = [
                it.name for it in mol_template.interface_template_list]

            if i == 0:
                mol = next(m for m in molecule_list if m.name == chain_id)
                mol.normal_point = [com_coords.x, com_coords.y, com_coords.z + 1]
                continue

            chain1 = [c for c in cg_model["chains"] if c.id == group[0]][0]
            chain2 = [c for c in cg_model["chains"] if c.id == chain_id][0]
            R, t = rigid_transform_chains(chain1, chain2)
            com_trans = apply_rigid_transform(R, t, np.array(
                [com_coords.x, com_coords.y, com_coords.z]))
            intf_coords_trans = [apply_rigid_transform(
                R, t, np.array([pt.x, pt.y, pt.z])) for pt in interface_coords]
            norm_pt = apply_rigid_transform(R, t, np.array(
                [com_coords.x, com_coords.y, com_coords.z + 1]))

            mol = next(m for m in molecule_list if m.name == chain_id)
            mol.coords = Coords(*com_trans)
            mol.normal_point = list(norm_pt)
            for j, interface in enumerate(mol.interface_list):
                for k, tid in enumerate(interface_ids):
                    if interface.my_template.name == tid:
                        interface.coords = Coords(*intf_coords_trans[k])
                        break


def regularize_repeated_chains(cg_model, chains_map, chains_groups,
                               params):
    """
    Aligns and regularizes all molecular chains so that homologous chains share 
    the same relative geometry. This method organizes molecule and interface objects 
    accordingly and sets up reaction objects.
    """
    if not chains_groups:
        raise ValueError("Invalid chains group!")
    for group in chains_groups:
        group.sort()
    chains_groups.sort()

    # ... [rest of the regularization logic remains mostly unchanged] ...

    # check if the structure has homologous chains
    # if any element in chains_group has more than one chain, then it has homologous chains
    has_repeated_chains = any(len(group) > 1 for group in chains_groups)
    if not has_repeated_chains:
        logger.debug("NO REPEATED CHAINS! OVERRIDE THRESHOLDS TO ZERO!")

        # try not to change the passed-in values via initializing a new hyperparameter class
        options = {"dist_threshold_intra": 0.0,
                   "dist_threshold_inter": 0.0,
                   "angle_threshold": 0.0}
        params = PDBModelHyperparameters(options=options)

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
            params, signature_to_template_map
        )
        all_interface_signatures.extend(sigs)
        interface_template_list.extend(templates)
        # debug log
        logger.debug("interface templates : %s \n", templates)
        logger.debug("All interface templates : %s \n",
                     interface_template_list)
        logger.debug("All interface coords : %s \n",
                     [it.coords for it in interface_template_list])

        interface_list.extend(interfaces)
        all_binding_pairs.extend(pairs)

    # Post-processing:
    # this function updates the regularized positions in-place
    # on `molecule_list` (and each molecule’s interface_list).
    # The original, raw coarse-grained data in cg_model remains
    # untouched—so the “original copy” is intact.
    update_molecule_interfaces(
        chains_groups, chains_map,
        molecule_template_list, molecule_list,
        cg_model
    )

    # update requirement list
    _update_interface_templates_free_required_list(chains_groups,
                                                   molecule_list,
                                                   cg_model["chains"],
                                                   chains_map, molecule_template_list)

    all_binding_pairs.sort()
    molecule_list.sort(key=lambda m: m.name)
    molecule_template_list.sort(key=lambda mt: mt.name)
    interface_list.sort(key=lambda i: i.name)
    interface_template_list.sort(key=lambda it: it.name)


    logger.debug(
        "Molecules Template and Reactions Template After Regularization:")
    for molecule_template in molecule_template_list:
        logger.debug(molecule_template)

    logger.debug("Molecules:")
    for molecule in molecule_list:
        logger.debug(molecule)

    molecule_types = _generate_molecule_types(molecule_template_list)

    # get updated cg_model
    updated_cg_model = build_updated_cg_model(cg_model, molecule_list)
    logger.debug("updated cg_model : \n%s", updated_cg_model)

   # Return structured model data
    return {
        "updated_cg_model": updated_cg_model,
        "molecule_templates": molecule_template_list,
        "molecules": molecule_list,
        "interface_templates": interface_template_list,
        "interfaces": interface_list,
        "binding_pairs": all_binding_pairs,
        "molecule_types": molecule_types,
    }

# This is a modularized scaffold for the original parsing loop
# The functions are broken down by responsibility

def build_updated_cg_model(cg_model, molecule_list):
    """
    Create a new cg_model dict whose COMs and interface_coords are replaced
    with the regularized positions from `molecule_list`. All other fields are
    copied from the original cg_model. The original cg_model is not modified.

    Returns
    -------
    dict
        A deepcopy of cg_model with updated 'COMs' and 'interface_coords'.
    """
    updated = deepcopy(cg_model)

    # Fast lookups
    mol_by_name = {m.name: m for m in molecule_list}
    chains = cg_model["chains"]                   # list of Bio.PDB Chain objects
    interfaces = cg_model["interfaces"]           # list[list[str]] partner IDs

    new_COMs = []
    new_interface_coords = []

    for idx, chain in enumerate(chains):
        chain_id = chain.id
        mol = mol_by_name.get(chain_id)

        # Fallback to original data if we don't have a regularized molecule
        if mol is None or mol.coords is None:
            new_COMs.append(cg_model["COMs"][idx])
            new_interface_coords.append(cg_model["interface_coords"][idx])
            continue

        # 1) COM from regularized molecule
        new_COMs.append(mol.coords)  # Coords

        # 2) Interface coords in the same partner order as cg_model['interfaces'][idx]
        partner_order = interfaces[idx]
        # Build dict partner_name -> Coords from the molecule's interfaces
        intf_coord_by_partner = {
            intf.name: intf.coords for intf in getattr(mol, "interface_list", [])}

        # Keep array length identical to original; fallback to old if missing
        old_chain_iface_coords = cg_model["interface_coords"][idx]
        coords_for_chain = []
        for j, partner in enumerate(partner_order):
            coords_for_chain.append(
                intf_coord_by_partner.get(partner, old_chain_iface_coords[j])
            )
        new_interface_coords.append(coords_for_chain)

    updated["COMs"] = new_COMs
    updated["interface_coords"] = new_interface_coords
    return updated


def _update_interface_templates_free_required_list(
    chains_group, molecule_list, all_chains, chains_map, molecule_template_list
):
    """
    Updates the `required_free_list` attribute for each interface template by checking 
    potential steric clashes among binding partners within the same molecule template.
    """
    logger.debug("\nStarting _update_interface_templates_free_required_list")
    logger.debug("---------------------------------------------------------------")
    logger.debug("Number of chain groups: %d", len(chains_group))
    for idx, group in enumerate(chains_group):
        logger.debug("  Group %d: %s", idx, group)
    logger.debug("Total molecules: %d", len(molecule_list))
    logger.debug("  Molecule names: %s", [m.name for m in molecule_list])
    logger.debug("Total chains in cg_model: %d", len(all_chains))
    logger.debug("  All chain IDs: %s", [chain.id for chain in all_chains])
    logger.debug("Total molecule templates: %d", len(molecule_template_list))
    for mt in molecule_template_list:
        logger.debug(
            "  MoleculeTemplate '%s' has interface templates: %s",
            mt.name, [it.name for it in mt.interface_template_list]
        )
    logger.debug("Chains map (chain_id -> molecule_template):")
    for k, v in chains_map.items():
        logger.debug("  %s -> %s", k, v)
    logger.debug("---------------------------------------------------------------\n")


    def extract_prefix(name):
        """Extracts the part before the last underscore."""
        return name.rsplit("_", 1)[0] if "_" in name else name

    def get_mol_template(name: str):
        return next(mt for mt in molecule_template_list if mt.name == name)

    def get_interface_template_by_partner_prefix(mol_template, partner_prefix: str):
        """
        Find an interface template on `mol_template` whose name prefix (before '_')
        matches `partner_prefix` (which is a molecule template name like 'A' or 'H').
        """
        matches = [it for it in mol_template.interface_template_list
                   if extract_prefix(it.name) == partner_prefix]
        if not matches:
            raise ValueError(
                f"Interface template with partner '{partner_prefix}' not found in molecule '{mol_template.name}'. "
                f"Available: {[it.name for it in mol_template.interface_template_list]}"
            )
        return matches[0]

    for group in chains_group:
        for i, chain_id in enumerate(group):
            if i == 0:
                continue

            molecule = next(
                mol for mol in molecule_list if mol.name == chain_id)

            for interface in molecule.interface_list:
                # Partner chain for THIS interface
                # partner chain id (e.g., 'A' or 'H')
                interface_id = interface.name

                # --- FIRST APPEARANCE CHECK (by partner prefix) ---
                # Determine the partner template prefix of THIS interface
                # e.g., 'A' or 'H'
                this_partner_prefix = chains_map[interface_id]

                first_appearance = True
                for j in range(i):
                    chain_id_2 = group[j]
                    molecule_2 = next(
                        m for m in molecule_list if m.name == chain_id_2)
                    for interface_2 in molecule_2.interface_list:
                        # Compare by partner template prefix (not by .my_template.name)
                        other_partner_prefix = chains_map[interface_2.name]
                        if other_partner_prefix == this_partner_prefix:
                            first_appearance = False
                            break
                    if not first_appearance:
                        break

                if not first_appearance:
                    continue

                # Get chains needed for clash check
                my_partner_chain = next(
                    c for c in all_chains if c.id == interface_id)
                my_chain = next(c for c in all_chains if c.id == chain_id)

                for j in range(i):
                    chain_id_2 = group[j]
                    molecule_2 = next(
                        m for m in molecule_list if m.name == chain_id_2)

                    for interface_2 in molecule_2.interface_list:
                        # Skip same-partner-prefix interfaces? (optional)
                        # if chains_map[interface_2.name] == this_partner_prefix:
                        #     continue

                        another_partner_chain = next(
                            c for c in all_chains if c.id == interface_2.name)
                        another_chain = next(
                            c for c in all_chains if c.id == chain_id_2)

                        R, t = rigid_transform_chains(my_chain, another_chain)

                        my_coords = [res['CA'].coords for res in my_partner_chain if is_aa(
                            res) and 'CA' in res]
                        my_coords_trans = [apply_rigid_transform(
                            R, t, coords) for coords in my_coords]
                        other_coords = [
                            res['CA'].coords for res in another_partner_chain if is_aa(res) and 'CA' in res]

                        if check_clashes_between_two_sets(np.array(my_coords_trans), np.array(other_coords)):
                            # Resolve the molecule template for THIS chain
                            # e.g., 'H'
                            mol_template_id = chains_map[chain_id]
                            mol_template = get_mol_template(mol_template_id)

                            # Resolve interface templates on this molecule template
                            # by their PARTNER prefixes.
                            # intf1: the template on `mol_template` that binds THIS interface’s partner
                            intf1 = get_interface_template_by_partner_prefix(
                                mol_template, this_partner_prefix
                            )
                            # intf2: the template on `mol_template` that binds interface_2’s partner
                            partner_prefix_2 = chains_map[interface_2.name]
                            intf2 = get_interface_template_by_partner_prefix(
                                mol_template, partner_prefix_2
                            )

                            # Update required_free_list symmetrically by full names
                            if intf1.name not in intf2.required_free_list:
                                intf2.required_free_list.append(intf1.name)
                            if intf2.name not in intf1.required_free_list:
                                intf1.required_free_list.append(intf2.name)
    return


def _generate_molecule_types(molecule_template_list):
    """Generates molecule types."""

    # Step 1: Generate molecule types
    molecule_types = []
    for mol_template in molecule_template_list:
        mol_name = mol_template.name
        mol_interfaces = []
        for intf_template in mol_template.interface_template_list:
            iface = MoleculeInterface(
                name=intf_template.name, coords=intf_template.coords)
            mol_interfaces.append(iface)
            logger.debug(
                "Template: %s", intf_template.name)
            logger.debug(
                "Coords: %s", intf_template.coords)
        molecule = MoleculeType(name=mol_name, interfaces=mol_interfaces, translational_diffusion_constant=mol_template.diffusion_translation,
                                rotational_diffusion_constant=mol_template.diffusion_rotation)
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

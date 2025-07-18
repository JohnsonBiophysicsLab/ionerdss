"""
repeated_chain_alignment.py

TODO: > change to regularize_repeats

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

from ionerdss.model.components import (
    MoleculeTemplate,
    CoarseGrainedMolecule,
    BindingInterface,
)
from ionerdss.math.diffusion_constant import compute_diffusion_constants_nm_us

def regularize_homologous_chains(chains_group,
                                 chains_map,
                                 all_chains,
                                 all_COM_chains_coords,
                                 all_chains_radius,
                                 dist_threshold_intra=3.5,
                                 dist_threshold_inter=3.5,
                                 angle_threshold=25.0,
                                 show_coarse_grained_structure=False,
                                 save_pymol_script=False,
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
    if not chains_group:
        raise ValueError("Invalid chains group!")
    for group in chains_group:
        group.sort()
    chains_group.sort()

    # check if the structure has homologous chains
    # if any element in chains_group has more than one chain, then it has homologous chains
    has_homologous_chains = any(len(group) > 1 for group in chains_group)
    if not has_homologous_chains:
        dist_threshold_intra = 0.0
        dist_threshold_inter = 0.0
        angle_threshold = 0.0

    molecule_list = []
    molecule_template_list = []
    interface_list = []
    interface_template_list = []
    interface_signatures = []

    for group in chains_group:
        # print(f"Start parsing chain group / molecule template {group}")
        molecule_template_name = chains_map[group[0]]
        is_existing_mol_temp, idx = _is_existing_mol_temp(molecule_template_name)
        if is_existing_mol_temp:
            # print(f"This is an existed mol template {molecule_template_name}")
            molecule_template = molecule_template_list[idx]
        else:
            molecule_template = MoleculeTemplate(molecule_template_name)
            # print(f"New mol template {molecule_template_name} is created.")
            molecule_template_list.append(molecule_template)

        for j, chain_id in enumerate(group):
            # print(f"Start parsing chain / molecule {chain_id}")
            molecule_name = chain_id
            is_existing_mol, mol_index = _is_existing_mol(molecule_name)
            if is_existing_mol:
                # print(f"This is an existing molecule {molecule_name}")
                molecule = molecule_list[mol_index]
                molecule.radius = all_chains_radius[all_chains.index([chain for chain in all_chains if chain.id == molecule_name][0])]
                molecule.diffusion_translation, molecule.diffusion_rotation = compute_diffusion_constants_nm_us(molecule.radius / 10.0)
                molecule.my_template.diffusion_translation, molecule.my_template.diffusion_rotation = molecule.diffusion_translation, molecule.diffusion_rotation
            else:
                molecule = CoarseGrainedMolecule(molecule_name)
                # print(f"New molecule {molecule_name} is created.")
                molecule.my_template = molecule_template
                molecule.coord = all_COM_chains_coords[all_chains.index([chain for chain in all_chains if chain.id == molecule_name][0])]
                molecule.radius = all_chains_radius[all_chains.index([chain for chain in all_chains if chain.id == molecule_name][0])]
                molecule.diffusion_translation, molecule.diffusion_rotation = compute_diffusion_constants_nm_us(molecule.radius / 10.0)
                molecule_list.append(molecule)
                molecule_template.radius = molecule.radius
                molecule_template.diffusion_translation, molecule_template.diffusion_rotation = molecule.diffusion_translation, molecule.diffusion_rotation
            
            # loop the interface of this chain (molecule)
            for i, interface_id in enumerate(all_interfaces[all_chains.index([chain for chain in all_chains if chain.id == molecule_name][0])]):
                A = molecule_name
                B = interface_id # this is the chain name of the partner
                partner_mol_template_name = chains_map[B]
                # print(f"Parsing the interface {interface_id} for molecule {molecule_name}; its binding partner is molecule {B} via its interface {A}")
                is_existing_mol_temp, idx = _is_existing_mol_temp(partner_mol_template_name)
                if is_existing_mol_temp:
                    # print(f"molecule {B} already has its template created.")
                    partner_molecule_template = molecule_template_list[idx]
                else:
                    partner_molecule_template = MoleculeTemplate(partner_mol_template_name)
                    # print(f"new mol template {partner_mol_template_name} created for molecule {B}.")
                    molecule_template_list.append(partner_molecule_template)

                is_existing_mol, partner_mol_index = _is_existing_mol(B)
                if is_existing_mol:
                    # print(f"molecule {B} is already created.")
                    partner_molecule = molecule_list[partner_mol_index]
                else:
                    partner_molecule = CoarseGrainedMolecule(B)
                    # print(f"New molecule {B} is created.")
                    partner_molecule.my_template = partner_molecule_template
                    partner_molecule.coord = all_COM_chains_coords[all_chains.index([chain for chain in all_chains if chain.id == B][0])]
                    molecule_list.append(partner_molecule)

                COM_A = all_COM_chains_coords[all_chains.index([chain for chain in all_chains if chain.id == A][0])]
                I_A = all_interfaces_coords[all_chains.index([chain for chain in all_chains if chain.id == A][0])][i]
                COM_B = all_COM_chains_coords[all_chains.index([chain for chain in all_chains if chain.id == B][0])]
                for k, partner_interface_id in enumerate(all_interfaces[all_chains.index([chain for chain in all_chains if chain.id == B][0])]):
                    if partner_interface_id == A:
                        I_B = all_interfaces_coords[all_chains.index([chain for chain in all_chains if chain.id == B][0])][k]
                        R_B = all_interfaces_residues[all_chains.index([chain for chain in all_chains if chain.id == B][0])][k]
                        E_B = all_interface_energies[all_chains.index([chain for chain in all_chains if chain.id == B][0])][k]
                        break

                signature = {
                    "dA": np.linalg.norm([(COM_A - I_A).x, (COM_A - I_A).y, (COM_A - I_A).z]),
                    "dB": np.linalg.norm([(COM_B - I_B).x, (COM_B - I_B).y, (COM_B - I_B).z]),
                    "dAB": np.linalg.norm([(I_A - I_B).x, (I_A - I_B).y, (I_A - I_B).z]),
                    "thetaA": _calc_angle(COM_A, I_A, I_B),
                    "thetaB": _calc_angle(COM_B, I_B, I_A)
                }

                # print the signature
                # print(f"Parsing signature: {signature}")

                is_existing_sig = False

                for existing_sig in interface_signatures:
                    if _sig_are_similar(signature, existing_sig, dist_threshold_intra, dist_threshold_inter, angle_threshold):
                        is_existing_sig = True
                        break

                if not is_existing_sig:
                    # print("this is a new signature. added to list.")
                    interface_signatures.append(signature)
                    signature_conjugated = {
                        "dA": signature["dB"],
                        "dB": signature["dA"],
                        "dAB": signature["dAB"],
                        "thetaA": signature["thetaB"],
                        "thetaB": signature["thetaA"]
                    }
                    interface_signatures.append(signature_conjugated)
                    # print(f"the conjugated signature: {signature_conjugated} is also added to the list.")

                    # build the interface template pairs for both molecule templates, need to check if this is homo dimerization or hetero
                    is_homo = False
                    if chains_map[A] != chains_map[B]:
                        pass
                    else:
                        if abs(signature["dA"] - signature["dB"]) > dist_threshold_intra or abs(signature["thetaA"] - signature["thetaB"]) > angle_threshold:
                            pass
                        else:
                            is_homo = True

                    if is_homo:
                        # only need to build the interface template once
                        interface_template_id_prefix = chains_map[A]

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
                            interface_template.coord = all_interfaces_coords[all_chains.index([chain for chain in all_chains if chain.id == chain_id][0])][i] - molecule.coord
                        else:
                            # align the current chain to the first chain in the group, then get the relative position of interface to COM
                            chain1 = all_chains[all_chains.index([chain for chain in all_chains if chain.id == group[0]][0])]
                            chain2 = all_chains[all_chains.index([chain for chain in all_chains if chain.id == chain_id][0])]
                            R, t = rigid_transform_chains(chain2, chain1)
                            Q = []
                            Q_COM_coord = all_COM_chains_coords[all_chains.index([chain for chain in all_chains if chain.id == chain_id][0])]
                            Q.append([Q_COM_coord.x, Q_COM_coord.y, Q_COM_coord.z])
                            temp_coord = all_interfaces_coords[all_chains.index([chain for chain in all_chains if chain.id == chain_id][0])][i]
                            Q.append([temp_coord.x, temp_coord.y, temp_coord.z])
                            Q2 = []
                            for point in Q:
                                transformed_point = apply_rigid_transform(R, t, np.array(point))
                                Q2.append(transformed_point)
                            interface_template.coord = Coords(Q2[1][0] - Q2[0][0], Q2[1][1] - Q2[0][1], Q2[1][2] - Q2[0][2])
                        molecule_template.interface_template_list.append(interface_template)
                        interface_template_list.append(interface_template)
                        partner_interface_template = interface_template
                        partner_molecule_template = molecule_template
                    else:
                        # add interface template 1
                        interface_template_id_prefix = chains_map[B]

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
                            interface_template.coord = all_interfaces_coords[all_chains.index([chain for chain in all_chains if chain.id == chain_id][0])][i] - molecule.coord
                        else:
                            # align the current chain to the first chain in the group, then get the relative position of interface to COM
                            chain1 = all_chains[all_chains.index([chain for chain in all_chains if chain.id == group[0]][0])]
                            chain2 = all_chains[all_chains.index([chain for chain in all_chains if chain.id == chain_id][0])]
                            R, t = rigid_transform_chains(chain2, chain1)
                            Q = []
                            Q_COM_coord = all_COM_chains_coords[all_chains.index([chain for chain in all_chains if chain.id == chain_id][0])]
                            Q.append([Q_COM_coord.x, Q_COM_coord.y, Q_COM_coord.z])
                            temp_coord = all_interfaces_coords[all_chains.index([chain for chain in all_chains if chain.id == chain_id][0])][i]
                            Q.append([temp_coord.x, temp_coord.y, temp_coord.z])
                            Q2 = []
                            for point in Q:
                                transformed_point = apply_rigid_transform(R, t, np.array(point))
                                Q2.append(transformed_point)
                            interface_template.coord = Coords(Q2[1][0] - Q2[0][0], Q2[1][1] - Q2[0][1], Q2[1][2] - Q2[0][2])
                        molecule_template.interface_template_list.append(interface_template)
                        interface_template_list.append(interface_template)

                        # add interface template 2
                        interface_template_id_prefix = chains_map[A]

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
                        for g in chains_group:
                            if B in g:
                                B_group = g

                        if B == B_group[0]:
                            partner_interface_template.coord = I_B - partner_molecule.coord
                        else:
                            # align the current chain to the first chain in the group, then get the relative position of interface to COM
                            chain1 = all_chains[all_chains.index([chain for chain in all_chains if chain.id == B_group[0]][0])]
                            chain2 = all_chains[all_chains.index([chain for chain in all_chains if chain.id == B][0])]
                            R, t = rigid_transform_chains(chain2, chain1)
                            Q = []
                            Q_COM_coord = all_COM_chains_coords[all_chains.index([chain for chain in all_chains if chain.id == B][0])]
                            Q.append([Q_COM_coord.x, Q_COM_coord.y, Q_COM_coord.z])
                            temp_coord = I_B
                            Q.append([temp_coord.x, temp_coord.y, temp_coord.z])
                            Q2 = []
                            for point in Q:
                                transformed_point = apply_rigid_transform(R, t, np.array(point))
                                Q2.append(transformed_point)
                            partner_interface_template.coord = Coords(Q2[1][0] - Q2[0][0], Q2[1][1] - Q2[0][1], Q2[1][2] - Q2[0][2])
                        partner_molecule_template.interface_template_list.append(partner_interface_template)
                        interface_template_list.append(partner_interface_template)

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
                    
                    # Find all matching templates for signature
                    matching_templates = []
                    for mol_temp in molecule_template_list:
                        for interface_temp in mol_temp.interface_template_list:
                            if _sig_are_similar(signature, interface_temp.signature, dist_threshold_intra, dist_threshold_inter, angle_threshold):
                                matching_templates.append((interface_temp, mol_temp))

                    # Check for errors in template matching
                    if len(matching_templates) == 0:
                        # Print all available signatures for debugging
                        available_sigs = []
                        for mol_temp in molecule_template_list:
                            for interface_temp in mol_temp.interface_template_list:
                                available_sigs.append({
                                    'template': f"{interface_temp.name} ({mol_temp.name})",
                                    'signature': interface_temp.signature
                                })
                        
                        print(f"Target signature: {signature}")
                        print("Available signatures:")
                        for sig_info in available_sigs:
                            print(f"  {sig_info['template']}: {sig_info['signature']}")
                        
                        raise ValueError(
                            f"No matching interface template found for signature: {signature}\n"
                            f"Available templates: {[(mt.name, len(mt.interface_template_list)) for mt in molecule_template_list]}\n"
                            f"Current thresholds: dist_intra={dist_threshold_intra}, dist_inter={dist_threshold_inter}, angle={angle_threshold}\n"
                            f"Please increase the thresholds to find a match."
                        )
                    elif len(matching_templates) > 1:
                        matching_templates.sort(key=lambda pair: _sig_difference(signature, pair[0].signature))
                        interface_template, molecule_template = matching_templates[0]
                        print(f"Multiple matches found. Using closest match: {interface_template.name} ({molecule_template.name})")
                    else:
                        interface_template, molecule_template = matching_templates[0]

                    # Find all matching templates for conjugated signature
                    matching_partner_templates = []
                    for mol_temp in molecule_template_list:
                        for interface_temp in mol_temp.interface_template_list:
                            if _sig_are_similar(signature_conjugated, interface_temp.signature, dist_threshold_intra, dist_threshold_inter, angle_threshold):
                                matching_partner_templates.append((interface_temp, mol_temp))

                    # Check for errors in partner template matching
                    if len(matching_partner_templates) == 0:
                        # Print all available signatures for debugging
                        available_sigs = []
                        for mol_temp in molecule_template_list:
                            for interface_temp in mol_temp.interface_template_list:
                                available_sigs.append({
                                    'template': f"{interface_temp.name} ({mol_temp.name})",
                                    'signature': interface_temp.signature
                                })
                        
                        print(f"Target conjugated signature: {signature_conjugated}")
                        print("Available signatures:")
                        for sig_info in available_sigs:
                            print(f"  {sig_info['template']}: {sig_info['signature']}")
                        
                        raise ValueError(
                            f"No matching partner interface template found for conjugated signature: {signature_conjugated}\n"
                            f"Available templates: {[(mt.name, len(mt.interface_template_list)) for mt in molecule_template_list]}\n"
                            f"Current thresholds: dist_intra={dist_threshold_intra}, dist_inter={dist_threshold_inter}, angle={angle_threshold}\n"
                            f"Please adjust the thresholds to find a match."
                        )
                    elif len(matching_partner_templates) > 1:
                        matching_partner_templates.sort(key=lambda pair: _sig_difference(signature_conjugated, pair[0].signature))
                        partner_interface_template, partner_molecule_template = matching_partner_templates[0]
                        print(f"Multiple conjugated matches found. Using closest match: {partner_interface_template.name} ({partner_molecule_template.name})")
                    else:
                        partner_interface_template, partner_molecule_template = matching_partner_templates[0]

                # build the interfaces for molecules, link the interface template to interface

                is_existing_interface, _ = _is_existing_interface(interface_id, molecule)

                if not is_existing_interface:
                    # print(f"Creating new interface {interface_id} for molecule {molecule_name}")
                    # create the interface
                    interface = BindingInterface(B)
                    interface.my_template = interface_template
                    interface.coord = all_interfaces_coords[all_chains.index([chain for chain in all_chains if chain.id == A][0])][i]
                    interface.my_residues = all_interfaces_residues[all_chains.index([chain for chain in all_chains if chain.id == A][0])][i]
                    interface.energy = all_interface_energies[all_chains.index([chain for chain in all_chains if chain.id == A][0])][i]
                    interface.my_template.energy = interface.energy
                    interface_list.append(interface)
                    molecule.interface_list.append(interface)

                    # print(f"Creating new interface {A} for partner molecule {B}")
                    # create the interface for the partner molecule
                    partner_interface = BindingInterface(A)
                    partner_interface.my_template = partner_interface_template
                    partner_interface.coord = I_B
                    partner_interface.my_residues = R_B
                    partner_interface.energy = E_B
                    partner_interface.my_template.energy = E_B
                    interface_list.append(partner_interface)
                    partner_molecule.interface_list.append(partner_interface)

                    # add the chains pair to binding_chains_pairs
                    if chain_id < interface_id:
                        binding_chains_pair = (chain_id, interface_id)
                    else:
                        binding_chains_pair = (interface_id, chain_id)
                    if binding_chains_pair not in binding_chains_pairs:
                        binding_chains_pairs.append(binding_chains_pair)
                else:
                    # print(f"Interface {interface_id} already exists for molecule {molecule_name}")
                    # print(f"Interface {A} already exists for molecule {B}")
                    pass

    # update the interfaces list of each molecule based on the molecule template
    for group in chains_group:
        for i, chain_id in enumerate(group):
            # determin the COM and interfaces of the corresponding molecule template
            molecule_template = [mol_template for mol_template in molecule_template_list if mol_template.name == chains_map[chain_id]][0]
            molecule_0 = [mol for mol in molecule_list if mol.name == group[0]][0]
            com_coord = molecule_0.coord
            interface_coords = [interface_template.coord + com_coord for interface_template in molecule_template.interface_template_list]
            interface_template_ids = [interface_template.name for interface_template in molecule_template.interface_template_list]

            # calculate the R and t for the rigid transformation
            if i == 0:
                # calculate the normal_point for this molecule
                molecule = [mol for mol in molecule_list if mol.name == chain_id][0]
                molecule.normal_point = [com_coord.x, com_coord.y, com_coord.z + 1] # normal_point - COM is [0,0,1]
                # no need to transform the first chain
                continue
            else:
                chain1 = all_chains[all_chains.index([chain for chain in all_chains if chain.id == group[0]][0])]
                chain2 = all_chains[all_chains.index([chain for chain in all_chains if chain.id == chain_id][0])]
                R, t = rigid_transform_chains(chain1, chain2)
                com_coord_transformed = apply_rigid_transform(R, t, np.array([com_coord.x, com_coord.y, com_coord.z]))
                interface_coords_transformed = []
                for interface_coord in interface_coords:
                    interface_coord_transformed = apply_rigid_transform(R, t, np.array([interface_coord.x, interface_coord.y, interface_coord.z]))
                    interface_coords_transformed.append(interface_coord_transformed)
                normal_point_transformed = apply_rigid_transform(R, t, np.array([com_coord.x, com_coord.y, com_coord.z + 1]))
                # update the COM and interfaces of the molecule
                molecule = [mol for mol in molecule_list if mol.name == chain_id][0]
                molecule.coord = Coords(com_coord_transformed[0], com_coord_transformed[1], com_coord_transformed[2])
                for j, interface in enumerate(molecule.interface_list):
                    # find the corresponding interface template
                    interface_template_id = interface.my_template.name
                    for k, intf_template in enumerate(interface_template_ids):
                        if interface_template_id == intf_template:
                            interface.coord = Coords(interface_coords_transformed[k][0], interface_coords_transformed[k][1], interface_coords_transformed[k][2])
                            break
                molecule.normal_point = [normal_point_transformed[0], normal_point_transformed[1], normal_point_transformed[2]]

    _update_interface_templates_free_required_list()

    binding_chains_pairs.sort()
    molecule_list.sort(key=lambda m: m.name)
    molecule_template_list.sort(key=lambda mt: mt.name)
    interface_list.sort(key=lambda i: i.name)
    interface_template_list.sort(key=lambda it: it.name)

    # print("binding chains pairs:")
    # for pair in binding_chains_pairs:
    #     print(pair)
    # print("molecule list:")
    # for molecule in molecule_list:
    #     print(molecule)
    # print("molecule template list:")
    # for molecule_template in molecule_template_list:
    #     print(molecule_template)
    # print("interface list:")
    # for interface in interface_list:
    #     print(interface)
    # print("interface template list:")
    # for interface_template in interface_template_list:
    #     print(interface_template)

    _build_reactions()

    _rescale_energies()

    if standard_output:
        print("Molecules Template and Reactions Template After Regularization:")
        for molecule_template in molecule_template_list:
            print(molecule_template)
        for reaction_template in reaction_template_list:
            print(reaction_template)

        print("Molecules and Reactions:")
        for molecule in molecule_list:
            print(molecule)
        for reaction in reaction_list:
            print(reaction)

    if show_coarse_grained_structure:
        plot_regularized_structure()

    if save_pymol_script:
        save_regularized_coarse_grained_structure()

    _generate_model_data()

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
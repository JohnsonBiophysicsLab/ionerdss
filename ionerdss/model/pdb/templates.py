"""
templates.py

Generates NERDSS modeling templates (Molecule, Interface, Reaction) from a regularized
coarse-grained protein structure.

This module transforms compact numerical representations of chains and interfaces into
template objects compatible with NERDSS's model definition system. It creates:

- Molecule templates for each unique chain or homologous group
- Interface templates positioned at contact sites
- Optional Reaction templates (e.g., docking events)

This is the final step in the coarse-graining pipeline before model export.

Functions
---------
- build_templates(model: dict) -> dict
    Create Molecule and Interface template objects from regularized geometry.
"""

from ionerdss.model.components import MoleculeTemplate, ReactionTemplate, InterfaceTemplate

def build_templates(model, collapse_homologs=True):
    """
    Convert regularized coarse-grain data into NERDSS-compatible molecule/interface templates.

    Parameters
    ----------
    model : dict
        Output from `regularize_model()`, with:
            'chain_ids', 'chain_labels', 'COMs', 'radii',
            'interfaces', 'interface_coords', 'interface_residues',
            'interface_energies', 'representatives'

    collapse_homologs : bool, optional
        If True, identical molecule types will only be created once (via representatives).

    Returns
    -------
    templates : dict
        {
            'molecules': list of MoleculeTemplate,
            'interfaces': list of InterfaceTemplate,
            'reactions': list of ReactionTemplate
        }
    """
    mol_templates = []
    iface_templates = []
    react_templates = []

    seen_representatives = set()
    chain_map = {}

    for i, rep_idx in enumerate(model["representatives"]):
        if collapse_homologs and rep_idx in seen_representatives:
            chain_map[i] = rep_idx  # reuse template
            continue

        # Build molecule template
        mol = MoleculeTemplate(
            name=f"MOL_{rep_idx}",
            com=model["COMs"][i],
            radius=model["radii"][i],
        )

        iface_list = []
        for j, neighbor_idx in enumerate(model["interfaces"][i]):
            iface = InterfaceTemplate(
                name=f"I{i}_{j}",
                coords=model["interface_coords"][i][j],
                residues=model["interface_residues"][i][j],
                energy=model["interface_energies"][i][j],
            )
            iface_templates.append(iface)
            mol.add_interface(iface)
            iface_list.append(iface)

        mol_templates.append(mol)
        seen_representatives.add(rep_idx)
        chain_map[i] = rep_idx

    # Optionally define reactions (example: 1 interface ↔ 1 interface binding)
    for mol in mol_templates:
        for iface in mol.interfaces:
            reaction = ReactionTemplate(
                reactants=[iface],
                products=[],  # would normally point to bound complex
                ka=1.0,      # placeholder rate
                kb=1.0
            )
            react_templates.append(reaction)

    return {
        "molecules": mol_templates,
        "interfaces": iface_templates,
        "reactions": react_templates
    }

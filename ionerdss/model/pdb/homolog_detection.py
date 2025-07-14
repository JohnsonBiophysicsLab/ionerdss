"""
homolog_detection.py

Identifies homologous protein chains in a biological assembly using either
header annotations (from mmCIF) or pairwise structural alignment. This module is used
in the NERDSS modeling pipeline to group chains into molecule templates and assign
consistent chain labels across symmetry groups.

Primary Goals
-------------
- Group chains that are structurally equivalent into equivalence classes
- Map all chains back to a canonical representative
- Enable template-based coarse-graining and alignment

Approaches
----------
1. **Header-based equivalence**: Uses `_pdbx_struct_assembly_gen` or `_entity_poly.pdbx_strand_id`
   from the mmCIF file to detect symmetric copies.
2. **Structure-based alignment** (fallback): Uses a scoring function (e.g., RMSD) to match chains.

Functions
---------
- identify_homologous_chains(pdb_id: str, structure) -> (chains_map, chains_group)
    Main interface that attempts header parsing, falls back to alignment if needed.
- _find_homologous_chains_by_alignment(chains) -> (chains_map, chains_group)
    Performs pairwise structure-based chain matching.
- _assign_original_chain_ids(chains_group, chains_map, structure)
    Ensures original ID labels are restored for compatibility.
- _validate_chain_mapping(chains_map, structure)
    Optional sanity check to ensure each chain maps to a representative.
"""

from collections import defaultdict
from Bio.PDB import Superimposer, Selection
from Bio.PDB.Polypeptide import is_aa

def identify_homologous_chains(pdb_id, structure):
    """
    Try to identify homologous chains from CIF header. If not found, fall back to structural alignment.

    Parameters
    ----------
    pdb_id : str
        ID or path of the CIF file (for header parsing).
    structure : Bio.PDB.Structure.Structure
        Parsed structure object.

    Returns
    -------
    chains_map : dict
        Mapping of each chain ID to a canonical chain ID (e.g., {'B': 'A', 'C': 'A'}).
    chains_group : list of lists
        List of groups of equivalent chain IDs (e.g., [['A', 'B', 'C'], ['D', 'E']]).
    """
    try:
        # Primary attempt: use entity_poly and assembly_gen fields from mmCIF
        mmcif_dict = structure.header.get('mmcif_dict', {})
        entity_map = defaultdict(list)

        # e.g. {'1': ['A', 'B'], '2': ['C']}
        if '_entity_poly.pdbx_strand_id' in mmcif_dict:
            strand_ids = mmcif_dict['_entity_poly.pdbx_strand_id']
            entity_ids = mmcif_dict['_entity_poly.entity_id']
            for eid, strands in zip(entity_ids, strand_ids):
                for c in strands.replace(' ', '').split(','):
                    entity_map[eid].append(c)
        else:
            raise KeyError("Missing entity_poly.pdbx_strand_id")

        chains_group = list(entity_map.values())
        chains_map = {c: group[0] for group in chains_group for c in group}
        return chains_map, chains_group

    except Exception:
        # Fall back to structural alignment
        print(f"Warning: Falling back to structural alignment for {pdb_id}")
        chains = [chain for chain in structure.get_chains() if any(is_aa(res) for res in chain)]
        return _find_homologous_chains_by_alignment(chains)


def _find_homologous_chains_by_alignment(chains, rmsd_threshold=1.0):
    """
    Fallback method: perform all-by-all structural alignment of chains using C-alpha atoms.

    Parameters
    ----------
    chains : list of Bio.PDB.Chain
        Chains to compare.

    Returns
    -------
    chains_map : dict
        Map each chain to its representative.
    chains_group : list of lists
        List of chain ID groups.
    """
    N = len(chains)
    groups = []
    visited = set()
    chain_ids = [chain.id for chain in chains]

    def extract_calpha(chain):
        return [res['CA'].get_vector() for res in chain if is_aa(res) and 'CA' in res]

    for i in range(N):
        if chain_ids[i] in visited:
            continue
        group = [chain_ids[i]]
        ref_coords = extract_calpha(chains[i])
        for j in range(i + 1, N):
            if chain_ids[j] in visited:
                continue
            mov_coords = extract_calpha(chains[j])
            if len(ref_coords) != len(mov_coords) or len(ref_coords) == 0:
                continue
            si = Superimposer()
            si.set_atoms(ref_coords, mov_coords)
            si.apply(mov_coords)
            if si.rms < rmsd_threshold:
                group.append(chain_ids[j])
                visited.add(chain_ids[j])
        visited.update(group)
        groups.append(group)

    chains_map = {c: group[0] for group in groups for c in group}
    return chains_map, groups


def _assign_original_chain_ids(chains_group, chains_map, structure):
    """
    Updates internal IDs to match original chain order from structure.

    Parameters
    ----------
    chains_group : list of lists
    chains_map : dict
    structure : Bio.PDB.Structure.Structure

    Returns
    -------
    (chains_map, chains_group) with reordered entries.
    """
    all_ids = [chain.id for chain in structure.get_chains()]
    id_order = {c: i for i, c in enumerate(all_ids)}

    # Sort each group and groups by their first appearance
    chains_group = sorted([sorted(g, key=lambda c: id_order.get(c, 1e9)) for g in chains_group],
                          key=lambda g: id_order.get(g[0], 1e9))
    chains_map = {c: g[0] for g in chains_group for c in g}
    return chains_map, chains_group


def _validate_chain_mapping(chains_map, structure):
    """
    Simple validation to ensure that every chain is assigned a representative.

    Parameters
    ----------
    chains_map : dict
    structure : Bio.PDB.Structure.Structure

    Raises
    ------
    ValueError if any chain is missing from map.
    """
    all_ids = {chain.id for chain in structure.get_chains()}
    missing = all_ids - set(chains_map.keys())
    if missing:
        raise ValueError(f"Missing mappings for chains: {missing}")

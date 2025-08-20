"""
repeated_chain_detection.py

Identifies repeated protein chains in a biological assembly using either
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
- identify_repeated_chains(pdb_id: str, structure) -> (chains_map, chains_group)
    Main interface that attempts header parsing, falls back to alignment if needed.
- _find_repeated_chains_by_alignment(chains) -> (chains_map, chains_group)
    Performs pairwise structure-based chain matching.
- _assign_original_chain_ids(chains_group, chains_map, structure)
    Ensures original ID labels are restored for compatibility.
- _validate_chain_mapping(chains_map, structure)
    Optional sanity check to ensure each chain maps to a representative.
"""
import warnings
from collections import defaultdict
from Bio.PDB import Superimposer
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner

from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters


def identify_repeated_chains(pdb_id, structure,
                             params: PDBModelHyperparameters,
                             matching_mode=None):
    """
    Identify repeated chains using mmCIF header, structural alignment, or sequence alignment.

    Parameters
    ----------
    pdb_id : str
        ID or path of the CIF file (used in warnings).
    structure : Bio.PDB.Structure.Structure
        Parsed structure object.
    matching_mode : str, optional
        One of ['default', 'structure', 'sequence'] to control matching strategy.
        Overrides matching_mode in params.
        However, if it is any values other than ['default', 'structure', 'sequence'],
        use params.matching_mode instead.
    params : PDBModelHyperparameters

    Returns
    -------
    chains_map : dict
        Mapping of each chain ID to a canonical chain ID (e.g., {'B': 'A', 'C': 'A'}).
    chains_group : list of lists
        List of groups of equivalent chain IDs (e.g., [['A', 'B', 'C'], ['D', 'E']]).
    """
    chains = [chain for chain in structure.get_chains() if any(is_aa(res)
                                                               for res in chain)]

    # falling back to matching_mode in params if matching_mode override is undefined
    if matching_mode not in {'default', 'structure', 'sequence'}:
        matching_mode = params.matching_mode

    if matching_mode == 'structure':
        return _find_repeated_chains_by_structure_superimposing(chains, params = params)

    if matching_mode == 'sequence':
        return _find_repeated_chains_by_sequence(chains, params = params)

    if matching_mode == 'default':
        mmcif_dict = structure.header.get('mmcif_dict', {})
        entity_map = defaultdict(list)

        if '_entity_poly.pdbx_strand_id' in mmcif_dict:
            strand_ids = mmcif_dict['_entity_poly.pdbx_strand_id']
            entity_ids = mmcif_dict['_entity_poly.entity_id']
            for eid, strands in zip(entity_ids, strand_ids):
                for c in strands.replace(' ', '').split(','):
                    entity_map[eid].append(c)
            chains_group = list(entity_map.values())
            chains_map = {c: group[0] for group in chains_group for c in group}
            return chains_map, chains_group
        else:
            warnings.warn(
                f"{pdb_id}: Missing mmCIF entity_poly info. Falling back to sequence alignment.")
            return _find_repeated_chains_by_sequence(chains, params=params)

    else:
        raise ValueError(
            f"Unknown mode '{matching_mode}'. Valid modes are: 'default', 'structure', 'sequence'.")


def _find_repeated_chains_by_structure_superimposing(chains,
                                                     params : PDBModelHyperparameters):
    """
    Fallback method: perform all-by-all structural alignment of chains using C-alpha atoms.

    Parameters
    ----------
    chains : list of Bio.PDB.Chain
        Chains to compare.
    params : PDBModelHyperparameters

    Returns
    -------
    chains_map : dict
        Map each chain to its representative.
    chains_group : list of lists
        List of chain ID groups.
    """
    num_aa = len(chains)
    groups = []
    visited = set()
    chain_ids = [chain.id for chain in chains]

    def extract_calpha(chain):
        return [res['CA'] for res in chain if is_aa(res) and 'CA' in res]

    for i in range(num_aa):
        if chain_ids[i] in visited:
            continue
        group = [chain_ids[i]]
        ref_coords = extract_calpha(chains[i])
        for j in range(i + 1, num_aa):
            if chain_ids[j] in visited:
                continue
            mov_coords = extract_calpha(chains[j])
            if len(ref_coords) != len(mov_coords) or len(ref_coords) == 0:
                continue
            si = Superimposer()
            si.set_atoms(ref_coords, mov_coords)
            si.apply(mov_coords)
            if si.rms < params.rmsd_threshold:
                group.append(chain_ids[j])
                visited.add(chain_ids[j])
        visited.update(group)
        groups.append(group)

    chains_map = {c: group[0] for group in groups for c in group}
    return chains_map, groups


def _find_repeated_chains_by_sequence(chains,
                                      params : PDBModelHyperparameters
                                      ):
    """
    Group chains by pairwise sequence identity using Bio.Align.PairwiseAligner.

    Parameters
    ----------
    chains : list of Bio.PDB.Chain
    params : PDBModelHyperparameters

    Returns
    -------
    chains_map : dict
    chains_group : list of lists
    """
    ppb = PPBuilder()
    sequences = {}
    chain_ids = [chain.id for chain in chains]

    for chain in chains:
        peptides = ppb.build_peptides(chain)
        sequences[chain.id] = str(
            peptides[0].get_sequence()) if peptides else ""
    if params.custom_aligner is None:
        aligner = PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.5
    else:
        aligner = params.custom_aligner

    chains_groups = []
    visited = set()

    for i in range(len(chains)):
        ci = chain_ids[i]
        if ci in visited:
            continue
        group = [ci]
        seq_i = sequences[ci]

        for j in range(i + 1, len(chains)):
            cj = chain_ids[j]
            if cj in visited:
                continue
            seq_j = sequences[cj]
            if not seq_i or not seq_j:
                continue
            score = aligner.align(seq_i, seq_j).score
            identity = score / max(len(seq_i), len(seq_j))
            if identity >= params.seq_threshold:
                group.append(cj)
                visited.add(cj)
        visited.update(group)
        chains_groups.append(group)

    chains_map = {c: group[0] for group in chains_groups for c in group}

    print(f"{len(chains_groups)} homologous chain groups identified:")
    print(chains_map)

    # sort output
    for group in chains_groups:
        group.sort()
    chains_groups.sort()

    return chains_map, chains_groups

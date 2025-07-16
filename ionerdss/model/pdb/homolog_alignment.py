"""
repeated_chain_alignment.py

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

def regularize_model(raw_data, chains_map=None):
    """
    Canonicalize chain/interface indexing from raw coarse-grained data.

    Parameters
    ----------
    raw_data : dict
        Dictionary returned from `coarse_grain_structure()`:
        {
            'chains': list of Chain objects,
            'COMs': list of Coords,
            'radii': list of float,
            'interfaces': list of neighbor IDs (e.g., ['B', 'C']),
            'interface_coords': list of Coords,
            'interface_residues': list of list of residue indices,
            'interface_energies': list of floats,
        }

    chains_map : dict or None
        Optional mapping from chain ID → canonical chain ID (e.g., 'B' → 'A').
        If provided, chain groups will be collapsed to representative chains.

    Returns
    -------
    dict
        A cleaned and regularized model:
        {
            'chain_ids': [0, 1, ..., N-1],
            'COMs': [Coords, ...],
            'radii': [float, ...],
            'interfaces': [[0, 2], [1], ...],  # neighbor chain indices
            'interface_coords': [[Coords, ...], ...],
            'interface_residues': [[res_ids], ...],
            'interface_energies': [[float, ...], ...],
            'chain_labels': ['A', 'B', ...],   # original chain IDs (for export)
            'representatives': [0, 0, 1, ...], # index of canonical representative per chain
        }
    """
    original_ids = [chain.id for chain in raw_data["chains"]]
    num_chains = len(original_ids)

    id_to_index = {cid: i for i, cid in enumerate(original_ids)}
    chain_labels = original_ids.copy()

    # Optionally collapse to canonical representatives
    if chains_map:
        reps = [chains_map.get(cid, cid) for cid in original_ids]
        unique_reps = sorted(set(reps), key=reps.index)
        rep_to_index = {rep: i for i, rep in enumerate(unique_reps)}
        representatives = [rep_to_index[r] for r in reps]
    else:
        representatives = list(range(num_chains))

    # Build regularized model
    model = {
        "chain_ids": list(range(num_chains)),
        "chain_labels": chain_labels,
        "COMs": raw_data["COMs"],
        "radii": raw_data["radii"],
        "interfaces": [],
        "interface_coords": [],
        "interface_residues": [],
        "interface_energies": [],
        "representatives": representatives
    }

    for i in range(num_chains):
        neighbor_ids = raw_data["interfaces"][i]
        neighbor_idx = [id_to_index[nid] for nid in neighbor_ids]
        model["interfaces"].append(neighbor_idx)
        model["interface_coords"].append(raw_data["interface_coords"][i])
        model["interface_residues"].append(raw_data["interface_residues"][i])
        model["interface_energies"].append(raw_data["interface_energies"][i])

    return model

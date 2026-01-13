"""
Graph-based complex parsing using NetworkX for general topologies.

This module provides an alternative implementation of complex parsing
that uses NetworkX graphs and the graph_based submodule to generate
all subcomplexes. This approach works for any topology (linear, cyclic,
branched, complete, etc.) unlike the original algorithm which was
optimized for linear structures.
"""

import networkx as nx
from collections import defaultdict
from typing import List, Tuple
from .complex import Complex
from .complex_to_graph import complex_to_networkx, networkx_to_complex
from ..graph_based.complexes.subcomplexes import get_unique_fully_connected_subgraphs


def build_pdb_model_graph(pdb_model):
    """
    Build a NetworkX graph from a PDB model's molecules and reactions.
    
    Args:
        pdb_model: PDBModel object containing molecule_list and reaction_list
        
    Returns:
        nx.Graph: Graph where nodes are molecules and edges are binding reactions
        Dict: Mapping from node ID to molecule object
        Dict: Mapping from edge to reaction object
    """
    G = nx.Graph()
    
    # Map molecules to node IDs
    mol_to_id = {mol: i for i, mol in enumerate(pdb_model.molecule_list)}
    id_to_mol = {i: mol for mol, i in mol_to_id.items()}
    
    # Add nodes with molecule template type
    for mol_id, mol in id_to_mol.items():
        mol_type = mol.my_template.name if hasattr(mol.my_template, 'name') else mol.name
        G.add_node(mol_id, type=mol_type, molecule=mol)
    
    # Add edges from reactions
    edge_to_reaction = {}
    for reaction in pdb_model.reaction_list:
        if not reaction.reactants or len(reaction.reactants) != 2:
            continue
        
        mol1, mol2 = reaction.reactants[0][0], reaction.reactants[1][0]
        if mol1 not in mol_to_id or mol2 not in mol_to_id:
            continue
        
        mol1_id = mol_to_id[mol1]
        mol2_id = mol_to_id[mol2]
        
        # Edge type from reaction expression
        edge_type = reaction.my_template.expression if hasattr(reaction.my_template, 'expression') else "binding"
        
        G.add_edge(mol1_id, mol2_id, type=edge_type, reaction=reaction)
        edge_to_reaction[(mol1_id, mol2_id)] = reaction
        edge_to_reaction[(mol2_id, mol1_id)] = reaction  # Bidirectional
    
    return G, id_to_mol, edge_to_reaction


def subgraph_to_complex(subgraph: nx.Graph, id_to_mol: dict, edge_to_reaction: dict) -> Complex:
    """
    Convert a NetworkX subgraph back to a Complex object with actual molecules and reactions.
    
    Args:
        subgraph: NetworkX subgraph from get_unique_fully_connected_subgraphs
        id_to_mol: Mapping from node ID to molecule object
        edge_to_reaction: Mapping from edge tuple to reaction object
        
    Returns:
        Complex: Complex object with proper molecule and reaction references
    """
    complex_obj = Complex()
    
    # Handle empty subgraph
    if len(subgraph.nodes) == 0:
        return complex_obj
    
    # Handle single-molecule complex
    if len(subgraph.nodes) == 1:
        node_id = list(subgraph.nodes)[0]
        mol = id_to_mol[node_id]
        complex_obj.add_interaction(mol, None, None)
        return complex_obj
    
    # Handle multi-molecule complex
    for u, v in subgraph.edges:
        mol_u = id_to_mol[u]
        mol_v = id_to_mol[v]
        
        # Get the reaction object
        reaction = edge_to_reaction.get((u, v))
        if reaction is None:
            # Fallback: try reverse direction
            reaction = edge_to_reaction.get((v, u))
        
        # Add bidirectional interactions
        complex_obj.add_interaction(mol_u, mol_v, reaction)
        complex_obj.add_interaction(mol_v, mol_u, reaction)
    
    return complex_obj


def parse_complexes_from_pdb_model_graphbased(pdb_model, max_complex_size=None) -> List[Complex]:
    """
    Parse all connected complexes from a PDB model using graph-based approach.
    
    This function uses NetworkX and get_unique_fully_connected_subgraphs() to
    generate all possible molecular complexes. This works for any topology including
    linear, cyclic, branched, and complete graphs.
    
    Args:
        pdb_model: The PDBModel object containing molecules and reactions.
        max_complex_size (int, optional): Maximum number of molecules in a complex.
            If None, no limit is applied. Defaults to None.

    Returns:
        List[Complex]: List of all possible complexes.
    """
    # Build graph representation of the PDB model
    G, id_to_mol, edge_to_reaction = build_pdb_model_graph(pdb_model)
    
    # Get all unique fully connected subgraphs
    subgraphs = get_unique_fully_connected_subgraphs(G)
    
    # Filter by max_complex_size if specified
    if max_complex_size is not None:
        subgraphs = [sg for sg in subgraphs if len(sg.nodes) <= max_complex_size]
    
    # Convert subgraphs to Complex objects
    complex_list = []
    for subgraph in subgraphs:
        complex_obj = subgraph_to_complex(subgraph, id_to_mol, edge_to_reaction)
        complex_list.append(complex_obj)
    
    return complex_list

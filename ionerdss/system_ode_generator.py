"""
System-compatible ODE model generator using graph_based functions.

This module generates ODE models from the System architecture using the actual
graph_based functions for proper species and reaction generation.
"""

import networkx as nx
from typing import List, Tuple
from ionerdss.model.components.system import System
from ionerdss.model.complex import ComplexReactionSystem
from ionerdss.model.complex_to_graph import generate_complex_name_from_graph


def generate_ode_model_from_system(system: System, max_complex_size: int = None, pdb_model=None, coarse_grainer=None) -> Tuple[List, ComplexReactionSystem]:
    """
    Generate ODE model from a System object using graph_based functions.
    
    This function:
    1. Builds the full assembly graph from the PDB model
    2. Uses get_unique_fully_connected_subgraphs to get all species
    3. Uses find_all_dimer_reactions and find_all_transformable_subgraph_pairs for reactions
    
    Args:
        system: System object containing molecule_types and interface_types registries
        max_complex_size: Maximum number of molecules in a complex (default: 20)
        pdb_model: PDB model object (optional, for compatibility)
        coarse_grainer: CoarseGrainer object with coarse-grained model data
    
    Returns:
        Tuple of (complex_names, reaction_system) where complex_names are string identifiers
        and reaction_system contains reactions between complexes
    """
    # Set default max size
    if max_complex_size is None:
        max_complex_size = 12
    
    if coarse_grainer is None:
        raise ValueError("coarse_grainer must be provided to generate ODE")
    
    # Step 1: Build the full assembly graph using build_simple_graph
    from ionerdss.model.graph_based.complexes.graphize import build_simple_graph
    
    # Build cg_model from coarse_grainer data
    chains_data = coarse_grainer.get_coarse_grained_chains()
    interfaces_data = coarse_grainer.get_interfaces()
    
    # Build the cg_model dict expected by build_simple_graph
    chains = list(chains_data.keys())
    
    # Build interfaces list: for each chain, list of partner chain IDs
    interfaces = [[] for _ in chains]
    chain_to_idx = {cid: i for i, cid in enumerate(chains)}
    
    for iface in interfaces_data:
        chain_i = iface.chain_i
        chain_j = iface.chain_j
        if chain_i in chain_to_idx and chain_j in chain_to_idx:
            idx_i = chain_to_idx[chain_i]
            idx_j = chain_to_idx[chain_j]
            # Add bidirectional connections
            if chain_j not in interfaces[idx_i]:
                interfaces[idx_i].append(chain_j)
            if chain_i not in interfaces[idx_j]:
                interfaces[idx_j].append(chain_i)
    
    cg_model = {
        'chains': chains,
        'interfaces': interfaces
    }
    
    G_full = build_simple_graph(cg_model)
    
    # Check if full assembly exceeds max_complex_size
    if len(G_full.nodes) > max_complex_size:
        raise ValueError(
            f"Assembly has {len(G_full.nodes)} molecules, exceeding max_complex_size_ode ({max_complex_size}). "
            f"Skipping ODE generation. Increase max_complex_size_ode parameter if you want to force calculation of ODE for this system, but be aware that this may take a long time."
        )
    
    # Step 2: Generate all unique fully connected subgraphs (species)
    from ionerdss.model.graph_based.complexes.subcomplexes import get_unique_fully_connected_subgraphs
    
    all_subgraphs_sets = get_unique_fully_connected_subgraphs(G_full)
    
    # Convert frozensets to NetworkX graphs
    # get_unique_fully_connected_subgraphs returns frozensets of node IDs, not graph objects
    import networkx as nx
    all_subgraphs = []
    for node_set in all_subgraphs_sets:
        # Create subgraph from the node set
        subgraph = G_full.subgraph(node_set).copy()
        all_subgraphs.append(subgraph)
    
    # Filter by max_complex_size
    subgraphs = [sg for sg in all_subgraphs if len(sg.nodes) <= max_complex_size]
    
    # Generate names for each subgraph using graph-based naming
    complex_names = []
    for subgraph in subgraphs:
        name = generate_complex_name_from_graph(subgraph, use_hash=True)
        complex_names.append(name)
    
    # Step 3: Build reaction system using graph_based functions
    from ionerdss.model.graph_based.reactions import find_all_dimer_reactions, find_all_transformable_subgraph_pairs
    
    reaction_system = ComplexReactionSystem()
    
    # Get dimer reactions (A + B -> AB)
    dimer_reactions = find_all_dimer_reactions(subgraphs, use_multiprocessing=False)
    
    # Get transformation reactions (bond formation/breaking)
    transformation_pairs = find_all_transformable_subgraph_pairs(G_full, subgraphs=subgraphs)
    
    # Convert graph reactions to reaction strings
    # Map subgraphs to their names using graph structure (node sets) instead of object IDs
    # because find_all_dimer_reactions creates new graph objects
    import networkx as nx
    
    subgraph_nodeset_to_name = {}
    for i, sg in enumerate(subgraphs):
        node_set = frozenset(sg.nodes())
        subgraph_nodeset_to_name[node_set] = complex_names[i]
    
    reaction_idx = 0
    
    # Process dimer reactions
    for reaction in dimer_reactions:
        # reaction format from find_all_dimer_reactions: (set1, set2, product_set) - sets not graphs!
        if len(reaction) >= 3:
            set1, set2, set_product = reaction[0], reaction[1], reaction[2]
            
            # Convert to frozensets for lookup
            nodeset1 = frozenset(set1)
            nodeset2 = frozenset(set2)
            nodeset_product = frozenset(set_product)
            
            name1 = subgraph_nodeset_to_name.get(nodeset1)
            name2 = subgraph_nodeset_to_name.get(nodeset2)
            name_product = subgraph_nodeset_to_name.get(nodeset_product)
            
            if name1 and name2 and name_product:
                rate_const_name = f"k_on_{reaction_idx}"
                reaction_expr = f"{name1} + {name2} -> {name_product}, {rate_const_name}"
                
                class SimpleReaction:
                    def __init__(self, expression, rate=1.0, rate_name=None):
                        self.expression = expression
                        self.rate = rate
                        self.rate_name = rate_name if rate_name else "k_on"
                
                rxn = SimpleReaction(reaction_expr, rate=1.0, rate_name=rate_const_name)
                reaction_system.reactions.append(rxn)
                reaction_idx += 1
    
    # Process transformation reactions
    for G1, G2, direction, edges_changed in transformation_pairs:
        # Find names using node sets
        nodeset1 = frozenset(G1.nodes()) if hasattr(G1, 'nodes') else frozenset(G1)
        nodeset2 = frozenset(G2.nodes()) if hasattr(G2,  'nodes') else frozenset(G2)
        name1 = subgraph_nodeset_to_name.get(nodeset1)
        name2 = subgraph_nodeset_to_name.get(nodeset2)
        
        if name1 and name2:
            rate_const_name = f"k_trans_{reaction_idx}"
            
            if direction == "forming":
                # Bond formation: G1 -> G2
                reaction_expr = f"{name1} -> {name2}, {rate_const_name}"
            else:
                # Bond breaking: G2 -> G1  
                reaction_expr = f"{name2} -> {name1}, {rate_const_name}"
            
            class SimpleReaction:
                def __init__(self, expression, rate=1.0, rate_name=None):
                    self.expression = expression
                    self.rate = rate
                    self.rate_name = rate_name if rate_name else "k_trans"
            
            rxn = SimpleReaction(reaction_expr, rate=1.0, rate_name=rate_const_name)
            reaction_system.reactions.append(rxn)
            reaction_idx += 1
    
    return complex_names, reaction_system

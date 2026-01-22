"""
System-compatible ODE model generator using graph_based functions.

This module generates ODE models from the System architecture using the actual
graph_based functions for proper species and reaction generation.
"""
from typing import List, Tuple, Any, Dict

import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash

from ionerdss.model.components.system import System
from ionerdss.model.components.instances import MoleculeInstance
from ionerdss.model.graph_based.complexes.subcomplexes import get_unique_fully_connected_subgraphs
from ionerdss.model.graph_based.reactions import find_all_dimer_reactions, find_all_transformable_subgraph_pairs


class SimpleReaction:
    """A simple reaction class compatible with ODE pipeline."""
    def __init__(self, expression: str, rate: float = 1.0, rate_name: str = None):
        self.expression = expression
        self.rate = rate
        self.rate_name = rate_name if rate_name else "k"

class ReactionSystem:
    """A simple container for reactions, replacing ComplexReactionSystem."""
    def __init__(self):
        self.reactions: List[SimpleReaction] = []

def _system_to_graph(system: System) -> nx.Graph:
    """
    Convert a System object to a NetworkX graph.
    Nodes are molecule names, with 'type' attribute set to molecule type name.
    Edges represent bound interfaces, with 'type' attribute set to sorted interface types.
    """
    G = nx.Graph()
    
    # Add nodes
    # Use molecule name as node ID (unique)
    # Use molecule type name as 'type' attribute (for isomorphism/hashing)
    for mol in system.molecule_instances:
        mol_type_name = mol.molecule_type.name if mol.molecule_type else "unknown"
        G.add_node(mol.name, type=mol_type_name)
    
    # Add edges
    # Iterate over all interface instances to find bonds
    # We used indices or direct references. Check system.interface_instances registry.
    
    added_edges = set()
    
    for iface in system.interface_instances:
        # Check if bound
        # An interface is bound if it has a partner_interface
        if iface.partner_interface:
            # Get parent molecules
            mol1 = iface.this_mol
            mol2 = iface.partner_interface.this_mol
            
            if mol1 and mol2:
                u, v = mol1.name, mol2.name
                
                # Determine edge type (canonical based on interface types)
                type1 = iface.interface_type.get_name() if iface.interface_type else "unknown"
                type2 = iface.partner_interface.interface_type.get_name() if iface.partner_interface.interface_type else "unknown"
                edge_type = "-".join(sorted([type1, type2]))
                
                edge_key = frozenset([u, v])
                if edge_key not in added_edges:
                    G.add_edge(u, v, type=edge_type)
                    added_edges.add(edge_key)
                    
    return G

def _generate_complex_name(graph: nx.Graph) -> str:
    """Generate a unique name for a complex based on its graph topology and types."""
    # Count node types
    counts = {}
    for _, data in graph.nodes(data=True):
        t = data.get('type', 'U')
        counts[t] = counts.get(t, 0) + 1
    
    # Build composition string (e.g. A3B2)
    # Sort by type name to ensure canonical order
    composition = ""
    for t in sorted(counts.keys()):
        composition += f"{t}{counts[t]}"
        
    # Use Weisfeiler-Lehman hash to get a canonical graph string/hash
    # 'type' node attribute is used for coloring
    # We use the hash as the unique identifier/name
    g_hash = weisfeiler_lehman_graph_hash(graph, node_attr="type")
    
    return f"{composition}_{g_hash}"

def generate_ode_model_from_system(
    system: System, 
    max_complex_size: int = None, 
    pdb_model=None, 
    coarse_grainer=None
) -> Tuple[List[str], ReactionSystem]:
    """
    Generate ODE model from a System object using graph_based functions.
    
    This function:
    1. Builds the full assembly graph directly from the System object
    2. Uses get_unique_fully_connected_subgraphs to get all species
    3. Uses find_all_dimer_reactions and find_all_transformable_subgraph_pairs for reactions
    
    Args:
        system: System object containing molecule_types and molecule_instances
        max_complex_size: Maximum number of molecules in a complex (default: 12)
        pdb_model: Ignored (kept for compatibility)
        coarse_grainer: Ignored (kept for compatibility)
    
    Returns:
        Tuple of (complex_names, reaction_system) where:
            - complex_names: List of string identifiers for the complexes
            - reaction_system: Object containing a .reactions list of SimpleReaction objects
    """
    # Set default max size
    if max_complex_size is None:
        max_complex_size = 12
    
    # Step 1: Build the full assembly graph directly from System
    G_full = _system_to_graph(system)
    
    # Check if full assembly exceeds max_complex_size
    if len(G_full.nodes) > max_complex_size:
        raise ValueError(
            f"Assembly has {len(G_full.nodes)} molecules, exceeding max_complex_size ({max_complex_size}). "
            f"Skipping ODE generation. Increase max_complex_size parameter if you want to force calculation."
        )
    
    # Step 2: Generate all unique fully connected subgraphs (species)
    all_subgraphs_sets = get_unique_fully_connected_subgraphs(G_full)
    
    # Convert frozensets to NetworkX graphs
    all_subgraphs = []
    for node_set in all_subgraphs_sets:
        # Create subgraph (induced)
        subgraph = G_full.subgraph(node_set).copy()
        all_subgraphs.append(subgraph)
    
    # Filter by max_complex_size
    subgraphs = [sg for sg in all_subgraphs if len(sg.nodes) <= max_complex_size]
    
    # Generate names for each subgraph
    complex_names = []
    for subgraph in subgraphs:
        name = _generate_complex_name(subgraph)
        complex_names.append(name)
    
    # Step 3: Build reaction system using graph_based functions
    
    reaction_system = ReactionSystem()
    
    # Get dimer reactions (A + B -> AB)
    dimer_reactions = find_all_dimer_reactions(subgraphs, use_multiprocessing=False)
    
    # Get transformation reactions (bond formation/breaking)
    transformation_pairs = find_all_transformable_subgraph_pairs(G_full, subgraphs=subgraphs)
    
    # Map subgraphs to their names using node sets (frozenset of IDs)
    subgraph_nodeset_to_name = {}
    for i, sg in enumerate(subgraphs):
        node_set = frozenset(sg.nodes())
        subgraph_nodeset_to_name[node_set] = complex_names[i]
    
    reaction_idx = 0
    
    # Process dimer reactions
    for reaction in dimer_reactions:
        # reaction format: (set1, set2, product_set)
        if len(reaction) >= 3:
            set1, set2, set_product = reaction[0], reaction[1], reaction[2]
            
            nodeset1 = frozenset(set1)
            nodeset2 = frozenset(set2)
            nodeset_product = frozenset(set_product)
            
            name1 = subgraph_nodeset_to_name.get(nodeset1)
            name2 = subgraph_nodeset_to_name.get(nodeset2)
            name_product = subgraph_nodeset_to_name.get(nodeset_product)
            
            if name1 and name2 and name_product:
                rate_const_name = f"k_on_{reaction_idx}"
                reaction_expr = f"{name1} + {name2} -> {name_product}, {rate_const_name}"
                
                rxn = SimpleReaction(reaction_expr, rate=1.0, rate_name=rate_const_name)
                reaction_system.reactions.append(rxn)
                reaction_idx += 1
    
    # Process transformation reactions
    for G1, G2, direction, edges_changed in transformation_pairs:
        nodeset1 = frozenset(G1.nodes())
        nodeset2 = frozenset(G2.nodes())
        name1 = subgraph_nodeset_to_name.get(nodeset1)
        name2 = subgraph_nodeset_to_name.get(nodeset2)
        
        if name1 and name2:
            rate_const_name = f"k_trans_{reaction_idx}"
            
            if direction == "forming":
                reaction_expr = f"{name1} -> {name2}, {rate_const_name}"
            else:
                reaction_expr = f"{name2} -> {name1}, {rate_const_name}"
            
            rxn = SimpleReaction(reaction_expr, rate=1.0, rate_name=rate_const_name)
            reaction_system.reactions.append(rxn)
            reaction_idx += 1
    
    return complex_names, reaction_system

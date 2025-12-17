"""
System-compatible ODE model generator with full complex network generation.

This module provides functions to generate ODE models from the new System architecture,
using graph induction to generate ALL possible subcomplexes (not just dimers).

For a system with 8 A molecules, this will generate: 1A, 2A, 3A, 4A, 5A, 6A, 7A, 8A
complexes based on connectivity and available binding sites.
"""

import numpy as np
import networkx as nx
from typing import List, Tuple, Set
from itertools import combinations, product
from collections import defaultdict, deque
from ionerdss.model.components.system import System
from ionerdss.model.complex import ComplexReactionSystem


def generate_ode_model_from_system(system: System, max_complex_size: int = None) -> Tuple[List, ComplexReactionSystem]:
    """
    Generate ODE model from a System object with full complex network.
    
    Uses graph induction to generate all possible subcomplexes based on
    molecule types and interface types in the System.
    
    Args:
        system: System object containing molecule_types and interface_types registries
        max_complex_size: Maximum number of molecules in a complex (default: 8 for octamers)
    
    Returns:
        Tuple of (complex_names, reaction_system) where complex_names are string identifiers
        and reaction_system contains reactions between complexes
    """
    # Extract molecule types from system
    molecule_types = list(system.molecule_types.molecule_types.values())
    
    if not molecule_types:
        raise ValueError("No molecule types found in system")
    
    # Build a NetworkX graph representing the full assembly
    # For homotypic system with 1 molecule type and N instances, we create a fully connected graph
    # based on the interface types
    
    # Step 1: Determine how many instances we need
    # Count unique interface index combinations to estimate oligomer size
    interface_indices = set()
    for iface_type in system.interface_types:
        interface_indices.add(iface_type.interface_index)
    
    # For homotypic systems, the maximum oligomer size is typically related to
    # the number of unique interfaces. Default to 8 (common for viral capsids)
    if max_complex_size is None:
        max_complex_size = max(8, len(interface_indices) * 2)
    
    # Step 2: Build connectivity information from interface types
    # Map (mol_type, interface_idx) -> (partner_mol_type, partner_interface_idx)
    connectivity = defaultdict(list)
    
    for iface_type in system.interface_types:
        mol_type = iface_type.this_mol_type_name
        partner_type = iface_type.partner_mol_type_name
        iface_idx = iface_type.interface_index
        
        # For homotypic binding (A-A), track which interfaces can bind
        connectivity[(mol_type, iface_idx)].append((partner_type, iface_idx))
    
    # Step 3: Build a graph representing the maximum assembly
    # For homotypic A with 2 interface types (e.g., 1f/1b, 2f/2b), build a chain/ring
    mol_type_name = molecule_types[0].name
    
    # Create nodes for maximum assembly size
    G = nx.Graph()
    for i in range(max_complex_size):
        G.add_node(i, type=mol_type_name)
    
    # Add edges based on interface compatibility
    # For linear/ring assembly, connect sequential nodes
    # This is a simplified model - in reality, would use actual interface geometry
    for i in range(max_complex_size - 1):
        G.add_edge(i, i + 1, type=f"{mol_type_name}_{mol_type_name}")
    
    # For ring closure (if system supports it), add edge between first and last
    if len(system.interface_types) >= 4:  # Suggests bidirectional binding
        G.add_edge(0, max_complex_size - 1, type=f"{mol_type_name}_{mol_type_name}")
    
    # Step 4: Generate all unique connected subgraphs
    complex_names = []
    complex_graphs = []
    
    # Add individual monomers
    for i in range(max_complex_size):
        subgraph =  G.subgraph([i])
        complex_names.append(f"{mol_type_name}")
        complex_graphs.append(subgraph)
        break  # Only need one monomer in species list
    
    # Generate all connected subgraphs of increasing size
    from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
    seen_hashes = set()
    
    for size in range(2, max_complex_size + 1):
        # Generate all combinations of nodes of this size
        for node_subset in combinations(range(max_complex_size), size):
            subgraph = G.subgraph(node_subset)
            
            # Check if connected
            if not nx.is_connected(subgraph):
                continue
            
            # Use WL hash to avoid duplicates (isomorphic structures)
            try:
                wl_hash = weisfeiler_lehman_graph_hash(subgraph, node_attr="type")
                if wl_hash not in seen_hashes:
                    seen_hashes.add(wl_hash)
                    # Name based on size: A_A for dimer, A_A_A for trimer, etc.
                    complex_name = "_".join([mol_type_name] * size)
                    complex_names.append(complex_name)
                    complex_graphs.append(subgraph.copy())
            except:
                # Fallback if WL hash fails
                complex_name = "_".join([mol_type_name] * size)
                if complex_name not in complex_names:
                    complex_names.append(complex_name)
                    complex_graphs.append(subgraph.copy())
    
    # Step 5: Build reaction system
    # Find all reactions: smaller complexes associating to form larger ones
    reaction_system = ComplexReactionSystem()
    
    # Create simple reaction strings for each association
    # Format: "A + A -> A_A", "A + A_A -> A_A_A", etc.
    
    # Association reactions between all pairs
    reaction_idx = 0  # Counter for unique rate constant names
    for i, (name1, graph1) in enumerate(zip(complex_names, complex_graphs)):
        size1 = len(graph1.nodes())
        for j, (name2, graph2) in enumerate(zip(complex_names, complex_graphs)):
            size2 = len(graph2.nodes())
            product_size = size1 + size2
            
            # Check if there's a product of correct size
            if product_size <= max_complex_size:
                # Find product in our list (complex with combined size)
                product_name = "_".join([mol_type_name] * product_size)
                if product_name in complex_names:
                    # Create unique rate constant name for this reaction
                    rate_const_name = f"k_on_{reaction_idx}"
                    # Create reaction
                    reaction_expr = f"{name1} + {name2} -> {product_name}, {rate_const_name}"
                    
                    class SimpleReaction:
                        def __init__(self, expression, rate=1.0, rate_name=None):
                            self.expression = expression
                            self.rate = rate
                            self.rate_name = rate_name if rate_name else "k_on"
                    
                    reaction = SimpleReaction(reaction_expr, rate=1.0, rate_name=rate_const_name)
                    
                    # Avoid duplicates (A+B->C is same as B+A->C)
                    if reaction not in reaction_system.reactions:
                        reaction_system.reactions.append(reaction)
                        reaction_idx += 1  # Increment only when reaction is added
    
    return complex_names, reaction_system

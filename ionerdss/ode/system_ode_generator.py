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
                
                # Get energy (average of both sides or just take one, assuming consistency)
                # Default to 0.0 if not specified (though usually negative for binding)
                e1 = iface.interface_type.energy if iface.interface_type and iface.interface_type.energy is not None else 0.0
                e2 = iface.partner_interface.interface_type.energy if iface.partner_interface.interface_type and iface.partner_interface.interface_type.energy is not None else 0.0
                # Use e1 if both present? Or average? 
                # In standard definition, InterfaceType stores deltaG.
                # Assuming symmetric definition or consistent typing.
                bond_energy = e1 
                
                edge_key = frozenset([u, v])
                if edge_key not in added_edges:
                    G.add_edge(u, v, type=edge_type, energy=bond_energy)
                    added_edges.add(edge_key)
                    
    return G

def _calculate_reaction_rates(
    reactants: List[nx.Graph], 
    products: List[nx.Graph], 
    default_kon: float = 120.0
) -> Tuple[float, float]:
    """
    Calculate forward and reverse rates for a reaction based on bond energetics.
    Uses 'energy' attribute on edges (representing deltaG/kT or similar).
    
    Formula:
    - delta_G_reaction = sum(E_bonds_products) - sum(E_bonds_reactants)
    - Keq = exp(-delta_G_reaction)
    - If forming bonds: k_fwd = max(kon_i) approx default_kon. k_rev = k_fwd / Keq
    - If breaking bonds: k_rev = max(kon_i) approx default_kon. k_fwd = k_rev * Keq
      (Effectively anchoring the association rate)
    """
    import numpy as np
    
    # Collect edge energies
    # Collect edge energies with unit normalization
    # Constants for normalization
    R_kJ = 0.008314
    T = 298.0
    RT = R_kJ * T
    
    def _sum_energy(graphs):
        total_E_RT = 0.0
        edges_seen = set()
        for G in graphs:
            for u, v, d in G.edges(data=True):
                key = frozenset([u, v])
                if key not in edges_seen:
                    # Get raw energy (stored in edge)
                    # Could be -1.0 (flag) or kJ/mol
                    raw_E = d.get('energy', 0.0)
                    
                    if raw_E == -1.0:
                        # Default strong binding: -16 RT
                        norm_E = -16.0
                    elif raw_E == 0.0:
                         # Assume 0 energy
                         norm_E = 0.0
                    else:
                        # Explicit energy in kJ/mol (e.g. from ProAffinity)
                        # Normalize by RT
                        norm_E = raw_E / RT
                        
                    total_E_RT += norm_E
                    edges_seen.add(key)
        return total_E_RT, edges_seen

    E_reactants, r_edges = _sum_energy(reactants)
    E_products, p_edges = _sum_energy(products)
    
    # Delta G of reaction in RT units
    # G_state = sum(E_bonds_RT). 
    # Delta G = G_final - G_initial
    # Typically E_bonds is negative for stability.
    delta_G = E_products - E_reactants
    
    # Unit Conversion Constants
    # 1. Convert kon from nm^3/us to uM^-1 s^-1
    #    Factor ~0.6022
    #    Derivation: 1 nm^3/us = 1e-15 cm^3/s. 
    #    N_A * 1e-15 * 1e-3 (to L) ... wait.
    #    Directly: 120 nm^3/us -> ~72 uM^-1 s^-1.
    #    Precise factor: 0.602214
    CONV_NM3_US_TO_UM_S = 0.602214
    
    # 2. Standard Concentration activity correction for uM units
    #    C0 = 1 M = 10^6 uM
    C0_uM = 1.0e6 
    
    k_fwd_val = default_kon * CONV_NM3_US_TO_UM_S
    
    # Determine reaction molecularity change (delta n)
    # 2 reactants -> 1 product : delta_n = -1 (Association)
    # 1 reactant -> 2 products : delta_n = +1 (Dissociation)
    # 1 -> 1 : delta_n = 0 (Isomerization)
    delta_n = len(products) - len(reactants)
    
    # Calculate Equilibrium Constant Kc in units of uM^delta_n
    # K_eq (dimensionless) = exp(-delta_G/RT) [activity based]
    # Kc = K_eq_activity * (C0_uM)^delta_n
    
    Keq_activity = np.exp(-delta_G) # delta_G assumed in RT units
    Kc = Keq_activity * (C0_uM ** delta_n)
    
    # Assign rates
    if delta_n < 0:
        # Association dominant (forming bonds)
        # Anchor forward rate (k_on in uM^-1 s^-1)
        k_fwd = k_fwd_val
        k_rev = k_fwd / Kc
    elif delta_n > 0:
        # Dissociation dominant (breaking bonds)
        # Anchor reverse association (if it were occurring)
        k_rev = k_fwd_val
        k_fwd = k_rev * Kc
    else:
        # Isomerization (1->1)
        k_rev = k_fwd_val # ~72 s^-1
        k_fwd = k_rev * Kc
        
    return k_fwd, k_rev

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
    g_hash = weisfeiler_lehman_graph_hash(graph, node_attr="type", edge_attr="type")
    
    return f"{composition}_{g_hash}"

def generate_ode_model_from_system(
    system: System, 
    max_complex_size: int = None, 
    pdb_model=None, 
    coarse_grainer=None,
    include_transformation_reactions: bool = False
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
        include_transformation_reactions: Whether to include bond rearrangement reactions (default: False)
    
    Returns:
        Tuple of (complex_names, reaction_system) where:
            - complex_names: List of string identifiers for the complexes
            - reaction_system: Object containing a .reactions list of SimpleReaction objects
    """
    # Set default max size
    if max_complex_size is None:
        max_complex_size = 12
        
    # Default kon used when no other rate information is available
    DEFAULT_KON = 120.0
    
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
    
    reaction_idx = 0
    
    # Process dimer reactions
    for reaction in dimer_reactions:
        # reaction format: (set1, set2, product_set)
        if len(reaction) >= 3:
            set1, set2, set_product = reaction[0], reaction[1], reaction[2]
            
            # Generate names on-the-fly because set1/set2 might not be the exact 
            # representative node sets stored in 'subgraphs' list, but are isomorphic.
            sub1 = G_full.subgraph(set1)
            sub2 = G_full.subgraph(set2)
            sub_prod = G_full.subgraph(set_product)
            
            name1 = _generate_complex_name(sub1)
            name2 = _generate_complex_name(sub2)
            name_product = _generate_complex_name(sub_prod)
            
            # Note: We don't filter by existence in 'complex_names' because
            # if the parts are valid induced subgraphs, they represent valid species.
            # However, we only care about reactions where reactants/products are within max_complex_size.
            # find_all_dimer_reactions usually ensures this if input subgraphs are filtered?
            # Actually, set_product is formed by union. If it exceeds max_size, we should skip?
            # 'subgraphs' input was filtered. 'set_product' is one of the 'subgraphs' (passed as specie).
            # So product is definitely within size.
            # Reactants are smaller. So they are within size.
            
            if name1 and name2 and name_product:
                # Calculate Rates
                k_on, k_off = _calculate_reaction_rates(
                    reactants=[sub1, sub2],
                    products=[sub_prod],
                    default_kon=DEFAULT_KON
                )
                
                # Forward Reaction (Association)
                rate_const_name_fwd = f"k_on_{reaction_idx}"
                reaction_expr_fwd = f"{name1} + {name2} -> {name_product}, {rate_const_name_fwd}"
                rxn_fwd = SimpleReaction(reaction_expr_fwd, rate=k_on, rate_name=rate_const_name_fwd)
                reaction_system.reactions.append(rxn_fwd)
                
                # Reverse Reaction (Dissociation)
                rate_const_name_rev = f"k_off_{reaction_idx}"
                reaction_expr_rev = f"{name_product} -> {name1} + {name2}, {rate_const_name_rev}"
                rxn_rev = SimpleReaction(reaction_expr_rev, rate=k_off, rate_name=rate_const_name_rev)
                reaction_system.reactions.append(rxn_rev)
                
                reaction_idx += 1
    
    # Process transformation reactions (optional)
    if include_transformation_reactions:
        transformation_pairs = find_all_transformable_subgraph_pairs(G_full, subgraphs=subgraphs)
        
        for G1, G2, direction, edges_changed in transformation_pairs:
            # G1 and G2 are graphs from the subgraphs list, so we can name them directly
            name1 = _generate_complex_name(G1)
            name2 = _generate_complex_name(G2)
            
            if name1 and name2:
                # Calculate Rates
                k_fwd, k_rev = _calculate_reaction_rates(
                    reactants=[G1],
                    products=[G2],
                    default_kon=DEFAULT_KON
                )
                
                # Forward
                rate_const_name_fwd = f"k_trans_fwd_{reaction_idx}"
                reaction_expr_fwd = f"{name1} -> {name2}, {rate_const_name_fwd}"
                rxn_fwd = SimpleReaction(reaction_expr_fwd, rate=k_fwd, rate_name=rate_const_name_fwd)
                reaction_system.reactions.append(rxn_fwd)
                
                # Reverse
                rate_const_name_rev = f"k_trans_rev_{reaction_idx}"
                reaction_expr_rev = f"{name2} -> {name1}, {rate_const_name_rev}"
                rxn_rev = SimpleReaction(reaction_expr_rev, rate=k_rev, rate_name=rate_const_name_rev)
                reaction_system.reactions.append(rxn_rev)
                
                reaction_idx += 1
    
    return complex_names, reaction_system

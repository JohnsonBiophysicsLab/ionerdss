"""
System-compatible ODE model generator using graph_based functions.

This module generates ODE models from the System architecture using the actual
graph_based functions for proper species and reaction generation.
"""
from typing import List, Tuple, Any, Dict
import numpy as np

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
                bond_energy = e1 
                
                # Dynamic sigma (bond length): distance between COMs
                # Assuming interfaces meet face-to-face perfectly, sigma is sum of magnitudes of local vectors.
                len1 = np.linalg.norm(iface.interface_type.local_coord) if iface.interface_type else 0.0
                len2 = np.linalg.norm(iface.partner_interface.interface_type.local_coord) if iface.partner_interface.interface_type else 0.0
                bond_sigma = float(len1 + len2)
                
                # Dynamic D_tot: sum of translational diffusion constants of the two molecules forming the initial bond
                dt1 = mol1.molecule_type.D_t_nm2_us if mol1.molecule_type else 0.0
                dt2 = mol2.molecule_type.D_t_nm2_us if mol2.molecule_type else 0.0
                bond_D_tot = float(dt1 + dt2)
                
                edge_key = frozenset([u, v])
                if edge_key not in added_edges:
                    G.add_edge(u, v, type=edge_type, energy=bond_energy, sigma=bond_sigma, D_tot=bond_D_tot)
                    added_edges.add(edge_key)
                    
    return G

def _calculate_macroscopic_kon(
    ka: float,
    sigma: float,
    D_tot: float,
) -> float:
    """
    Compute the macroscopic (diffusion-influenced) association rate constant using the
    Smoluchowski-Collins-Kimball (SCK) formula:

        k_on = (1/k_a + 1/(4*pi*sigma*D_tot))^{-1}

    All quantities must be in consistent units **before** calling this helper.
    Recommended unit system: nm^3/us for rates/diffusion products.

    Args:
        ka    : Microscopic (activation-limited) on-rate in nm^3/us.
        sigma : Binding (contact) radius in nm.
        D_tot : Total relative diffusion coefficient D1 + D2 in nm^2/us.

    Returns:
        Macroscopic on-rate in nm^3/us.
    """
    import math
    k_diff = 4.0 * math.pi * sigma * D_tot   # diffusion-limited rate, nm^3/us
    return 1.0 / (1.0 / ka + 1.0 / k_diff)


def _calculate_reaction_rates(
    reactants: List[nx.Graph],
    products: List[nx.Graph],
    default_ka: float = 120.0,
    m_fwd: int = 1,
    m_rev: int = 1,
) -> Tuple[float, float]:
    """
    Calculate forward and reverse rates for a reaction.

    The macroscopic on-rate uses the SCK formula:

        k_on = (1/k_a + 1/(4*pi*sigma*D_tot))^{-1}   [nm^3/us -> uM^-1 s^-1]

    Keq is derived from bond energetics (sum of edge deltaG in RT units):

        Keq_activity = exp(-delta_G_RT)
        Kc = Keq_activity * C0^(delta_n)     [C0 = 1 M = 10^6 uM]

    The macroscopic off-rate is then:

        k_off = k_on_macro / Kc      [s^-1]

    so that detailed balance is maintained: Keq = k_on / k_off = Kc.

    Args:
        reactants  : List of NetworkX graphs representing reactant complexes.
        products   : List of NetworkX graphs representing product complexes.
        default_ka : Microscopic activation-limited on-rate in nm^3/us. Default: 120.0.
        sigma      : Binding contact radius in nm. Default: 1.0.
        D_tot      : Total relative diffusion coefficient D1+D2 in nm^2/us. Default: 10.0.

    Returns:
        Tuple (k_fwd, k_rev) in (uM^-1 s^-1, s^-1) for association,
        or (s^-1, s^-1) for isomerization / dissociation.
    """
    import numpy as np

    # ------------------------------------------------------------------
    # Unit conversion: nm^3/us -> uM^-1 s^-1
    #   1 nm^3/us = N_A * 1e-24 L * 1e6 /s = 0.602214 uM^-1 s^-1
    # ------------------------------------------------------------------
    CONV_NM3_US_TO_UM_S = 0.602214

    # ------------------------------------------------------------------
    # Bond-energy and dynamic SCK properties from the newly formed bond(s)
    # ------------------------------------------------------------------
    R_kJ = 0.008314
    T = 298.0
    RT = R_kJ * T

    def _get_reaction_edges(reactants, products):
        r_edges = set()
        for G in reactants:
            for u, v in G.edges():
                r_edges.add(frozenset([u, v]))
        p_edges = set()
        for G in products:
            for u, v in G.edges():
                p_edges.add(frozenset([u, v]))
        return p_edges - r_edges

    new_edges = _get_reaction_edges(reactants, products)
    broken_edges = _get_reaction_edges(products, reactants)
    delta_n = len(products) - len(reactants)
    C0_uM = 1.0e6

    def _get_edge_properties(edges, source_graphs):
        types = []
        energies = []
        kons = []
        for edge in edges:
            u, v = tuple(edge)
            rxn_sigma, rxn_D_tot = 1.0, 10.0
            edge_energy = -16.0
            edge_type = 'unknown'
            for G in source_graphs:
                if G.has_edge(u, v):
                    data = G.edges[u, v]
                    rxn_sigma = data.get('sigma', 1.0)
                    rxn_D_tot = data.get('D_tot', 10.0)
                    edge_type = data.get('type', 'unknown')
                    val = data.get('energy', -16.0)
                    # Convert default fallback values directly
                    if val == -1.0: val = -16.0
                    edge_energy = float(val)
                    break
            types.append(edge_type)
            energies.append(edge_energy)
            kons.append(_calculate_macroscopic_kon(ka=default_ka, sigma=rxn_sigma, D_tot=rxn_D_tot))
        return types, energies, kons

    safe_m_fwd = max(1, m_fwd)
    safe_m_rev = max(1, m_rev)

    if delta_n < 0:
        # Association
        types, energies, kons = _get_edge_properties(new_edges, products)
        if not types:
            return 0.0, 0.0
            
        all_same_type = (len(set(types)) == 1)
        
        if all_same_type:
            raw_kon = kons[0]
            single_energy_kT = energies[0]
            
            k_fwd_nm3 = raw_kon * safe_m_fwd
            k_fwd = k_fwd_nm3 * CONV_NM3_US_TO_UM_S
            
            Kc_single = np.exp(-single_energy_kT) / C0_uM  # exp(-E) * C0^-1
            baseline_k_off = (raw_kon * CONV_NM3_US_TO_UM_S) / Kc_single
            k_rev = baseline_k_off * safe_m_rev
        else:
            fastest_kon = max(kons)
            n_bonds = len(kons)
            effective_kon_nm3 = fastest_kon * n_bonds
            k_fwd = effective_kon_nm3 * CONV_NM3_US_TO_UM_S
            
            total_energy_kT = sum(energies)
            Kc_total = np.exp(-total_energy_kT) / C0_uM
            
            k_rev = k_fwd / Kc_total

    elif delta_n > 0:
        # Dissociation
        types, energies, kons = _get_edge_properties(broken_edges, reactants)
        if not types:
            return 0.0, 0.0
            
        all_same_type = (len(set(types)) == 1)
        
        if all_same_type:
            raw_kon = kons[0]
            single_energy_kT = energies[0]
            
            k_rev_nm3 = raw_kon * safe_m_rev
            k_rev = k_rev_nm3 * CONV_NM3_US_TO_UM_S
            
            Kc_single = np.exp(-single_energy_kT) / C0_uM
            baseline_k_off = (raw_kon * CONV_NM3_US_TO_UM_S) / Kc_single
            k_fwd = baseline_k_off * safe_m_fwd
        else:
            fastest_kon = max(kons)
            n_bonds = len(kons)
            effective_kon_nm3 = fastest_kon * n_bonds
            k_rev = effective_kon_nm3 * CONV_NM3_US_TO_UM_S
            
            total_energy_kT = sum(energies)
            Kc_total = np.exp(-total_energy_kT) / C0_uM
            
            k_fwd = k_rev / Kc_total

    else:
        # Isomerization
        def _sum_energy_kT(graphs):
            total_E = 0.0
            seen = set()
            for G in graphs:
                for u, v, d in G.edges(data=True):
                    k = frozenset([u, v])
                    if k not in seen:
                        e = float(d.get('energy', -16.0))
                        if e == -1.0: e = -16.0
                        total_E += e
                        seen.add(k)
            return total_E
            
        delta_G = _sum_energy_kT(products) - _sum_energy_kT(reactants)
        Kc_total = np.exp(-delta_G)
        
        baseline_kon_nm3 = _calculate_macroscopic_kon(ka=default_ka, sigma=1.0, D_tot=10.0)
        baseline_rate = baseline_kon_nm3 * CONV_NM3_US_TO_UM_S
        
        k_fwd = baseline_rate * Kc_total * safe_m_fwd
        k_rev = baseline_rate * safe_m_rev

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

def _node_match(n1, n2):
    return n1.get('type') == n2.get('type')

def _edge_match(e1, e2):
    return e1.get('type') == e2.get('type')

def _compute_multiplicity_dimer(G1: nx.Graph, G2: nx.Graph, G_prod: nx.Graph) -> Tuple[int, int]:
    from itertools import combinations
    m_rev = 0
    nodes = list(G_prod.nodes())
    
    seen_partitions = set()
    for size in [len(G1), len(G2)]:
        for A_nodes in combinations(nodes, size):
            A_set = frozenset(A_nodes)
            B_set = frozenset(set(nodes) - A_set)
            
            partition_key = frozenset([A_set, B_set])
            if partition_key in seen_partitions:
                continue
            seen_partitions.add(partition_key)
            
            c1 = G_prod.subgraph(A_set)
            c2 = G_prod.subgraph(B_set)
            
            if not nx.is_connected(c1) or not nx.is_connected(c2):
                continue
                
            if (nx.is_isomorphic(c1, G1, node_match=_node_match, edge_match=_edge_match) and 
                nx.is_isomorphic(c2, G2, node_match=_node_match, edge_match=_edge_match)):
                m_rev += 1
            elif len(G1) != len(G2) and (nx.is_isomorphic(c1, G2, node_match=_node_match, edge_match=_edge_match) and 
                  nx.is_isomorphic(c2, G1, node_match=_node_match, edge_match=_edge_match)):
                m_rev += 1

    m_fwd = max(1, m_rev)
    
    return m_fwd, m_rev

def _compute_multiplicity_transformation(G1: nx.Graph, G2: nx.Graph) -> Tuple[int, int]:
    m_fwd = 0
    types_in_g2 = set(d.get('type') for u, v, d in G2.edges(data=True))
    # Number of ways to add an edge in G1 to get G2
    for u in G1.nodes():
        for v in G1.nodes():
            if u != v and not G1.has_edge(u, v):
                for t in types_in_g2:
                    G1_temp = G1.copy()
                    G1_temp.add_edge(u, v, type=t)
                    if nx.is_isomorphic(G1_temp, G2, node_match=_node_match, edge_match=_edge_match):
                        m_fwd += 1
    m_fwd = m_fwd // 2  # Each edge is checked twice due to undirected node pairs

    m_rev = 0
    # Number of ways to remove an edge in G2 to get G1
    for u, v, d in list(G2.edges(data=True)):
        G2_temp = G2.copy()
        G2_temp.remove_edge(u, v)
        if nx.is_isomorphic(G2_temp, G1, node_match=_node_match, edge_match=_edge_match):
            m_rev += 1

    return m_fwd, m_rev

def generate_ode_model_from_system(
    system: System,
    max_complex_size: int = None,
    pdb_model=None,
    coarse_grainer=None,
    include_transformation_reactions: bool = False,
    default_ka: float = 120.0,
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
        default_ka: Microscopic activation-limited on-rate in nm^3/us. Default: 120.0.
               Keq is derived from bond energetics; k_off = k_on_macro / Kc.
               sigma and D_tot are now dynamically fetched from the molecule properties.

    Returns:
        Tuple of (complex_names, reaction_system) where:
            - complex_names: List of string identifiers for the complexes
            - reaction_system: Object containing a .reactions list of SimpleReaction objects
    """
    # Set default max size
    if max_complex_size is None:
        max_complex_size = 12

    # Default microscopic activation rate ka (nm^3/us)
    DEFAULT_KA = default_ka
    
    # Step 1: Build the full assembly graph directly from System
    G_full = _system_to_graph(system)
    
    # Check if full assembly exceeds max_complex_size
    if len(G_full.nodes) > max_complex_size:
        raise ValueError(
            f"Assembly has {len(G_full.nodes)} molecules, exceeding max_complex_size ({max_complex_size}). "
            f"Skipping ODE generation. Increase max_complex_size parameter if you want to force calculation (RISKY)."
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
                m_fwd, m_rev = _compute_multiplicity_dimer(sub1, sub2, sub_prod)
                
                # Calculate Rates using SCK macroscopic formula
                k_on, k_off = _calculate_reaction_rates(
                    reactants=[sub1, sub2],
                    products=[sub_prod],
                    default_ka=DEFAULT_KA,
                    m_fwd=m_fwd,
                    m_rev=m_rev,
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
                m_fwd, m_rev = _compute_multiplicity_transformation(G1, G2)
                
                # Calculate Rates using SCK macroscopic formula
                k_fwd, k_rev = _calculate_reaction_rates(
                    reactants=[G1],
                    products=[G2],
                    default_ka=DEFAULT_KA,
                    m_fwd=m_fwd,
                    m_rev=m_rev,
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

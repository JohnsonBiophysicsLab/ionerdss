import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from itertools import combinations
from multiprocessing import Pool, cpu_count

def is_connected(G, nodes):
    """Check if nodes induce a connected subgraph in G."""
    return nx.is_connected(G.subgraph(nodes))

def typed_wl_hash(Gsub):
    """Get WL graph hash using node and edge types."""
    return weisfeiler_lehman_graph_hash(Gsub, node_attr="type", edge_attr="type")

def all_unique_induced_splits(G):
    """Yield all unique ways to split G into two connected, typed-isomorphic-aware induced subgraphs."""
    nodes = list(G.nodes)
    n = len(nodes)
    seen = set()

    for r in range(1, n // 2 + 1):
        for A in combinations(nodes, r):
            A = set(A)
            B = set(nodes) - A

            if not is_connected(G, A) or not is_connected(G, B):
                continue

            h1 = typed_wl_hash(G.subgraph(A))
            h2 = typed_wl_hash(G.subgraph(B))
            key = tuple(sorted((h1, h2)))

            if key in seen:
                continue
            seen.add(key)
            yield A, B

def deduplicate_species(species):
    """Filter out isomorphic species using WL hash."""
    seen_hashes = set()
    unique_species = []
    for sp in species:
        h = typed_wl_hash(sp)
        if h not in seen_hashes:
            seen_hashes.add(h)
            unique_species.append(sp)
    return unique_species

def get_broken_edges(G, part1, part2):
    """Return edges between part1 and part2, with type attributes."""
    cut_edges = []
    for u in part1:
        for v in G.neighbors(u):
            if v in part2:
                edge_type = G[u][v].get("type", None)
                cut_edges.append((u, v, edge_type))
    return cut_edges

def is_single_bond_change(G1, G2):
    """Check if G1 and G2 differ by exactly one edge (added or removed), not both."""
    edges1 = set((min(u, v), max(u, v), d['type']) for u, v, d in G1.edges(data=True))
    edges2 = set((min(u, v), max(u, v), d['type']) for u, v, d in G2.edges(data=True))

    diff1 = edges1 - edges2
    diff2 = edges2 - edges1

    # Only bond(s) formed or broken, not both
    return (len(diff1) >= 1 and len(diff2) == 0) or (len(diff1) == 0 and len(diff2) >= 1)

def find_single_bond_transformations(species):
    """Find pairs of species that differ by exactly one bond (added or removed)."""
    transformations = []
    seen = set()
    n = len(species)
    for i in range(n):
        for j in range(i+1, n):
            G1 = species[i]
            G2 = species[j]
            if set(G1.nodes) != set(G2.nodes):  # must be same node set
                continue
            if is_single_bond_change(G1, G2):
                h1, h2 = typed_wl_hash(G1), typed_wl_hash(G2)
                key = tuple(sorted((h1, h2)))
                if key not in seen:
                    seen.add(key)
                    transformations.append((G1, G2))
    return transformations

def compute_reactions_for_species(specie):
    """Compute reactions (split pairs) for a single species."""
    reactions = []
    for part1, part2 in all_unique_induced_splits(specie):
        reactions.append((part1, part2, specie))
    return reactions

def find_all_dimer_reactions(species, use_multiprocessing=False):
    """Compute all reactions across a list of species with optional multiprocessing."""
    species = deduplicate_species(species)

    if use_multiprocessing:
        with Pool(cpu_count()) as pool:
            results = pool.map(compute_reactions_for_species, species)
        reactions = [r for group in results for r in group]
    else:
        reactions = []
        for specie in species:
            reactions.extend(compute_reactions_for_species(specie))
    return reactions

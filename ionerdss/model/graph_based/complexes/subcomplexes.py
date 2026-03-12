import networkx as nx
from itertools import combinations
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from networkx.algorithms.components import connected_components

def powerset_connected_nodes(nodes):
    """Generate all non-empty subsets of nodes, up to len(nodes)"""
    nodes = list(nodes)
    for r in range(1, len(nodes) + 1):
        yield from combinations(nodes, r)

def _enumerate_cis_reverse_search(comp_G):
    """Enumerate all connected induced subgraphs using reverse search methodology."""
    nodes = list(comp_G.nodes)
    try:
        nodes.sort()
    except TypeError:
        nodes.sort(key=str)
    
    node_to_idx = {n: i for i, n in enumerate(nodes)}
    
    def _is_connected_induced(S_nodes):
        """Connectivity of induced subgraph comp_G[S_nodes] without constructing full graph objects."""
        # Fast path: single node
        it = iter(S_nodes)
        start = next(it, None)
        if start is None:
            return False
        if len(S_nodes) == 1:
            return True

        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for nb in comp_G.adj[v]:
                if nb in S_nodes and nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        return len(seen) == len(S_nodes)
    
    def parent_of(S_set):
        """Return parent set P(S) by removing the largest removable vertex (keeping connectivity)."""
        if len(S_set) <= 1:
            return None

        # Try removing vertices from largest to smallest.
        # In practice, the max one is often removable quickly.
        ordered = sorted(S_set, key=lambda n: node_to_idx[n], reverse=True)

        for v in ordered:
            T = S_set - {v}
            # connectivity check on induced subgraph T, cheap BFS
            if _is_connected_induced(T):
                return T
        return None  # should not happen for connected S_set

    def explore(S_set):
        N_S = set()
        for v in S_set:
            for nb in comp_G.adj[v]:
                if nb not in S_set:
                    N_S.add(nb)
                    
        for y in N_S:
            S_prime = S_set | {y}
            if parent_of(S_prime) == S_set:
                yield S_prime
                yield from explore(S_prime)
                
    for v in nodes:
        yield {v}
        yield from explore({v})

def get_unique_fully_connected_subgraphs(G, use_reverse_search=False):
    seen_hashes = set()
    unique_subgraphs = []

    for component_nodes in connected_components(G):
        comp_nodes = list(component_nodes)
        node_to_int = {n: i for i, n in enumerate(comp_nodes)}
        int_to_node = {i: n for i, n in enumerate(comp_nodes)}
        
        # Build an integer-labeled component graph once
        comp_Gi = nx.relabel_nodes(G.subgraph(comp_nodes), node_to_int, copy=True)
        
        # cache degrees on the integer component graph
        deg_i = dict(comp_Gi.degree())

        if use_reverse_search:
            # now yields sets of ints
            subset_generator = _enumerate_cis_reverse_search(comp_Gi)
        else:
            subset_generator = powerset_connected_nodes(comp_Gi.nodes)

        for node_subset in subset_generator:
            H = comp_Gi.subgraph(node_subset)

            if len(H) == 1:
                node = next(iter(node_subset))
                if deg_i[node] == 0:
                    continue
            else:
                if (not use_reverse_search) and (not nx.is_connected(H)):
                    continue

            # Compute hash directly on the integer subgraph
            wl_hash = weisfeiler_lehman_graph_hash(H, node_attr="type", edge_attr="type")

            if wl_hash not in seen_hashes:
                seen_hashes.add(wl_hash)
                # Remap the integer nodes back to original nodes for the resulting subgraph
                original_nodes = [int_to_node[n] for n in node_subset]
                H_orig = G.subgraph(original_nodes)
                unique_subgraphs.append(H_orig)

    return unique_subgraphs

def all_nonempty_proper_subsets(s):
    """All non-empty subsets of s that are not equal to s itself."""
    s = list(s)
    return (set(combo) for r in range(1, len(s)) for combo in combinations(s, r))


import networkx as nx
from itertools import combinations, chain
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from networkx.algorithms.components import connected_components

def powerset_connected_nodes(nodes):
    """Generate all non-empty subsets of nodes, up to len(nodes)"""
    nodes = list(nodes)
    for r in range(1, len(nodes) + 1):
        yield from combinations(nodes, r)

def get_unique_fully_connected_subgraphs(G):
    seen_hashes = set()
    unique_subgraphs = []

    full_degrees = dict(G.degree())  # cache full graph degrees

    for component_nodes in connected_components(G):
        for node_subset in powerset_connected_nodes(component_nodes):
            #H = G.subgraph(node_subset).copy() <- Avoid copy graph, 20% speed up
            H = G.subgraph(node_subset)

            if len(H) == 1:
                # Allow size-1 subgraphs only if connected in G
                node = node_subset[0]
                if full_degrees[node] == 0:
                    continue
            else:
                # Require every node in subgraph to have degree >= 1 in H
                if not nx.is_connected(H):
                    continue

            # Use canonical hash for deduplication
            H_relabel = nx.convert_node_labels_to_integers(H)
            wl_hash = weisfeiler_lehman_graph_hash(H_relabel, node_attr="type", edge_attr="type")

            if wl_hash not in seen_hashes:
                seen_hashes.add(wl_hash)
                unique_subgraphs.append(H)

    return unique_subgraphs

def all_nonempty_proper_subsets(s):
    """All non-empty subsets of s that are not equal to s itself."""
    s = list(s)
    return (set(combo) for r in range(1, len(s)) for combo in combinations(s, r))

if __name__ == "__main__":
    # Example graph (your complex test case)
    G = nx.Graph()
    G.add_node(0, type="A")
    G.add_node(1, type="A")
    G.add_node(2, type="A")
    G.add_node(3, type="A")

    G.add_edge(0, 1, type="ab")
    G.add_edge(1, 2, type="bc")
    G.add_edge(2, 3, type="ca")
    G.add_edge(3, 1, type="ab")

    get_unique_fully_connected_subgraphs(G)

# Final fix: relabel both graphs canonically FIRST, then compare edges and test isomorphism

import networkx as nx
from collections import defaultdict
from itertools import combinations
from networkx.algorithms.isomorphism import GraphMatcher
from ..complexes.subcomplexes import get_unique_fully_connected_subgraphs

# Match functions
# We define that two graphs are isomorphic if there exists some one-on-one mapping
# for nodes, edges, node types, and edge types
def node_match(n1, n2):
    return n1["type"] == n2["type"]

def edge_match(e1, e2):
    return e1["type"] == e2["type"]

def are_type_isomorphic(G1, G2):
    gm = GraphMatcher(G1, G2, node_match=node_match, edge_match=edge_match)
    return gm.is_isomorphic()

# Canonical relabeling that preserves type attributes
def relabel_graph_by_type(G):
    type_counts = defaultdict(int)
    mapping = {}
    new_attrs = {}

    for n, d in G.nodes(data=True):
        t = d["type"]
        label = f"{t}{type_counts[t]}"
        mapping[n] = label
        new_attrs[label] = {"type": t}
        type_counts[t] += 1

    G_relabel = nx.relabel_nodes(G, mapping, copy=True)
    nx.set_node_attributes(G_relabel, new_attrs)
    return G_relabel

# Main transformation logic
def is_transformable_by_forming_or_breaking_canonically(G1, G2):
    # Relabel both graphs canonically by type
    G1c = relabel_graph_by_type(G1)
    G2c = relabel_graph_by_type(G2)

    if sorted([d["type"] for _, d in G1c.nodes(data=True)]) != sorted([d["type"] for _, d in G2c.nodes(data=True)]):
        return False, None, None

    def edge_set(G):
        return set((min(u, v), max(u, v), d["type"]) for u, v, d in G.edges(data=True))

    e1 = edge_set(G1c)
    e2 = edge_set(G2c)

    # Try G1 -> G2 by forming edges
    if e1 < e2:
        edges_added = list(e2 - e1)
        G1_aug = G1c.copy()
        for u, v, t in edges_added:
            G1_aug.add_edge(u, v, type=t)
        if are_type_isomorphic(G1_aug, G2c):
            return True, "added", edges_added

    # Try G1 -> G2 by removing edges
    if e1 > e2:
        edges_removed = list(e1 - e2)
        G1_red = G1c.copy()
        for u, v, t in edges_removed:
            if G1_red.has_edge(u, v) and G1_red[u][v].get("type") == t:
                G1_red.remove_edge(u, v)
        if are_type_isomorphic(G1_red, G2c):
            return True, "removed", edges_removed

    return False, None, None

def find_all_transformable_subgraph_pairs(G, subgraphs = None):
    """
    This function calls `molecule_gen.get_unique_fully_connected_subgraphs` to
    get all the subgraphs first and the use the function
    `is_transformable_by_forming_or_breaking_canonically` above to check 
    for the pairs.
    """
    # get all subgraphs if subgraphs is not given
    if subgraphs is None:
        subgraphs = get_unique_fully_connected_subgraphs(G)

    # iterate over all non-repeated unordered pairs of subgraphs and
    # get transformable pairs
    transformable_pairs = []
    for G1, G2 in combinations(subgraphs, 2):
        is_transformable, direction, list_of_edges_changed =\
            is_transformable_by_forming_or_breaking_canonically(G1, G2)
        if is_transformable:
            transformable_pairs.append((G1, G2, direction, list_of_edges_changed))
    
    return transformable_pairs

if __name__ == "__main__":
    # Example graph
    G = nx.Graph()
    G.add_node(0, type="A")
    G.add_node(1, type="B")
    G.add_node(2, type="C")
    G.add_node(3, type="A")
    G.add_edge(0, 1, type="ab")
    G.add_edge(1, 2, type="bc")
    G.add_edge(2, 3, type="ca")
    G.add_edge(3, 1, type="ab")

    transformable_pairs = find_all_transformable_subgraph_pairs(G)

    for t1, t2, dir, list_of_edges_changed in transformable_pairs:
        print(t1.nodes)
        print(t2.nodes)
        print(dir)
        print(list_of_edges_changed)
        
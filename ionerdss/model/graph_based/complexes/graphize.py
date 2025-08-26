# graphize.py
# simple_graph.py
from collections import defaultdict
import networkx as nx
import numpy as np

def _coord_tuple(p):
    """Return (x,y,z) from a Coords-like object or array-like."""
    if p is None:
        return None
    if hasattr(p, "to_numpy"):
        arr = p.to_numpy()
        return (float(arr[0]), float(arr[1]), float(arr[2]))
    # Coords(x=..., y=..., z=...) or plain tuple/list/np array
    if hasattr(p, "x") and hasattr(p, "y") and hasattr(p, "z"):
        return (float(p.x), float(p.y), float(p.z))
    p = np.asarray(p, dtype=float)
    return (float(p[0]), float(p[1]), float(p[2]))

def build_simple_graph(cg_model, chains_map=None):
    """
    Build a minimal NetworkX graph from a coarse-grained model.

    Nodes:
      - 'type': molecule name (chain id, or chains_map[chain id] if provided)

    Edges:
      - 'type': interface template name (e.g., "A_1", "E_2") if available,
                otherwise falls back to an arbitrary placeholder.

    Parameters
    ----------
    cg_model : dict
        Expects:
          - 'chains' : list of chain objects (with .id) or plain ids
          - 'interfaces' : list[list[str]] partner chain IDs per chain
          - optionally 'interface_templates' or molecules with interface_list
    chains_map : dict, optional
        Map from chain_id -> molecule name (e.g., template/type). If provided,
        node 'name' will use this mapping; otherwise it uses the chain id.

    Returns
    -------
    nx.Graph
    """
    chains = [getattr(c, "id", c) for c in cg_model["chains"]]
    interfaces = cg_model["interfaces"]
    idx_of = {cid: i for i, cid in enumerate(chains)}

    def mol_name(cid):
        return chains_map.get(cid, cid) if chains_map else cid

    G = nx.Graph()
    for i, cid in enumerate(chains):
        G.add_node(i, type=mol_name(cid))

    # If model has molecule objects with interface_list → use template names
    template_lookup = {}
    if "molecules" in cg_model:
        for mol in cg_model["molecules"]:
            for iface in getattr(mol, "interface_list", []):
                if iface.my_template and hasattr(iface.my_template, "name"):
                    template_lookup[(mol.name, iface.name)] = iface.my_template.name

    for i, partners in enumerate(interfaces):
        for iface_idx, partner_cid in enumerate(partners):
            j = idx_of[partner_cid]
            if i >= j:  # undirected dedupe
                continue

            # Default label
            edge_label = f"{mol_name(chains[i])}_{iface_idx+1}"

            # If we know the interface template name, use it
            key = (chains[i], partner_cid)
            if key in template_lookup:
                edge_label = template_lookup[key]

            G.add_edge(i, j, type=edge_label)

    return G


def annotate_edges_by_cycles(G):
    """
    Label edges by cycle membership:
      - 'hex' if the edge belongs to any 6-cycle
      - 'tri' if the edge belongs to any 3-cycle (and not already 'hex')
      - 'di' otherwise

    This mimics the 5l93 labeling style (hex ring edges vs triangle edges vs
    leftover dimer-like edges).
    """
    # reset existing tag if present
    for u, v in G.edges:
        G[u][v].pop("type", None)

    # find simple cycles
    cycles = nx.cycle_basis(G)
    # Mark hex first (takes precedence)
    for cyc in cycles:
        L = len(cyc)
        edges = [(cyc[i], cyc[(i + 1) % L]) for i in range(L)]
        if L == 6:
            for u, v in edges:
                if G.has_edge(u, v):
                    G[u][v]["type"] = "hex"
    # Then triangles
    for cyc in cycles:
        L = len(cyc)
        if L == 3:
            edges = [(cyc[i], cyc[(i + 1) % L]) for i in range(L)]
            for u, v in edges:
                if G.has_edge(u, v) and "type" not in G[u][v]:
                    G[u][v]["type"] = "tri"
    # Remaining = 'di'
    for u, v in G.edges:
        if "type" not in G[u][v]:
            G[u][v]["type"] = "di"
    return G

def networkx_graph_to_string(G: nx.Graph) -> str:
    """
    Convert a NetworkX graph into a human-readable string with nodes and edges.

    Parameters
    ----------
    G : nx.Graph
        The NetworkX graph to summarize.

    Returns
    -------
    str
        A string containing the list of nodes with attributes and
        the list of edges with attributes.
    """
    lines = []
    lines.append("Nodes:")
    for n, attrs in G.nodes(data=True):
        lines.append(f"  {n}: {attrs}")

    lines.append("\nEdges:")
    for u, v, attrs in G.edges(data=True):
        lines.append(f"  {u} -- {v}, attrs={attrs}")

    return "\n".join(lines)

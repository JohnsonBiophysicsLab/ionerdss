'''
This script benchmarks the performance of the subcomplex enumeration methods.
It measures the time taken to enumerate the subcomplexes.
Also, when `verify_equivalence = True`, it compares the results of the powerset enumeration and Avis-Fukuda reverse search methods.

Usage:
python benchmark_subcomplexes.py
'''

import timeit
import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from ionerdss.model.graph_based.complexes.examples import graph_registry
from ionerdss.model.graph_based.complexes.subcomplexes import get_unique_fully_connected_subgraphs

def verify_equivalence(res1, res2):
    if len(res1) != len(res2):
        return False, f"Length mismatch: {len(res1)} vs {len(res2)}"
    
    # Hash all subgraphs canonically
    def get_hashes(subgraphs):
        hashes = set()
        for H in subgraphs:
            H_relabel = nx.convert_node_labels_to_integers(H)
            wl_hash = weisfeiler_lehman_graph_hash(H_relabel, node_attr="type", edge_attr="type")
            hashes.add(wl_hash)
        return hashes
        
    hash1 = get_hashes(res1)
    hash2 = get_hashes(res2)
    
    if hash1 != hash2:
        return False, "Hashes do not match completely."
        
    return True, "Identical"

def main():
    verify_equivalence = False
    print("=" * 60)
    print(f"{'Graph':<15} | {'Unique Subgraphs'} | {'Equivalence':<10} | {'Reverse=False (s)':<17} | {'Reverse=True (s)':<17}")
    print("-" * 60)
    
    for graph_name, graph_func in graph_registry.items():
        G = graph_func()
        
        # Verify correctness
        if graph_name in ('ring_25', 'complex_30mer'):
            # Skip powerset for massive graphs as iterations take > 20 mins
            res_reverse = get_unique_fully_connected_subgraphs(G, use_reverse_search=True)
            num_subgraphs = len(res_reverse)
            is_equiv = None
            equiv_msg = "Skipped"
        else:
            res_powerset = get_unique_fully_connected_subgraphs(G, use_reverse_search=False)
            res_reverse = get_unique_fully_connected_subgraphs(G, use_reverse_search=True)
            
            if verify_equivalence:
                is_equiv, equiv_msg = verify_equivalence(res_powerset, res_reverse)
                num_subgraphs = len(res_powerset)
            else:
                is_equiv = None
                equiv_msg = "N/A"
                num_subgraphs = len(res_powerset)
        
        # Benchmark time
        def run_powerset():
            get_unique_fully_connected_subgraphs(G, use_reverse_search=False)
            
        def run_reverse():
            get_unique_fully_connected_subgraphs(G, use_reverse_search=True)
            
        if graph_name == '5l93':
            num_runs = 1
            time_powerset = timeit.timeit(run_powerset, number=num_runs) / num_runs
            time_reverse = timeit.timeit(run_reverse, number=num_runs) / num_runs
        elif graph_name == 'ring_25':
            num_runs = 10
            # User set `time_powerset` to basically simulate N/A skipping the calculation
            time_powerset = timeit.timeit(run_reverse, number=num_runs) / num_runs
            time_reverse = timeit.timeit(run_reverse, number=num_runs) / num_runs
        elif graph_name == 'complex_30mer':
            num_runs = 1
            time_powerset = 0
            time_reverse = timeit.timeit(run_reverse, number=num_runs) / num_runs
        else:
            num_runs = 50
            time_powerset = timeit.timeit(run_powerset, number=num_runs) / num_runs
            time_reverse = timeit.timeit(run_reverse, number=num_runs) / num_runs

        powerset_str = f"{time_powerset:.5f}"
        
        print(f"{graph_name:<15} | {num_subgraphs:<16} | {equiv_msg:<10} | {powerset_str:<17} | {time_reverse:<17.5f}")

if __name__ == '__main__':
    main()

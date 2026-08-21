"""
test_subcomplexes_extraction.py

Unit tests for `get_unique_fully_connected_subgraphs` from
`ionerdss.model.graph_based.complexes.subcomplexes`.

This test suite validates the correctness and performance of subgraph enumeration
in various synthetic graph configurations:

1. `TestFullyConnectedSubgraphDetection8y7s`
   - A small symmetric 4-node complete graph with edge types
   - Expected to return 6 unique fully connected subgraphs

2. `TestFullyConnectedSubgraphDetectionHetero8mer`
   - A ring of 8 uniquely typed nodes and edges
   - Expected to return 57 unique subgraphs (sum of 1 to 8)

3. `TestFullyConnectedSubgraphDetectionAsymmetric`
   - A 4-node asymmetric graph with repeated node and edge types
   - Expected to return 10 unique subgraphs

Each test measures performance and emits a warning if runtime exceeds 0.1s,
and fails if it exceeds 10s.

Author: yying7@jh.edu
-------
Auto-generated and maintained as part of the ioNERDSS framework
for molecular complex modeling and reaction network generation.
"""

import unittest
import networkx as nx
import time
import warnings

# Replace with actual import path
from ionerdss.model.graph_based.complexes.subcomplexes import get_unique_fully_connected_subgraphs

class TestFullyConnectedSubgraphDetection8y7s(unittest.TestCase):
    def setUp(self):
        self.G = nx.Graph()
        self.G.add_node(0, type="X")
        self.G.add_node(1, type="X")
        self.G.add_node(2, type="X")
        self.G.add_node(3, type="X")

        self.G.add_edge(0, 1, type="a")
        self.G.add_edge(0, 2, type="b")
        self.G.add_edge(0, 3, type="c")
        self.G.add_edge(2, 3, type="a")
        self.G.add_edge(1, 3, type="b")
        self.G.add_edge(1, 2, type="c")

    def test_unique_subgraphs(self):
        start = time.time()
        result = get_unique_fully_connected_subgraphs(self.G)
        elapsed = time.time() - start

        if elapsed > 10.0:
            self.fail(f"ERROR: Test terminated — elapsed time {elapsed:.2f}s exceeds 10 seconds.")
        elif elapsed > 0.01:
            warnings.warn(f"WARNING: Elapsed time {elapsed:.2f}s exceeds 0.1 second.")

        self.assertEqual(
            len(result), 6,
            f"ERROR: Expected 6 unique subgraphs, but got {len(result)}."
        )

        print(f"Test passed in {elapsed:.4f} seconds.")

class TestFullyConnectedSubgraphDetectionHetero8mer(unittest.TestCase):
    def setUp(self):
        self.G = nx.Graph()
        self.G.add_node(0, type="A")
        self.G.add_node(1, type="B")
        self.G.add_node(2, type="C")
        self.G.add_node(3, type="D")
        self.G.add_node(4, type="E")
        self.G.add_node(5, type="F")
        self.G.add_node(6, type="G")
        self.G.add_node(7, type="H")

        self.G.add_edge(0, 1, type="ab")
        self.G.add_edge(1, 2, type="bc")
        self.G.add_edge(2, 3, type="cd")
        self.G.add_edge(3, 4, type="de")
        self.G.add_edge(4, 5, type="ef")
        self.G.add_edge(5, 6, type="fg")
        self.G.add_edge(6, 7, type="gh")
        self.G.add_edge(7, 0, type="ha")

    def test_unique_subgraphs(self):
        start = time.time()
        result = get_unique_fully_connected_subgraphs(self.G)
        elapsed = time.time() - start

        if elapsed > 10.0:
            self.fail(f"ERROR: Test terminated — elapsed time {elapsed:.2f}s exceeds 10 seconds.")
        elif elapsed > 0.1:
            warnings.warn(f"WARNING: Elapsed time {elapsed:.2f}s exceeds 0.1 second.")

        self.assertEqual(
            len(result), 57,
            f"ERROR: Expected 57 unique subgraphs, but got {len(result)}."
        ) # 7 * 8 + 1 = 57

        print(f"Test passed in {elapsed:.4f} seconds.")
        
class TestFullyConnectedSubgraphDetectionAsymmetric(unittest.TestCase):
    def setUp(self):
        self.G = nx.Graph()
        # Add typed nodes
        self.G.add_node(0, type="A")
        self.G.add_node(1, type="B")
        self.G.add_node(2, type="C")
        self.G.add_node(3, type="A")

        # Add typed edges
        self.G.add_edge(0, 1, type="ab")
        self.G.add_edge(1, 2, type="bc")
        self.G.add_edge(2, 3, type="ca")
        self.G.add_edge(1, 3, type="ab")

    def test_unique_subgraphs(self):
        start = time.time()
        result = get_unique_fully_connected_subgraphs(self.G)
        elapsed = time.time() - start

        if elapsed > 10.0:
            self.fail(f"ERROR: Test terminated — elapsed time {elapsed:.2f}s exceeds 10 seconds.")
        elif elapsed > 0.1:
            warnings.warn(f"WARNING: Elapsed time {elapsed:.2f}s exceeds 0.1 second.")

        self.assertEqual(
            len(result), 10,
            f"ERROR: Expected 10 unique subgraphs, but got {len(result)}."
        ) # 3 + 3 + 3 + 1 = 10

        print(f"Test passed in {elapsed:.4f} seconds.")

if __name__ == '__main__':
    unittest.main()

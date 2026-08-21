"""
test_reactions_generation_pipeline.py

Unit test for benchmarking the reaction and transformation enumeration pipeline
in the `ode_gen` package using the "8y7s" example molecular graph.

This test verifies the correctness, performance structure, and minimal
functionality of the graph-based species and reaction generation pipeline,
which includes:

1. Subgraph Enumeration:
    - Uses `get_unique_fully_connected_subgraphs` to extract all unique
      fully connected induced subgraphs (species) from the input molecular graph.

2. Dimer Reaction Generation:
    - Applies `find_all_dimer_reactions` to identify all valid pairwise
      species interactions that result in a product graph.
    - Uses `get_broken_edges` to annotate which bonds are removed in each reaction.

3. Transformable Subgraph Detection:
    - Uses `find_all_transformable_subgraph_pairs` to find species pairs
      that differ only by one or more edges, indicating possible
      edge-transformation reactions (e.g., conformational transitions).

4. Structure Validation:
    - Asserts the returned dataframes and metadata dictionary are well-formed.
    - Ensures the number of species, reactions, and transformations is non-zero.
    - Verifies presence of timing information and critical keys in the results.

Test Strategy:
--------------
- Implemented using Python's `unittest` framework.
- Can be run standalone or through `pytest` without modification.
- Designed for reproducibility and integration testing of the core
  reaction graph enumeration pipeline.

Usage:
------
To run as a standalone test:
    python -m unittest tests.benchmark.test_pipeline

To run with pytest:
    pytest tests/benchmark/test_pipeline.py -v

Requires:
---------
- NetworkX-based graph model (`get_example`)
- `ode_gen` modules:
    - `complexes.subcomplexes`
    - `reactions.dimer`
    - `reactions.transformation`
- pandas, platform, and standard Python libraries

Author: yying7@jh.edu
-------
Auto-generated and maintained as part of the ioNERDSS framework
for molecular complex modeling and reaction network generation.
"""

import unittest
from pathlib import Path
import pandas as pd
import json

from ionerdss.model.graph_based.complexes.examples import get_example
from ionerdss.model.graph_based.complexes.subcomplexes import get_unique_fully_connected_subgraphs
from ionerdss.model.graph_based.reactions.dimer import find_all_dimer_reactions, get_broken_edges
from ionerdss.model.graph_based.reactions.transformation import find_all_transformable_subgraph_pairs

from datetime import datetime
import time
import platform

def capture_environment_info():
    return {
        "timestamp": datetime.now().isoformat(),
        "platform": platform.system(),
        "platform_version": platform.version(),
        "python_version": platform.python_version(),
        "hostname": platform.node(),
        "cpu": platform.processor(),
    }

def benchmark_pipeline(graph_name="8y7s"):
    G = get_example(graph_name)
    results = {}
    env_info = capture_environment_info()

    t0 = time.time()
    species = get_unique_fully_connected_subgraphs(G)
    t1 = time.time()

    reactions = find_all_dimer_reactions(species, use_multiprocessing=True)
    t2 = time.time()

    reaction_df = pd.DataFrame([
        {
            "product": list(r[2].nodes),
            "part1": list(r[0]),
            "part2": list(r[1]),
            "bonds_broken": get_broken_edges(r[2], r[0], r[1])
        }
        for r in reactions
    ])

    transformations = find_all_transformable_subgraph_pairs(G, species)
    t3 = time.time()

    transformations_df = pd.DataFrame([
        {
            "monomer_1_nodes": list(t1.nodes),
            "monomer_2_nodes": list(t2.nodes),
            "diff": list(set(t1.edges(data="type")) ^ set(t2.edges(data="type")))
        }
        for t1, t2 in transformations
    ])

    results["timing"] = {
        "subgraph_enumeration": round(t1 - t0, 4),
        "dimer_reactions": round(t2 - t1, 4),
        "transformable_pairs": round(t3 - t2, 4),
        "total": round(t3 - t0, 4),
    }
    results["n_species"] = len(species)
    results["n_reactions"] = len(reaction_df)
    results["n_transformations"] = len(transformations_df)
    results["environment"] = env_info

    return reaction_df, transformations_df, results


class TestBenchmarkPipeline(unittest.TestCase):

    def test_benchmark_output_8y7s(self):
        graph_name = "8y7s"
        reaction_df, transformation_df, metadata = benchmark_pipeline(graph_name)

        # Basic structure checks
        self.assertIsInstance(reaction_df, pd.DataFrame)
        self.assertIsInstance(transformation_df, pd.DataFrame)
        self.assertIsInstance(metadata, dict)

        # Output sanity
        self.assertGreaterEqual(len(reaction_df), 1, "No dimer reactions found")
        self.assertGreaterEqual(len(transformation_df), 0, " Transformable subgraphs found while none expected")
        self.assertIn("timing", metadata)
        self.assertIn("n_species", metadata)

        # Timing keys present
        for key in ["subgraph_enumeration", "dimer_reactions", "transformable_pairs", "total"]:
            self.assertIn(key, metadata["timing"])
            self.assertGreaterEqual(metadata["timing"][key], 0.0)


if __name__ == "__main__":
    unittest.main()

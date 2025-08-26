"""
tests/test_visualize.py

Unit tests for `visualize.py`.

Purpose
-------
Ensure the visualization utilities for coarse-grained models:
- plot_interfaces
- plot_coarse_grain_model
- save_coarse_grained_structure

behave correctly on simple synthetic inputs, without requiring a real PDB or
full pipeline.

Testing strategy
----------------
1) plot_interfaces
   - Creates a 3D matplotlib axis.
   - Plots two dummy interface coordinates.
   - Verifies a scatter object is added to the axis.

2) plot_coarse_grain_model
   - Constructs a minimal model dict with two chains, two COMs, and one interface each.
   - Calls the plotting function inside a non-interactive backend.
   - Confirms the figure/axes are generated without error.

3) save_coarse_grained_structure
   - Uses a temporary directory.
   - Writes a mock cg_model with one chain, COM, and interface.
   - Checks that both the `.cif` and `.pml` files are created and contain expected keywords.

Assumptions
-----------
- `Coords` class has `.x`, `.y`, `.z`, and `.to_numpy()` methods.
- matplotlib can run in headless mode with Agg backend.
- File I/O is permitted in temporary directories.
"""

import unittest
import tempfile
import os
import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

from ionerdss.model.pdb import visualize


class DummyCoords:
    """Minimal stand-in for ionerdss.Coords."""
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z

    def to_numpy(self):
        return np.array([self.x, self.y, self.z])
    
    def __iter__(self):
        # Makes zip() work
        return iter((self.x, self.y, self.z))


class TestVisualize(unittest.TestCase):
    def test_plot_interfaces_adds_points(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        coords = [DummyCoords(1, 2, 3), DummyCoords(-1, -2, -3)]

        visualize.plot_interfaces(ax, coords, color="blue", label="test", size=20)

        # Axes3D stores scatter collections in ax.collections
        self.assertTrue(any(hasattr(col, "get_offsets") or hasattr(col, "_offsets3d")
                            for col in ax.collections))

    def test_plot_coarse_grain_model_runs(self):
        model = {
            "COMs": [DummyCoords(0, 0, 0), DummyCoords(1, 1, 1)],
            "interface_coords": [[DummyCoords(0.5, 0, 0)], [DummyCoords(1.5, 1, 1)]],
            "interfaces": [[1], [0]],
            "chains": ["A", "B"],
        }

        # Should run without raising and return a figure
        visualize.plot_coarse_grain_model(model, show_reactions=True)

    def test_save_coarse_grained_structure_creates_files(self):
        cg_model = {
            "chains": [type("C", (), {"id": "A"})()],  # dummy chain with id attribute
            "COMs": [DummyCoords(0, 0, 0)],
            "interface_coords": [[DummyCoords(1, 0, 0)]],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_file = os.path.join(tmpdir, "dummy.pdb")
            with open(pdb_file, "w", encoding="utf-8") as fh:
                fh.write("HEADER DUMMY\n")

            visualize.save_coarse_grained_structure(cg_model, tmpdir, pdb_file)

            cif_path = os.path.join(tmpdir, "updated_coarse_grained_structure.cif")
            pml_path = os.path.join(tmpdir, "updated_visualize_coarse_grained.pml")
            self.assertTrue(os.path.exists(cif_path))
            self.assertTrue(os.path.exists(pml_path))

            with open(cif_path, encoding="utf-8") as f:
                content = f.read()
                self.assertIn("Coarse-grained structure CIF", content)

            with open(pml_path, encoding="utf-8") as f:
                content = f.read()
                self.assertIn("PyMOL script", content)
                self.assertIn("pseudoatom", content)


if __name__ == "__main__":
    unittest.main()

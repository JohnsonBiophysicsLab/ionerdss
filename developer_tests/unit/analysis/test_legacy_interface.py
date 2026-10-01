"""
Unit tests for ionerdss.analysis.legacy.interface

LegacyPlotInterface.plot_figure routes on simulations, x and y, and passes every
other kwarg on to the analyzer.plot method it picks, and from there to Matplotlib.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
import numpy as np
from ionerdss.analysis import Analyzer, LegacyPlotInterface

# Final running totals of two runs with different size distributions
FINAL_MATRICES = {
    "1": [[40, 2, 0], [4, 30, 0], [1, 1, 20]],
    "2": [[10, 1, 0], [2, 30, 1], [0, 1, 60]],
}


class TestLegacyPlotInterface(unittest.TestCase):
    """plot_figure calls in the old style, for every mapping it supports."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        for sim_id, matrix in FINAL_MATRICES.items():
            data_dir = Path(self.temp_dir) / sim_id / "DATA"
            data_dir.mkdir(parents=True)
            with open(data_dir / "transition_matrix_time.dat", 'w') as f:
                f.write("time: 0\n")
                f.write("transition matrix for each mol type: \n")
                f.write("A\n")
                for _ in matrix:
                    f.write(" 0 0 0\n")
                f.write("time: 0.01\n")
                f.write("transition matrix for each mol type: \n")
                f.write("A\n")
                for row in matrix:
                    f.write("".join(f" {x}" for x in row) + "\n")

        self.analyzer = Analyzer(self.temp_dir)
        # os.walk order is not sorted; index the runs by ID
        self.analyzer.simulations.sort(key=lambda sim: sim.id)
        self.legacy = LegacyPlotInterface(self.analyzer)

    def tearDown(self):
        plt.close("all")
        shutil.rmtree(self.temp_dir)

    def test_line_free_energy(self):
        """line / size / free_energy draws the free energy of the chosen simulation."""
        ax = self.legacy.plot_figure("line", x="size", y="free_energy", simulations=[1])

        self.assertIsInstance(ax, Axes)
        line = ax.get_lines()[0]
        self.assertEqual(line.get_label(), "2")
        expected = self.analyzer.compute_free_energy(self.analyzer.get_simulation(1))
        np.testing.assert_allclose(line.get_ydata(), expected["free_energy"])

    def test_line_probability(self):
        """line / size / *probability* draws growth and shrinkage probabilities."""
        ax = self.legacy.plot_figure(
            "line", x="size", y="symmetric_association_probability", simulations=[0]
        )

        self.assertIsInstance(ax, Axes)
        self.assertEqual(
            [line.get_label() for line in ax.get_lines()],
            ["Growth (Assoc)", "Shrinkage (Dissoc)"],
        )

    def test_hist_count(self):
        """hist / size / *count* draws the size distribution."""
        ax = self.legacy.plot_figure("hist", x="size", y="complex_count", simulations=[0])

        self.assertIsInstance(ax, Axes)
        heights = [bar.get_height() for bar in ax.patches]
        expected = self.analyzer.compute_size_distribution(self.analyzer.get_simulation(0))
        np.testing.assert_allclose(heights, expected["probability"])

    def test_heatmap_with_axis_kwargs(self):
        """heatmap ignores x and y rather than passing them to plot_heatmap."""
        ax = self.legacy.plot_figure("heatmap", x="size", y="size", simulations=[0])
        self.assertIsInstance(ax, Axes)

    def test_default_simulation_is_first(self):
        """Without simulations, the first simulation is plotted."""
        ax = self.legacy.plot_figure("line", x="size", y="free_energy")
        self.assertEqual(ax.get_lines()[0].get_label(), "1")

    def test_other_kwargs_reach_matplotlib(self):
        """Styling kwargs still arrive at the Matplotlib call."""
        ax = self.legacy.plot_figure("line", x="size", y="free_energy", color="red")
        self.assertEqual(to_rgba(ax.get_lines()[0].get_color()), to_rgba("red"))

    def test_unmapped_combination_returns_none(self):
        """An old combination with no mapping warns and returns None."""
        self.assertIsNone(self.legacy.plot_figure("line", x="time", y="species"))


if __name__ == '__main__':
    unittest.main(verbosity=2)

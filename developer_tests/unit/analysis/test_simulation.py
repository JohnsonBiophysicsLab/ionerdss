"""
Unit tests for ionerdss.analysis.core.simulation

NERDSS never resets its transition matrix or lifetime lists (write_transition.cpp
writes the running totals), so every time point in transition_matrix_time.dat
holds everything recorded since the start of the run.
"""

import unittest
import tempfile
import shutil
from pathlib import Path
import numpy as np
from ionerdss.analysis import Analyzer

# Running totals at each time point, as NERDSS writes them: every matrix entry
# only grows, and every lifetime list extends the one before it.
SNAPSHOTS = [
    (0.0, [[0, 0, 0], [0, 0, 0], [0, 0, 0]], {1: [], 2: [], 3: []}),
    (0.01, [[10, 1, 0], [2, 20, 0], [0, 1, 30]], {1: [0.5, 0.25], 2: [0.125], 3: []}),
    (0.02, [[15, 2, 0], [3, 50, 1], [1, 1, 80]], {1: [0.5, 0.25, 0.75], 2: [0.125], 3: []}),
]


def _write_nerdss_transition_file(path: Path) -> None:
    """Writes SNAPSHOTS in the layout of NERDSS's write_transition.cpp."""
    with open(path, 'w') as f:
        for time, matrix, lifetimes in SNAPSHOTS:
            f.write(f"time: {time}\n")
            f.write("transition matrix for each mol type: \n")
            f.write("A\n")
            for row in matrix:
                f.write("".join(f" {x}" for x in row) + "\n")
            f.write("lifetime for each mol type: \n")
            f.write("A\n")
            for size, values in lifetimes.items():
                f.write(f"size of the cluster:{size}\n")
                f.write("".join(f" {x}" for x in values) + "\n")


class TestCumulativeTransitionData(unittest.TestCase):
    """Reading NERDSS's running totals without adding them up again."""

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        data_dir = Path(self.temp_dir) / "1" / "DATA"
        data_dir.mkdir(parents=True)
        _write_nerdss_transition_file(data_dir / "transition_matrix_time.dat")

        self.analyzer = Analyzer(self.temp_dir)
        self.sim = self.analyzer.get_simulation(0)

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_whole_run_is_last_time_point(self):
        """The whole run's counts are the last time point, not the sum of all of them."""
        matrix = self.sim.get_transition_matrix()

        np.testing.assert_array_equal(matrix, SNAPSHOTS[-1][1])

        # A copy: changing it leaves the loaded data alone
        matrix[0, 0] = -1
        self.assertEqual(self.sim.data.transitions[-1]["matrix"][0, 0], 15)

    def test_time_range_is_difference_of_its_end_points(self):
        """A time range counts only the transitions between its first and last time points."""
        window = self.sim.get_transition_matrix(time_range=(0.01, 0.02))
        expected = np.array(SNAPSHOTS[2][1]) - np.array(SNAPSHOTS[1][1])
        np.testing.assert_array_equal(window, expected)

        # From time 0 the start is all zeros, so the whole run comes back
        whole = self.sim.get_transition_matrix(time_range=(0.0, 0.02))
        np.testing.assert_array_equal(whole, SNAPSHOTS[-1][1])

    def test_time_range_with_fewer_than_two_time_points_is_empty(self):
        """One time point in range gives nothing to take a difference against."""
        with self.assertLogs("ionerdss.analysis.core.simulation", level="WARNING"):
            one_point = self.sim.get_transition_matrix(time_range=(0.005, 0.015))
        self.assertEqual(one_point.size, 0)

        with self.assertLogs("ionerdss.analysis.core.simulation", level="WARNING"):
            no_points = self.sim.get_transition_matrix(time_range=(0.03, 0.04))
        self.assertEqual(no_points.size, 0)

    def test_lifetimes_are_each_counted_once(self):
        """Each recorded lifetime appears once, not once per later time point."""
        self.assertEqual(self.sim.get_lifetimes(1), [0.5, 0.25, 0.75])
        self.assertEqual(self.sim.get_lifetimes(2), [0.125])
        self.assertEqual(self.sim.get_lifetimes(3), [])
        self.assertEqual(self.sim.get_lifetimes(4), [])

    def test_size_distribution_uses_final_counts(self):
        """The size distribution, and so the free energy, comes from the final totals."""
        df = self.analyzer.compute_size_distribution(self.sim)

        final_counts = np.array(SNAPSHOTS[-1][1]).sum(axis=1)
        np.testing.assert_array_equal(df["count"], final_counts)
        np.testing.assert_allclose(df["probability"], final_counts / final_counts.sum())


if __name__ == '__main__':
    unittest.main(verbosity=2)

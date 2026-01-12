"""
test_ode_solver.py

Unit tests for reaction ODE solver utilities in the `ionerdss` package.

Tested Functions
----------------
- `calculate_macroscopic_reaction_rates`
- `reaction_dydt`
- `solve_reaction_ode`

Test Scenarios
--------------
1. `test_calculate_macroscopic_reaction_rates`
   - Tests basic macroscopic rate calculation with a simple linear system.

2. `test_reaction_dydt`
   - Tests time derivative calculation for a model enzymatic reaction:
        E + S --> ES       (k1)
        ES    --> E + P    (k2)
   - Checks that the rate of change `dy/dt` matches expected stoichiometric flow.

3. `test_solve_reaction_ode`
   - Ensures the ODE system solves without error on a simple bidirectional 2-species system.

Assumptions
-----------
- Reactions follow mass-action kinetics.
- Stoichiometry is represented with `reactant_matrix` and `product_matrix`.
- Reaction rate constants are provided as a flat list `k`.

This test suite does not generate plots, nor does it test solver accuracy over long time spans.

Author: yying7@jh.edu
-------
Auto-generated and maintained as part of the ioNERDSS framework
for molecular complex modeling and reaction network generation.
"""

import unittest
import numpy as np
from ionerdss import calculate_macroscopic_reaction_rates, reaction_dydt, solve_reaction_ode

class TestODESolver(unittest.TestCase):

    def test_calculate_macroscopic_reaction_rates(self):
        y = [1.0, 0.5]
        reactant_matrix = np.array([[1, 0], [0, 1]])
        k = [0.1, 0.2]

        result = calculate_macroscopic_reaction_rates(y, reactant_matrix, k)
        expected_result = np.array([0.1, 0.1])
        np.testing.assert_allclose(result, expected_result)

    def test_reaction_dydt(self):
        # Model system:
        # E + S -> ES, k1
        # ES -> E + P, k2
        #
        # y = [[E], [S], [ES], [P]]
        t = 0.1
        y = [0.1, 50.0, 0.1, 1.5]
        reactant_matrix = np.array([[1, 1, 0, 0], [0, 0, 1, 0]])
        product_matrix = np.array([[0, 0, 1, 0], [1, 0, 0, 1]])
        k = [100.0, 1.0]

        result = reaction_dydt(t, y, reactant_matrix, product_matrix, k)
        print(result)
        expected_result = np.array([-499.9, -500.0, 499.9, 0.1])
        np.testing.assert_allclose(result, expected_result)

    def test_solve_reaction_ode(self):
        t_span = (0, 10)
        y_initial = [1.0, 0.5]
        reactant_matrix = np.array([[1, 0], [0, 1]])
        product_matrix = np.array([[0, 1], [1, 0]])
        k = [0.1, 0.2]

        with self.subTest(msg="Check if solve_reaction_ode runs without errors"):
            solve_reaction_ode(reaction_dydt, t_span, y_initial, reactant_matrix, product_matrix, k, plotting=False)

if __name__ == '__main__':
    unittest.main()

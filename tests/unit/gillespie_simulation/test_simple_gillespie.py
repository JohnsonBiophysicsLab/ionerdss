"""
test_simple_gillespie.py

Unit tests for `SimpleGillespieSimulator` in the `ionerdss` package.

Tested Components
-----------------
- `convert_to_microscopic_rate_constants`:
    Converts macroscopic rate constants to molecule-based (microscopic) rates
    based on Avogadro's number and reaction order.

- `calculate_propensity`:
    Calculates the reaction propensities for a set of discrete molecular species.

- `gillespie_simulation`:
    Runs the full stochastic simulation algorithm (SSA) using the Gillespie direct method.

Test Model
----------
The test reactions represent the following system:
    R1: A + B → C        with macroscopic rate k1 = 1.0e6
    R2: 2C   → A + 2B    with macroscopic rate k2 = 5.0

Initial species concentrations:
    y = [A, B, C] = [10, 5, 3]
Simulation volume:
    1e-18 L

Author: yying7@jh.edu
-------
Part of the ioNERDSS modeling framework.
"""

import unittest
import numpy as np
from ionerdss.gillespie_simulation import simple_gillespie as SimpleGillespie
from ionerdss.gillespie_simulation import adaptive_rates as AdaptiveRates

class TestReactionGillespie(unittest.TestCase):

    def setUp(self):
        # Example parameters for testing
        self.macroscopic_rate_constants = np.array([1.0e6, 5.0])
        self.reactant_matrix = np.array([[1, 1, 0], [0, 0, 2]])
        self.product_matrix = np.array([[0, 0, 1], [1, 2, 1]])
        self.volume = 1.0e-18 # Litre!
        self.y = np.array([10, 5, 3])
        self.sgs = SimpleGillespie

    def test_convert_to_microscopic_rate_constants(self):
        microscopic_rate_constants = self.sgs.convert_to_microscopic_rate_constants(
            self.macroscopic_rate_constants, self.reactant_matrix, self.volume
        )
        expected_result = np.array([1.66053928e+00, 1.66053928e-05])
        np.testing.assert_allclose(microscopic_rate_constants, expected_result)

    def test_calculate_propensity(self):
        microscopic_rate_constants = np.array([0.1 * 90, 0.05 * 5])
        propensities = self.sgs.calculate_propensity(self.y, self.reactant_matrix, microscopic_rate_constants)
        expected_result = np.array([450, 0.75])
        np.testing.assert_allclose(propensities, expected_result)

    def test_gillespie_simulation(self):
        max_time = 10.0
        y_init = np.array([10, 5, 3])
        
        record_interval = 10
        full_update_scheme = True

        y_record, t_record = self.sgs.gillespie_simulation(
            max_time, y_init, self.reactant_matrix, self.product_matrix,
            self.macroscopic_rate_constants, record_interval, full_update_scheme
        )

    def test_gillespie_simulation_1D_adaptive_rate(self):

        ratelist = np.array([10, 0.2])
        y_init = np.array([10, 3])
        reactant_matrix = np.array([[2, 0], [0, 1]])
        product_matrix = np.array([[0, 1], [2, 0]])
        Length = 50
        diffusion_constants = np.array([0.1, 0.05])
        sigmalist = np.array([1, 0])
        reverse_reaction_pairs = {0:1}

        max_time = 10.0
        record_interval = 0.01
        full_update_scheme = True

        rate_update_rule = AdaptiveRates.get_rule(
            AdaptiveRates.adaptive_bimolecular_rate_1D,
            Length, diffusion_constants, sigmalist, reverse_reaction_pairs
        )

        y_record, t_record = self.sgs.run_Gillespie(
            max_time, y_init, reactant_matrix, product_matrix,
            ratelist, Length, macroscopic=True, 
            rate_update_rules=rate_update_rule,
            record_interval=record_interval, 
            full_update_scheme=full_update_scheme
        )

        # You can add more assertions based on the expected behavior of your simulation

if __name__ == '__main__':
    unittest.main()

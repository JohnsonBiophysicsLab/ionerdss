"""
Unit tests for interface processing functions in `regularize_repeats.py`.

This module tests:
1. `determine_homo_dimerization`: logic for deciding whether a chain-chain interaction is a symmetric homodimer.
2. `process_interfaces_for_chain`: processing and registering interface templates between chain pairs.

These functions are used in the context of symmetry-aware coarse-graining of molecular structures, where chain-chain interfaces
are abstracted into geometric templates and reused or created based on spatial signatures.

Mocks are used for molecule and template objects, and a minimal cg_model-like structure is provided for test input.

Author: yying7@jh.edu
in
"""

import unittest
from unittest.mock import MagicMock

from ionerdss.model.pdb.regularize_repeats import (
    determine_homo_dimerization,
    process_interfaces_for_chain
)
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

class TestInterfaceProcessing(unittest.TestCase):
    
    def setUp(self):
        self.params = PDBModelHyperparameters()

    def test_determine_homo_dimerization_symmetric(self):
        signature = {"dA": 1.0, "dB": 1.01, "thetaA": 45.0, "thetaB": 45.1}
        chains_map = {"A": "mol", "B": "mol"}
        self.assertTrue(
            determine_homo_dimerization(signature, "A", "B", chains_map, params=self.params)
        )

    def test_determine_homo_dimerization_asymmetric(self):
        signature = {"dA": 1.0, "dB": 1.3, "thetaA": 45.0, "thetaB": 75.0}
        chains_map = {"A": "mol", "B": "mol"}
        self.assertFalse(
            determine_homo_dimerization(signature, "A", "B", chains_map, params=self.params)
        )

    def test_determine_homo_dimerization_different_templates(self):
        signature = {"dA": 1.0, "dB": 1.0, "thetaA": 45.0, "thetaB": 45.0}
        chains_map = {"A": "mol1", "B": "mol2"}
        self.assertFalse(
            determine_homo_dimerization(signature, "A", "B", chains_map, params=self.params)
        )

    def test_process_interfaces_for_chain_registers_interface(self):
        # Set up mocks and minimal test data
        cg_model = {
            "interfaces": [["B"], ["A"]],
            "COMs": [(0, 0, 0), (10, 0, 0)],
            "interface_coords": [[(1, 0, 0)], [(9, 0, 0)]],
            "interface_residues": [[list(range(5))], [list(range(5))]],
            "interface_energies": [[-1.0], [-1.0]]
        }

        def mock_get_chain_data(name, model, key):
            index = 0 if name == "A" else 1
            return model[key][index]

        # Patch functions used in process_interfaces_for_chain
        import ionerdss.model.pdb.regularize_repeats as rr
        rr.get_chain_data = mock_get_chain_data
        rr.get_or_create_molecule_template = lambda name, _: MagicMock(name=name)
        rr.get_or_create_molecule = lambda name, template, _, __: MagicMock(name=name, template=template)
        rr.build_signature = lambda c1, i1, c2, i2: {"dA": 1.0, "dB": 1.0, "thetaA": 45.0, "thetaB": 45.0}
        rr.invert_signature = lambda sig: sig
        rr.signature_hash = lambda sig: str(sig)
        rr.register_interfaces = MagicMock()
        rr.build_new_interface_templates = lambda *args, **kwargs: [MagicMock(name="I1"), MagicMock(name="I2")]

        # Inputs
        chains_map = {"A": "A", "B": "A"}
        molecule_list = []
        molecule_template_list = []
        interface_signatures = []
        interface_template_list = []
        interface_list = []
        binding_chains_pairs = []
        signature_to_template_map = {}

        molecule = MagicMock(name="A")

        process_interfaces_for_chain(
            chain_id="A",
            group=["A", "B"],
            j=0,
            molecule=molecule,
            cg_model=cg_model,
            chains_map=chains_map,
            molecule_list=molecule_list,
            molecule_template_list=molecule_template_list,
            interface_signatures=interface_signatures,
            params=self.params,
            interface_template_list=interface_template_list,
            interface_list=interface_list,
            binding_chains_pairs=binding_chains_pairs,
            signature_to_template_map=signature_to_template_map
        )

        # Confirm interface was registered and template added
        self.assertIn(str({'dA': 1.0, 'dB': 1.0, 'thetaA': 45.0, 'thetaB': 45.0}), signature_to_template_map)


if __name__ == "__main__":
    unittest.main()

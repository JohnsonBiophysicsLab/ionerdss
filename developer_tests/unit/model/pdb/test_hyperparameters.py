"""
Unit tests for ionerdss.model.pdb.hyperparameters

Tests the PDBModelHyperparameters class and its configuration management.

"""

import json
import unittest
from types import SimpleNamespace

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices

from ionerdss.model.pdb.api import set_hyperparameters
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

# The twelve gap scores to_dict writes for an aligner, one by one
GAP_SCORE_KEYS = [
    f"{side}_{position}_{kind}_gap_score"
    for side in ("target", "query")
    for position in ("internal", "left", "right")
    for kind in ("open", "extend")
]


class TestPDBModelHyperparameters(unittest.TestCase):
    """Test cases for PDBModelHyperparameters class."""

    def test_default_initialization(self):
        """Test default hyperparameter initialization."""
        params = PDBModelHyperparameters()

        # Check default values
        self.assertEqual(params.interface_detect_distance_cutoff, 0.9)
        self.assertEqual(params.interface_detect_n_residue_cutoff, 2)
        self.assertEqual(params.chain_grouping_rmsd_threshold, 2.0)
        self.assertEqual(params.chain_grouping_seq_threshold, 0.5)
        self.assertEqual(params.chain_grouping_matching_mode, "default")
        self.assertEqual(params.steric_clash_mode, "off")
        self.assertEqual(params.signature_precision, 6)
        self.assertEqual(params.homodimer_distance_threshold, 0.5)
        self.assertEqual(params.homodimer_angle_threshold, 0.5)

        self.assertEqual(params.is_on_sphere, False)
        self.assertEqual(params.pdb_file_format, "bioassembly1")

        # Check that custom_aligner is created
        self.assertIsInstance(params.chain_grouping_custom_aligner, PairwiseAligner)

    def test_custom_initialization(self):
        """Test hyperparameter initialization with custom values."""
        custom_aligner = PairwiseAligner()
        custom_aligner.mode = "local"

        params = PDBModelHyperparameters(
            interface_detect_distance_cutoff=0.8,
            interface_detect_n_residue_cutoff=5,
            chain_grouping_rmsd_threshold=1.5,
            chain_grouping_seq_threshold=0.8,
            chain_grouping_custom_aligner=custom_aligner,
            chain_grouping_matching_mode="sequence",
            steric_clash_mode="auto",
            signature_precision=4,
            homodimer_distance_threshold=0.2,
            homodimer_angle_threshold=0.15,

            is_on_sphere=True
        )

        # Check custom values
        self.assertEqual(params.interface_detect_distance_cutoff, 0.8)
        self.assertEqual(params.interface_detect_n_residue_cutoff, 5)
        self.assertEqual(params.chain_grouping_rmsd_threshold, 1.5)
        self.assertEqual(params.chain_grouping_seq_threshold, 0.8)
        self.assertEqual(params.chain_grouping_custom_aligner, custom_aligner)
        self.assertEqual(params.chain_grouping_matching_mode, "sequence")
        self.assertEqual(params.steric_clash_mode, "auto")
        self.assertEqual(params.signature_precision, 4)
        self.assertEqual(params.homodimer_distance_threshold, 0.2)
        self.assertEqual(params.homodimer_angle_threshold, 0.15)

        self.assertEqual(params.is_on_sphere, True)

    def test_post_init_default_aligner_creation(self):
        """Test that __post_init__ creates default aligner when None provided."""
        params = PDBModelHyperparameters(chain_grouping_custom_aligner=None)

        # Should create default aligner
        self.assertIsInstance(params.chain_grouping_custom_aligner, PairwiseAligner)
        self.assertEqual(params.chain_grouping_custom_aligner.mode, "global")
        self.assertEqual(params.chain_grouping_custom_aligner.match_score, 1.0)
        self.assertEqual(params.chain_grouping_custom_aligner.mismatch_score, 0.0)
        self.assertEqual(params.chain_grouping_custom_aligner.open_gap_score, -0.5)
        self.assertEqual(params.chain_grouping_custom_aligner.extend_gap_score, -0.5)

    def test_create_default_aligner(self):
        """Test _create_default_aligner method."""
        params = PDBModelHyperparameters()
        aligner = params._create_default_aligner()

        # Check aligner configuration
        self.assertIsInstance(aligner, PairwiseAligner)
        self.assertEqual(aligner.mode, "global")
        self.assertEqual(aligner.match_score, 1.0)
        self.assertEqual(aligner.mismatch_score, 0.0)
        self.assertEqual(aligner.open_gap_score, -0.5)
        self.assertEqual(aligner.extend_gap_score, -0.5)

    def test_to_dict_with_default_aligner(self):
        """Test to_dict method with default aligner."""
        params = PDBModelHyperparameters()
        result = params.to_dict()

        # Check basic fields
        self.assertEqual(result['interface_detect_distance_cutoff'], 0.9)
        self.assertEqual(result['interface_detect_n_residue_cutoff'], 2)
        self.assertEqual(result['chain_grouping_rmsd_threshold'], 2.0)
        self.assertEqual(result['chain_grouping_seq_threshold'], 0.5)
        self.assertEqual(result['chain_grouping_matching_mode'], "default")
        self.assertEqual(result['steric_clash_mode'], "off")

        # Check aligner serialization: every setting, each gap score on its own
        self.assertIn('chain_grouping_custom_aligner', result)
        aligner_dict = result['chain_grouping_custom_aligner']
        self.assertEqual(aligner_dict, {
            'mode': 'global',
            'match_score': 1.0,
            'mismatch_score': 0.0,
            **{key: -0.5 for key in GAP_SCORE_KEYS},
            'wildcard': None,
            'epsilon': 1e-06,
        })

    def test_to_dict_with_none_aligner(self):
        """Test to_dict method with None aligner."""
        # Create params and manually set aligner to None (bypassing __post_init__)
        params = PDBModelHyperparameters.__new__(PDBModelHyperparameters)
        params.interface_detect_distance_cutoff = 0.6
        params.interface_detect_n_residue_cutoff = 3
        params.chain_grouping_rmsd_threshold = 2.0
        params.chain_grouping_seq_threshold = 0.5
        params.chain_grouping_custom_aligner = None
        params.chain_grouping_matching_mode = "default"
        params.steric_clash_mode = "off"
        params.signature_precision = 6
        params.homodimer_distance_threshold = 0.5
        params.homodimer_angle_threshold = 0.5
        params.homotypic_detection = "auto"
        params.homotypic_detection_residue_similarity_threshold = 0.7
        params.homotypic_detection_interface_radius = 8.0
        params.homotypic_detection_interface_radius = 8.0
        params.is_on_sphere = False
        params.template_regularization_strength = 0.0
        params.generate_visualizations = True
        params.generate_nerdss_files = True
        params.nerdss_water_box = [100.0, 100.0, 100.0]
        params.predict_affinity = False
        params.adfr_path = None
        params.pdb_file_format = "bioassembly1"
        params.ode_enabled = False
        params.ode_time_span = (0.0, 10.0)
        params.ode_solver_method = "BDF"
        params.ode_atol = 1e-4
        params.ode_plot = True
        params.ode_save_csv = True
        params.ode_initial_concentrations = None
        params.count_transition = False
        params.transition_matrix_size = 500
        params.transition_write = None
        from ionerdss.model.components.units import Units
        params.units = Units()

        result = params.to_dict()
        self.assertIsNone(result['chain_grouping_custom_aligner'])

    def test_to_dict_with_custom_aligner(self):
        """Test to_dict method with custom aligner."""
        custom_aligner = PairwiseAligner()
        custom_aligner.mode = "local"
        custom_aligner.match_score = 2.0
        custom_aligner.mismatch_score = -1.0

        params = PDBModelHyperparameters(chain_grouping_custom_aligner=custom_aligner)
        result = params.to_dict()

        # Check custom aligner serialization
        aligner_dict = result['chain_grouping_custom_aligner']
        self.assertEqual(aligner_dict['mode'], 'local')
        self.assertEqual(aligner_dict['match_score'], 2.0)
        self.assertEqual(aligner_dict['mismatch_score'], -1.0)

    def test_from_dict_empty(self):
        """Test from_dict with empty dictionary."""
        params = PDBModelHyperparameters.from_dict({})

        # Should create default instance
        self.assertEqual(params.interface_detect_distance_cutoff, 0.9)
        self.assertEqual(params.interface_detect_n_residue_cutoff, 2)
        self.assertIsInstance(params.chain_grouping_custom_aligner, PairwiseAligner)

    def test_from_dict_none(self):
        """Test from_dict with None input."""
        params = PDBModelHyperparameters.from_dict(None)

        # Should create default instance
        self.assertEqual(params.interface_detect_distance_cutoff, 0.9)
        self.assertEqual(params.interface_detect_n_residue_cutoff, 2)

    def test_from_dict_basic_fields(self):
        """Test from_dict with basic field values."""
        data = {
            'interface_detect_distance_cutoff': 0.8,
            'interface_detect_n_residue_cutoff': 5,
            'chain_grouping_rmsd_threshold': 1.5,
            'chain_grouping_seq_threshold': 0.8,
            'chain_grouping_matching_mode': 'sequence',
            'steric_clash_mode': 'auto',
            'signature_precision': 4
        }

        params = PDBModelHyperparameters.from_dict(data)

        self.assertEqual(params.interface_detect_distance_cutoff, 0.8)
        self.assertEqual(params.interface_detect_n_residue_cutoff, 5)
        self.assertEqual(params.chain_grouping_rmsd_threshold, 1.5)
        self.assertEqual(params.chain_grouping_seq_threshold, 0.8)
        self.assertEqual(params.chain_grouping_matching_mode, 'sequence')
        self.assertEqual(params.steric_clash_mode, 'auto')
        self.assertEqual(params.signature_precision, 4)

    def test_from_dict_with_aligner_dict(self):
        """Test from_dict with aligner dictionary."""
        data = {
            'interface_detect_distance_cutoff': 0.7,
            'chain_grouping_custom_aligner': {
                'mode': 'local',
                'match_score': 2.0,
                'mismatch_score': -1.0,
                'open_gap_score': -1.0,
                'extend_gap_score': -0.1
            }
        }

        params = PDBModelHyperparameters.from_dict(data)

        self.assertEqual(params.interface_detect_distance_cutoff, 0.7)
        self.assertIsInstance(params.chain_grouping_custom_aligner, PairwiseAligner)
        self.assertEqual(params.chain_grouping_custom_aligner.mode, 'local')
        self.assertEqual(params.chain_grouping_custom_aligner.match_score, 2.0)
        self.assertEqual(params.chain_grouping_custom_aligner.mismatch_score, -1.0)
        self.assertEqual(params.chain_grouping_custom_aligner.open_gap_score, -1.0)
        self.assertEqual(params.chain_grouping_custom_aligner.extend_gap_score, -0.1)

    def test_from_dict_with_none_aligner(self):
        """Test from_dict with None aligner."""
        data = {
            'interface_detect_distance_cutoff': 0.7,
            'chain_grouping_custom_aligner': None
        }

        params = PDBModelHyperparameters.from_dict(data)

        self.assertEqual(params.interface_detect_distance_cutoff, 0.7)
        # Should still create default aligner due to __post_init__
        self.assertIsInstance(params.chain_grouping_custom_aligner, PairwiseAligner)

    def test_from_dict_unknown_fields(self):
        """Test from_dict ignores unknown fields."""
        data = {
            'interface_detect_distance_cutoff': 0.8,
            'unknown_field': 'should_be_ignored',
            'another_unknown': 123
        }

        params = PDBModelHyperparameters.from_dict(data)

        self.assertEqual(params.interface_detect_distance_cutoff, 0.8)
        self.assertFalse(hasattr(params, 'unknown_field'))
        self.assertFalse(hasattr(params, 'another_unknown'))

    def test_validate_valid_parameters(self):
        """Test validate method with valid parameters."""
        params = PDBModelHyperparameters()
        errors = params.validate()

        self.assertEqual(len(errors), 0)

    def test_validate_invalid_distance_cutoff(self):
        """Test validate method with invalid distance_cutoff."""
        params = PDBModelHyperparameters(interface_detect_distance_cutoff=0.0)
        errors = params.validate()

        self.assertIn("distance_cutoff must be positive", errors)

        params = PDBModelHyperparameters(interface_detect_distance_cutoff=-0.5)
        errors = params.validate()

        self.assertIn("distance_cutoff must be positive", errors)

    def test_validate_invalid_residue_cutoff(self):
        """Test validate method with invalid residue_cutoff."""
        params = PDBModelHyperparameters(interface_detect_n_residue_cutoff=0)
        errors = params.validate()

        self.assertIn("residue_cutoff must be at least 1", errors)

    def test_validate_invalid_rmsd_threshold(self):
        """Test validate method with invalid rmsd_threshold."""
        params = PDBModelHyperparameters(chain_grouping_rmsd_threshold=-1.0)
        errors = params.validate()

        self.assertIn("rmsd_threshold must be non-negative", errors)

    def test_validate_invalid_seq_threshold(self):
        """Test validate method with invalid seq_threshold."""
        params = PDBModelHyperparameters(chain_grouping_seq_threshold=-0.1)
        errors = params.validate()

        self.assertIn("seq_threshold must be between 0 and 1", errors)

        params = PDBModelHyperparameters(chain_grouping_seq_threshold=1.5)
        errors = params.validate()

        self.assertIn("seq_threshold must be between 0 and 1", errors)

    def test_validate_invalid_signature_precision(self):
        """Test validate method with invalid signature_precision."""
        params = PDBModelHyperparameters(signature_precision=-1)
        errors = params.validate()

        self.assertIn("signature_precision must be non-negative", errors)

    def test_validate_invalid_homodimer_thresholds(self):
        """Test validate method with invalid homodimer thresholds."""
        params = PDBModelHyperparameters(homodimer_distance_threshold=-0.1)
        errors = params.validate()

        self.assertIn(
            "homodimer_distance_threshold must be non-negative", errors)

        params = PDBModelHyperparameters(homodimer_angle_threshold=-0.1)
        errors = params.validate()

        self.assertIn("homodimer_angle_threshold must be non-negative", errors)

    def test_validate_multiple_errors(self):
        """Test validate method with multiple invalid parameters."""
        params = PDBModelHyperparameters(
            interface_detect_distance_cutoff=-0.5,
            interface_detect_n_residue_cutoff=0,
            chain_grouping_seq_threshold=2.0,
            signature_precision=-1
        )
        errors = params.validate()

        # Should have multiple error messages
        self.assertGreaterEqual(len(errors), 4)
        self.assertIn("distance_cutoff must be positive", errors)
        self.assertIn("residue_cutoff must be at least 1", errors)
        self.assertIn("seq_threshold must be between 0 and 1", errors)
        self.assertIn("signature_precision must be non-negative", errors)

    def test_str_representation(self):
        """Test string representation."""
        params = PDBModelHyperparameters()
        str_repr = str(params)

        self.assertIn("PDBModelHyperparameters", str_repr)
        self.assertIn("distance_cutoff=0.9", str_repr)
        self.assertIn("residue_cutoff=2", str_repr)
        self.assertIn("matching_mode='default'", str_repr)
        self.assertIn("steric_clash_mode='off'", str_repr)

    def test_repr_representation(self):
        """Test repr representation."""
        params = PDBModelHyperparameters()
        repr_str = repr(params)

        # Should be same as __str__
        self.assertEqual(repr_str, str(params))

    def test_round_trip_serialization(self):
        """Test round-trip serialization (to_dict -> from_dict)."""
        # Create params with custom values
        original_params = PDBModelHyperparameters(
            interface_detect_distance_cutoff=0.8,
            interface_detect_n_residue_cutoff=5,
            chain_grouping_rmsd_threshold=1.5,
            chain_grouping_seq_threshold=0.8,
            chain_grouping_matching_mode="sequence",
            steric_clash_mode="auto",
            signature_precision=4,
            homodimer_distance_threshold=0.2,
            homodimer_angle_threshold=0.15
        )

        # Serialize to dict
        params_dict = original_params.to_dict()

        # Deserialize from dict
        restored_params = PDBModelHyperparameters.from_dict(params_dict)

        # Check that all values are preserved
        self.assertEqual(restored_params.interface_detect_distance_cutoff, 0.8)
        self.assertEqual(restored_params.interface_detect_n_residue_cutoff, 5)
        self.assertEqual(restored_params.chain_grouping_rmsd_threshold, 1.5)
        self.assertEqual(restored_params.chain_grouping_seq_threshold, 0.8)
        self.assertEqual(restored_params.chain_grouping_matching_mode, "sequence")
        self.assertEqual(restored_params.steric_clash_mode, "auto")
        self.assertEqual(restored_params.signature_precision, 4)
        self.assertEqual(restored_params.homodimer_distance_threshold, 0.2)
        self.assertEqual(restored_params.homodimer_angle_threshold, 0.15)

        # Check aligner parameters
        self.assertEqual(restored_params.chain_grouping_custom_aligner.mode,
                         original_params.chain_grouping_custom_aligner.mode)
        self.assertEqual(restored_params.chain_grouping_custom_aligner.match_score,
                         original_params.chain_grouping_custom_aligner.match_score)

    def test_literal_type_constraints(self):
        """Test that literal type constraints are properly defined."""
        # Test valid matching_mode values
        for mode in ["default", "sequence", "structure"]:
            params = PDBModelHyperparameters(chain_grouping_matching_mode=mode)
            self.assertEqual(params.chain_grouping_matching_mode, mode)

        # Test valid steric_clash_mode values
        for mode in ["off", "auto", "custom"]:
            params = PDBModelHyperparameters(steric_clash_mode=mode)
            self.assertEqual(params.steric_clash_mode, mode)

    def test_sphere_regularization_parameters(self):
        """Test sphere regularization parameters."""
        params = PDBModelHyperparameters(
            is_on_sphere=True
        )

        self.assertEqual(params.is_on_sphere, True)

    def test_to_dict_rejects_aligner_of_another_type(self):
        """Test to_dict raises for an aligner that is not a PairwiseAligner."""
        # It used to write default scores in place of the ones it could not read
        params = PDBModelHyperparameters(chain_grouping_custom_aligner="None")

        with self.assertRaisesRegex(TypeError, "must be a Bio.Align.PairwiseAligner"):
            params.to_dict()

    def test_from_dict_invalid_aligner_params(self):
        """Test from_dict skips aligner keys that are not settings, with a warning."""
        data = {
            'chain_grouping_custom_aligner': {
                'mode': 'local',
                'invalid_param': 'should_be_ignored',
                'match_score': 2.0
            }
        }

        with self.assertWarnsRegex(UserWarning, "ignored 'invalid_param'"):
            params = PDBModelHyperparameters.from_dict(data)

        # Should create aligner and set valid parameters
        self.assertEqual(params.chain_grouping_custom_aligner.mode, 'local')
        self.assertEqual(params.chain_grouping_custom_aligner.match_score, 2.0)
        # Invalid parameter should be ignored (no error)
        self.assertFalse(hasattr(params.chain_grouping_custom_aligner, 'invalid_param'))


class TestPDBModelHyperparametersIntegration(unittest.TestCase):
    """Integration tests for PDBModelHyperparameters."""

    def test_realistic_configuration(self):
        """Test with realistic configuration values."""
        # High-resolution structure parameters
        high_res_params = PDBModelHyperparameters(
            interface_detect_distance_cutoff=0.5,  # Tight contacts
            interface_detect_n_residue_cutoff=5,     # Substantial interfaces
            chain_grouping_rmsd_threshold=1.0,   # Strict structural similarity
            chain_grouping_seq_threshold=0.9,    # High sequence identity
            chain_grouping_matching_mode="structure",
            steric_clash_mode="auto",
            signature_precision=8  # High precision
        )

        errors = high_res_params.validate()
        self.assertEqual(len(errors), 0)

        # Low-resolution structure parameters
        low_res_params = PDBModelHyperparameters(
            interface_detect_distance_cutoff=1.2,  # Loose contacts
            interface_detect_n_residue_cutoff=3,     # Minimal interfaces
            chain_grouping_rmsd_threshold=5.0,   # Permissive structural similarity
            chain_grouping_seq_threshold=0.3,    # Low sequence identity
            chain_grouping_matching_mode="default"
        )

        errors = low_res_params.validate()
        self.assertEqual(len(errors), 0)

    def test_configuration_serialization_workflow(self):
        """Test complete configuration save/load workflow."""
        # Create configuration
        config = PDBModelHyperparameters(
            interface_detect_distance_cutoff=0.7,
            interface_detect_n_residue_cutoff=4,
            chain_grouping_matching_mode="sequence",
            steric_clash_mode="auto"
        )

        # Serialize
        config_dict = config.to_dict()

        # Simulate saving/loading (e.g., JSON)
        import json
        json_str = json.dumps(config_dict)
        loaded_dict = json.loads(json_str)

        # Deserialize
        restored_config = PDBModelHyperparameters.from_dict(loaded_dict)

        # Verify configuration is preserved
        self.assertEqual(restored_config.interface_detect_distance_cutoff, 0.7)
        self.assertEqual(restored_config.interface_detect_n_residue_cutoff, 4)
        self.assertEqual(restored_config.chain_grouping_matching_mode, "sequence")
        self.assertEqual(restored_config.steric_clash_mode, "auto")


class TestAlignerSerialization(unittest.TestCase):
    """The custom aligner must come back from a saved configuration unchanged."""

    @staticmethod
    def json_round_trip(aligner):
        """Save and load an aligner as export/import_hyperparameters do.

        Returns the aligner's JSON entry and the aligner from_dict rebuilds from it.
        """
        data = json.loads(json.dumps(
            PDBModelHyperparameters(chain_grouping_custom_aligner=aligner).to_dict()))
        restored = PDBModelHyperparameters.from_dict(data).chain_grouping_custom_aligner
        return data['chain_grouping_custom_aligner'], restored

    def assertScoresAlike(self, expected, restored, pairs):
        for target, query in pairs:
            self.assertEqual(restored.score(target, query),
                             expected.score(target, query), (target, query))

    def test_blosum62_aligner(self):
        """Test a BLOSUM62 aligner is saved by matrix name and scores the same."""
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10.0
        aligner.extend_gap_score = -0.5

        entry, restored = self.json_round_trip(aligner)

        self.assertEqual(entry['substitution_matrix'], 'BLOSUM62')
        self.assertNotIn('match_score', entry)
        self.assertEqual(restored.substitution_matrix.alphabet,
                         aligner.substitution_matrix.alphabet)
        # as plain arrays: the Array subclass cannot take the boolean masks
        # assert_array_equal indexes with
        np.testing.assert_array_equal(np.asarray(restored.substitution_matrix),
                                      np.asarray(aligner.substitution_matrix))
        self.assertScoresAlike(aligner, restored,
                               [("KEVLAPWQ", "EVLAKWQ"), ("MKTAYIAKQR", "KTAYIQR")])

    def test_free_end_gap_aligner(self):
        """Test an aligner with free end gaps and penalized internal gaps."""
        aligner = PairwiseAligner()
        aligner.end_gap_score = 0.0
        aligner.internal_gap_score = -1.0

        entry, restored = self.json_round_trip(aligner)

        self.assertEqual(entry['target_left_open_gap_score'], 0.0)
        self.assertEqual(entry['query_internal_extend_gap_score'], -1.0)
        # The first pair needs end gaps and the second an internal one, so losing
        # either kind of score changes a score
        self.assertScoresAlike(aligner, restored,
                               [("GGGKEVLGGG", "KEVL"), ("KEVLA", "KVLA")])

    def test_custom_matrix_with_codon_alphabet(self):
        """Test a matrix that is not a standard one is saved as alphabet and values."""
        codons = ("AAA", "AAG", "GAA", "GAG")  # Lys, Lys, Glu, Glu
        values = np.array([[5.0, 2.0, -1.0, -1.0],
                           [2.0, 5.0, -1.0, -1.0],
                           [-1.0, -1.0, 5.0, 2.0],
                           [-1.0, -1.0, 2.0, 5.0]])
        aligner = PairwiseAligner()
        aligner.substitution_matrix = substitution_matrices.Array(codons, dims=2, data=values)
        aligner.open_gap_score = -3.0
        aligner.extend_gap_score = -1.0

        entry, restored = self.json_round_trip(aligner)

        self.assertEqual(entry['substitution_matrix'],
                         {'alphabet': list(codons), 'values': values.tolist()})
        self.assertEqual(restored.substitution_matrix.alphabet, codons)
        target, query = ["AAA", "GAA", "GAG", "AAG"], ["AAG", "GAG", "AAG"]
        self.assertEqual(restored.score(target, query), aligner.score(target, query))

    def test_legacy_five_key_aligner_still_loads(self):
        """Test the aligner entry ionerdss 2.2.5 and earlier exported still loads."""
        data = json.loads(
            '{"chain_grouping_custom_aligner": {"mode": "local", "match_score": 2.0, '
            '"mismatch_score": -1.0, "open_gap_score": -2.0, "extend_gap_score": -0.5}}')

        restored = PDBModelHyperparameters.from_dict(data).chain_grouping_custom_aligner

        expected = PairwiseAligner()
        expected.mode = "local"
        expected.match_score = 2.0
        expected.mismatch_score = -1.0
        expected.open_gap_score = -2.0
        expected.extend_gap_score = -0.5
        self.assertScoresAlike(expected, restored,
                               [("KEVLAPWQ", "EVLAKWQ"), ("GGGKEVLGGG", "KEVL")])

    def test_legacy_matrix_aligner_entry_asks_for_the_matrix(self):
        """Test the entry 2.2.5 wrote for a matrix aligner raises a clear error.

        That version wrote match_score and mismatch_score as null and left the
        matrix out, so the file does not say which matrix was used.
        """
        settings = {'mode': 'global', 'match_score': None, 'mismatch_score': None,
                    'open_gap_score': -10.0, 'extend_gap_score': -0.5}

        with self.assertRaisesRegex(ValueError, '"substitution_matrix": "BLOSUM62"'):
            PDBModelHyperparameters.from_dict({'chain_grouping_custom_aligner': settings})

        # Naming the matrix, as the message says, makes the entry load
        settings['substitution_matrix'] = 'BLOSUM62'
        restored = PDBModelHyperparameters.from_dict(
            {'chain_grouping_custom_aligner': settings}).chain_grouping_custom_aligner
        np.testing.assert_array_equal(np.asarray(restored.substitution_matrix),
                                      np.asarray(substitution_matrices.load("BLOSUM62")))

    def test_single_gap_score_overrides_combined_one_in_any_order(self):
        """Test a gap score of one side and end wins over open_gap_score wherever it is."""
        for settings in ({'open_gap_score': -2.0, 'target_left_open_gap_score': 0.0},
                         {'target_left_open_gap_score': 0.0, 'open_gap_score': -2.0}):
            with self.subTest(order=list(settings)):
                params = PDBModelHyperparameters.from_dict(
                    {'chain_grouping_custom_aligner': settings})

                entry = params.to_dict()['chain_grouping_custom_aligner']
                self.assertEqual(entry['target_left_open_gap_score'], 0.0)
                self.assertEqual(entry['target_right_open_gap_score'], -2.0)

    def test_gap_score_function_cannot_be_serialized(self):
        """Test to_dict raises a clear error for an aligner scoring gaps with a function."""
        aligner = PairwiseAligner()
        aligner.gap_score = lambda start, length: -2.0 - 0.5 * length
        params = PDBModelHyperparameters(chain_grouping_custom_aligner=aligner)

        with self.assertRaisesRegex(ValueError, "scores gaps with a function"):
            params.to_dict()

    def test_unknown_substitution_matrix_name(self):
        """Test from_dict only loads matrices Bio.Align.substitution_matrices provides."""
        data = {'chain_grouping_custom_aligner': {'substitution_matrix': 'BLOSUM63'}}

        with self.assertRaisesRegex(ValueError, "unknown substitution matrix 'BLOSUM63'"):
            PDBModelHyperparameters.from_dict(data)

    def test_from_dict_keeps_aligner_object(self):
        """Test from_dict keeps a PairwiseAligner given in place of its settings."""
        aligner = PairwiseAligner()

        params = PDBModelHyperparameters.from_dict({'chain_grouping_custom_aligner': aligner})

        self.assertIs(params.chain_grouping_custom_aligner, aligner)

    def test_from_dict_rejects_aligner_entry_of_another_type(self):
        """Test from_dict raises for an aligner entry it cannot turn into an aligner."""
        data = {'chain_grouping_custom_aligner': 'BLOSUM62'}

        with self.assertRaisesRegex(TypeError, "dictionary of aligner settings"):
            PDBModelHyperparameters.from_dict(data)

    def test_set_hyperparameters_keeps_custom_aligner(self):
        """Test set_hyperparameters on a configured builder keeps the aligner it is given."""
        # It used to replace the aligner with a default PairwiseAligner
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 3.0
        builder = SimpleNamespace(hyperparams=PDBModelHyperparameters())

        set_hyperparameters(builder, chain_grouping_custom_aligner=aligner,
                            interface_detect_distance_cutoff=0.8)

        self.assertIs(builder.hyperparams.chain_grouping_custom_aligner, aligner)
        self.assertEqual(builder.hyperparams.interface_detect_distance_cutoff, 0.8)


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

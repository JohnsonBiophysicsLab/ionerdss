"""
Unit tests for ionerdss.model.pdb.hyperparameters

Tests the PDBModelHyperparameters class and its configuration management.

"""

import unittest
from unittest.mock import Mock
from Bio.Align import PairwiseAligner

from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters


class TestPDBModelHyperparameters(unittest.TestCase):
    """Test cases for PDBModelHyperparameters class."""

    def test_default_initialization(self):
        """Test default hyperparameter initialization."""
        params = PDBModelHyperparameters()

        # Check default values
        self.assertEqual(params.interface_detect_distance_cutoff, 0.9)
        self.assertEqual(params.interface_detect_n_residue_cutoff, 3)
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
        self.assertEqual(result['interface_detect_n_residue_cutoff'], 3)
        self.assertEqual(result['chain_grouping_rmsd_threshold'], 2.0)
        self.assertEqual(result['chain_grouping_seq_threshold'], 0.5)
        self.assertEqual(result['chain_grouping_matching_mode'], "default")
        self.assertEqual(result['steric_clash_mode'], "off")

        # Check aligner serialization
        self.assertIn('chain_grouping_custom_aligner', result)
        aligner_dict = result['chain_grouping_custom_aligner']
        self.assertEqual(aligner_dict['mode'], 'global')
        self.assertEqual(aligner_dict['match_score'], 1.0)
        self.assertEqual(aligner_dict['mismatch_score'], 0.0)
        self.assertEqual(aligner_dict['open_gap_score'], -0.5)
        self.assertEqual(aligner_dict['extend_gap_score'], -0.5)

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
        self.assertEqual(params.interface_detect_n_residue_cutoff, 3)
        self.assertIsInstance(params.chain_grouping_custom_aligner, PairwiseAligner)

    def test_from_dict_none(self):
        """Test from_dict with None input."""
        params = PDBModelHyperparameters.from_dict(None)

        # Should create default instance
        self.assertEqual(params.interface_detect_distance_cutoff, 0.9)
        self.assertEqual(params.interface_detect_n_residue_cutoff, 3)

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
        self.assertIn("residue_cutoff=3", str_repr)
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

    def test_aligner_parameter_robustness(self):
        """Test robustness of aligner parameter handling."""
        # Test with aligner missing some attributes
        mock_aligner = Mock(spec=PairwiseAligner)
        mock_aligner.mode = "local"
        # Deliberately don't set other attributes to test getattr defaults

        # Remove the automatic Mock creation for missing attributes
        del mock_aligner.match_score
        del mock_aligner.mismatch_score
        del mock_aligner.open_gap_score
        del mock_aligner.extend_gap_score

        params = PDBModelHyperparameters(chain_grouping_custom_aligner=mock_aligner)
        params_dict = params.to_dict()

        # Should handle missing attributes gracefully
        aligner_dict = params_dict['chain_grouping_custom_aligner']
        self.assertEqual(aligner_dict['mode'], 'local')
        # Missing attributes should get default values from getattr
        # Default from getattr
        self.assertEqual(aligner_dict['match_score'], 1.0)
        # Default from getattr
        self.assertEqual(aligner_dict['mismatch_score'], 0.0)
        # Default from getattr
        self.assertEqual(aligner_dict['open_gap_score'], -0.5)
        # Default from getattr
        self.assertEqual(aligner_dict['extend_gap_score'], -0.5)

    def test_from_dict_invalid_aligner_params(self):
        """Test from_dict with invalid aligner parameters."""
        data = {
            'chain_grouping_custom_aligner': {
                'mode': 'local',
                'invalid_param': 'should_be_ignored',
                'match_score': 2.0
            }
        }

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


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

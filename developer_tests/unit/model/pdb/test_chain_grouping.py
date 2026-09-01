"""
Unit tests for ionerdss.model.pdb.chain_grouping

Tests the ChainGrouper class and its various grouping strategies.
"""

import unittest
from pathlib import Path
from unittest.mock import Mock, patch
import numpy as np

from ionerdss.model.pdb.chain_grouping import ChainGroup, ChainGrouper
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters

# Multi-chain assembly fixture: 8 copies of the same chain.
MULTI_CHAIN_CIF = Path(__file__).resolve().parents[3] / "data" / "test_6BNO.cif"


class TestChainGroup(unittest.TestCase):
    """Test cases for ChainGroup class."""

    def test_chain_group_creation(self):
        """Test ChainGroup initialization."""
        group = ChainGroup("A", ["A", "B", "C"], "header")

        self.assertEqual(group.representative, "A")
        self.assertEqual(group.members, ["A", "B", "C"])
        self.assertEqual(group.grouping_method, "header")
        self.assertEqual(len(group), 3)

    def test_chain_group_contains(self):
        """Test ChainGroup membership checking."""
        group = ChainGroup("A", ["A", "B", "C"], "header")

        self.assertIn("A", group)
        self.assertIn("B", group)
        self.assertIn("C", group)
        self.assertNotIn("D", group)

    def test_chain_group_sorting(self):
        """Test that members are sorted deterministically."""
        group = ChainGroup("A", ["C", "A", "B"], "header")
        self.assertEqual(group.members, ["A", "B", "C"])


class TestChainGrouper(unittest.TestCase):
    """Test cases for ChainGrouper class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock parser
        self.mock_parser = Mock()
        self.mock_parser.get_chain_ids.return_value = ["A", "B", "C", "D"]

        # Create mock coarse grainer
        self.mock_coarse_grainer = Mock()
        self.mock_coarse_grainer.get_chain_interfaces.return_value = []

        # Create hyperparameters
        self.hyperparams = PDBModelHyperparameters(
            chain_grouping_matching_mode="default",
            chain_grouping_seq_threshold=0.8,
            chain_grouping_rmsd_threshold=2.0
        )

        # Mock aligner
        mock_aligner = Mock()
        mock_alignment = Mock()
        mock_alignment.score = 80
        mock_aligner.align.return_value = mock_alignment
        self.hyperparams.chain_grouping_custom_aligner = mock_aligner

    def test_header_based_grouping_success(self):
        """Test successful header-based grouping."""
        # Setup mock data
        self.mock_parser.get_strand_ids.return_value = {
            "1": ["A", "B"],
            "2": ["C", "D"]
        }

        # Create grouper
        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        # Verify groups were created
        self.assertEqual(len(grouper.groups), 2)

        # Check first group
        group1 = next(g for g in grouper.groups if g.representative == "A")
        self.assertEqual(group1.members, ["A", "B"])
        self.assertEqual(group1.grouping_method, "header")

        # Check second group
        group2 = next(g for g in grouper.groups if g.representative == "C")
        self.assertEqual(group2.members, ["C", "D"])
        self.assertEqual(group2.grouping_method, "header")

        # Check chain to group mapping
        self.assertEqual(grouper.chain_to_group["A"], "A")
        self.assertEqual(grouper.chain_to_group["B"], "A")
        self.assertEqual(grouper.chain_to_group["C"], "C")
        self.assertEqual(grouper.chain_to_group["D"], "C")

    def test_header_based_grouping_fallback(self):
        """Test fallback to sequence grouping when header fails."""
        # Setup mock data - no strand IDs available
        self.mock_parser.get_strand_ids.return_value = {}

        # Mock chain data for sequence grouping
        chain_data = {
            "A": {"sequence": "HCGK"},
            "B": {"sequence": "HCGK"},
            "C": {"sequence": "KGCH"},
            "D": {"sequence": "KGCH"}
        }
        self.mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]

        # Create grouper
        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        # Should have fallen back to sequence grouping
        self.assertTrue(
            any(g.grouping_method == "sequence" for g in grouper.groups))

    def test_sequence_based_grouping(self):
        """Test sequence-based grouping."""
        # Set hyperparameters for sequence grouping
        self.hyperparams.chain_grouping_matching_mode = "sequence"

        # Mock chain data
        chain_data = {
            "A": {"sequence": "HCGKHCGK"},  # 8 chars
            "B": {"sequence": "HCGKHCGK"},  # identical
            "C": {"sequence": "KCGHKCGH"},  # 8 chars, different
            "D": {"sequence": "KCGHKCGH"}   # identical to C
        }
        self.mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]

        # Mock aligner to return high score for identical sequences
        def mock_align(seq1, seq2):
            mock_alignment = Mock()
            if seq1 == seq2:
                mock_alignment.score = len(seq1)  # Perfect match
            else:
                mock_alignment.score = 0  # No match
            return mock_alignment

        self.hyperparams.chain_grouping_custom_aligner.align.side_effect = mock_align

        # Create grouper
        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        # Should have 2 groups based on sequence similarity
        self.assertEqual(len(grouper.groups), 2)

        # Check that identical sequences are grouped together
        groups_by_rep = {g.representative: g.members for g in grouper.groups}

        # A and B should be in one group, C and D in another
        if "A" in groups_by_rep:
            self.assertIn("B", groups_by_rep["A"])
        if "C" in groups_by_rep:
            self.assertIn("D", groups_by_rep["C"])

    def test_structure_based_grouping(self):
        """Test structure-based grouping."""
        # Set hyperparameters for structure grouping
        self.hyperparams.chain_grouping_matching_mode = "structure"
        self.hyperparams.chain_grouping_rmsd_threshold = 1.0

        # Mock chain data with coordinates
        coords_similar = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        coords_different = np.array([[10, 10, 10], [11, 10, 10], [12, 10, 10]])

        chain_data = {
            "A": {"ca_coords": coords_similar},
            "B": {"ca_coords": coords_similar + 0.1},  # Very similar
            "C": {"ca_coords": coords_different},
            "D": {"ca_coords": coords_different + 0.1}  # Very similar to C
        }
        self.mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]

        # Create grouper
        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            with patch.object(ChainGrouper, '_are_structures_similar_coords') as mock_similar:
                # Mock similarity function
                def mock_similarity(coords1, coords2):
                    # Similar if both are "similar" type or both are "different" type
                    is_coords1_similar = np.allclose(
                        coords1, coords_similar, atol=1.0)
                    is_coords2_similar = np.allclose(
                        coords2, coords_similar, atol=1.0)
                    return is_coords1_similar == is_coords2_similar

                mock_similar.side_effect = mock_similarity

                grouper = ChainGrouper(
                    self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        # Should have groups based on structural similarity
        self.assertTrue(
            any(g.grouping_method == "structure" for g in grouper.groups))

    def test_structure_similarity_coords(self):
        """Test _are_structures_similar_coords method."""
        # Create grouper without running full initialization
        self.mock_parser.get_strand_ids.return_value = {}  # Force fallback to sequence

        # Mock chain data for sequence grouping (to avoid header grouping)
        chain_data = {
            "A": {"sequence": "HCGK"},
            "B": {"sequence": "HCGK"},
            "C": {"sequence": "KCGH"},
            "D": {"sequence": "KCGH"}
        }
        self.mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]

        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        # Test identical coordinates (RMSD 0, below the 2.0 threshold)
        coords1 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        coords2 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])

        result = grouper._are_structures_similar_coords(coords1, coords2)
        self.assertTrue(result)

        # Test coordinates far beyond the threshold once superimposed
        coords_bent = np.array([[0, 0, 0], [1, 0, 0], [0, 10, 0]])
        result = grouper._are_structures_similar_coords(coords1, coords_bent)
        self.assertFalse(result)

        # Test different length coordinates
        coords3 = np.array([[0, 0, 0], [1, 0, 0]])  # Different length
        result = grouper._are_structures_similar_coords(coords1, coords3)
        self.assertFalse(result)

        # Test empty coordinates
        coords_empty = np.array([]).reshape(0, 3)
        result = grouper._are_structures_similar_coords(coords1, coords_empty)
        self.assertFalse(result)

    def test_structure_similarity_short_chains(self):
        """Test structure similarity for very short chains."""
        # Create grouper without running full initialization
        self.mock_parser.get_strand_ids.return_value = {}  # Force fallback to sequence

        # Mock chain data for sequence grouping
        chain_data = {
            "A": {"sequence": "HCGK"},
            "B": {"sequence": "HCGK"}
        }
        self.mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]
        self.mock_parser.get_chain_ids.return_value = ["A", "B"]

        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        # Test short chains
        coords1 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
        coords2 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])

        result = grouper._are_structures_similar_coords(coords1, coords2)
        # Should use distance-based comparison for short chains
        self.assertTrue(isinstance(result, bool))

    def test_ensure_all_chains_grouped(self):
        """Test that all chains are assigned to groups."""
        # Create grouper with some ungrouped chains
        grouper = ChainGrouper.__new__(ChainGrouper)  # Skip __init__
        grouper.parser = self.mock_parser
        grouper.groups = []
        grouper.chain_to_group = {"A": "A", "B": "A"}  # Only A and B grouped

        # Call the method
        grouper._ensure_all_chains_grouped()

        # Check that C and D got singleton groups
        self.assertIn("C", grouper.chain_to_group)
        self.assertIn("D", grouper.chain_to_group)
        self.assertEqual(grouper.chain_to_group["C"], "C")
        self.assertEqual(grouper.chain_to_group["D"], "D")

        # Check that singleton groups were created
        singleton_groups = [
            g for g in grouper.groups if g.grouping_method == "singleton"]
        self.assertEqual(len(singleton_groups), 2)

    def test_get_methods(self):
        """Test getter methods."""
        # Setup a simple grouper
        self.mock_parser.get_strand_ids.return_value = {
            "1": ["A", "B"], "2": ["C"]}

        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        # Test get_groups
        groups = grouper.get_groups()
        self.assertIsInstance(groups, list)
        self.assertTrue(len(groups) > 0)

        # Test get_group_for_chain
        group_a = grouper.get_group_for_chain("A")
        self.assertIsNotNone(group_a)
        self.assertIn("A", group_a.members)

        # Test get_representative
        rep_a = grouper.get_representative("A")
        self.assertEqual(rep_a, "A")

        # Test with non-existent chain
        group_x = grouper.get_group_for_chain("X")
        self.assertIsNone(group_x)

        rep_x = grouper.get_representative("X")
        self.assertIsNone(rep_x)

    def test_get_summary(self):
        """Test get_summary method."""
        # Setup a simple grouper
        self.mock_parser.get_strand_ids.return_value = {"1": ["A", "B"]}

        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        summary = grouper.get_summary()

        # Check summary structure
        self.assertIn("num_groups", summary)
        self.assertIn("grouping_method", summary)
        self.assertIn("groups", summary)

        # Check summary content
        self.assertIsInstance(summary["num_groups"], int)
        self.assertEqual(summary["grouping_method"], "default")
        self.assertIsInstance(summary["groups"], list)

        # Check group details
        if summary["groups"]:
            group_info = summary["groups"][0]
            self.assertIn("representative", group_info)
            self.assertIn("members", group_info)
            self.assertIn("size", group_info)
            self.assertIn("method", group_info)

    def _quasi_equivalent_chain_data(self):
        """Four chains of one sequence: A/B in one conformation, C/D in another."""
        def helix(twist):
            # a short backbone whose shape depends on `twist`
            return [
                {"id": i, "name": name,
                 "ca_coord": np.array([i * 3.8, twist * i * i * 0.4, 0.0])}
                for i, name in enumerate("HCGKHCGK")
            ]

        return {
            "A": {"sequence": "HCGKHCGK", "residues": helix(0.0)},
            "B": {"sequence": "HCGKHCGK", "residues": helix(0.0)},
            "C": {"sequence": "HCGKHCGK", "residues": helix(1.0)},
            "D": {"sequence": "HCGKHCGK", "residues": helix(1.0)},
        }

    def test_sequence_structure_grouping_splits_conformers(self):
        """Identical sequence but different shape must land in different groups."""
        self.hyperparams.chain_grouping_matching_mode = "sequence_structure"
        self.hyperparams.chain_grouping_custom_aligner = None
        self.hyperparams.chain_grouping_rmsd_threshold = 1.0

        chain_data = self._quasi_equivalent_chain_data()
        self.mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]

        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        self.assertEqual(len(grouper.groups), 2)
        groups_by_rep = {g.representative: g.members for g in grouper.groups}
        self.assertIn("B", groups_by_rep["A"])
        self.assertIn("D", groups_by_rep["C"])
        self.assertEqual(grouper.groups[0].grouping_method, "sequence_structure")

    def test_sequence_structure_grouping_merges_when_rmsd_permits(self):
        """A loose RMSD threshold collapses the conformers back into one group."""
        self.hyperparams.chain_grouping_matching_mode = "sequence_structure"
        self.hyperparams.chain_grouping_custom_aligner = None
        self.hyperparams.chain_grouping_rmsd_threshold = 100.0

        chain_data = self._quasi_equivalent_chain_data()
        self.mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]

        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        self.assertEqual(len(grouper.groups), 1)

    def test_sequence_structure_grouping_requires_matching_sequence(self):
        """Same shape is not enough when the sequences differ."""
        self.hyperparams.chain_grouping_matching_mode = "sequence_structure"
        self.hyperparams.chain_grouping_custom_aligner = None
        self.hyperparams.chain_grouping_rmsd_threshold = 100.0

        chain_data = self._quasi_equivalent_chain_data()
        chain_data["C"]["sequence"] = "KCGHKCGH"
        chain_data["D"]["sequence"] = "KCGHKCGH"
        self.mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]

        with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
            grouper = ChainGrouper(
                self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        self.assertEqual(len(grouper.groups), 2)

    def test_aligned_ca_pairs_tolerates_unequal_lengths(self):
        """Correspondence comes from the residues, not from array position."""
        long_chain = {"residues": [
            {"id": i, "name": name, "ca_coord": np.array([i * 3.8, 0.0, 0.0])}
            for i, name in enumerate("XYHCGKHCGK")]}
        short_chain = {"residues": [
            {"id": i, "name": name, "ca_coord": np.array([i * 3.8, 0.0, 0.0])}
            for i, name in enumerate("HCGKHCGK")]}

        coords_i, coords_j = ChainGrouper._aligned_ca_pairs(long_chain, short_chain)

        # the shared 8-residue block is matched despite the 2-residue overhang
        self.assertEqual(len(coords_i), 8)
        self.assertEqual(len(coords_j), 8)
        self.assertAlmostEqual(ChainGrouper._kabsch_rmsd(coords_i, coords_j), 0.0, places=6)

    def test_invalid_matching_mode(self):
        """Test handling of invalid matching mode."""
        self.hyperparams.chain_grouping_matching_mode = "invalid_mode"

        with self.assertRaises(ValueError) as context:
            ChainGrouper(self.mock_parser,
                         self.mock_coarse_grainer, self.hyperparams)

        self.assertIn("Unknown matching_mode", str(context.exception))

    def test_sequence_grouping_with_alignment_errors(self):
        """Test sequence grouping handles alignment errors gracefully."""
        self.hyperparams.chain_grouping_matching_mode = "sequence"

        # Mock chain data
        chain_data = {
            "A": {"sequence": "HCGK"},
            "B": {"sequence": "KCGH"}
        }
        self.mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]
        self.mock_parser.get_chain_ids.return_value = ["A", "B"]

        # Mock aligner to raise exception
        self.hyperparams.chain_grouping_custom_aligner.align.side_effect = Exception(
            "Alignment failed")

        # Should handle the error gracefully
        with patch('builtins.print'):  # Suppress warning prints
            with patch.object(ChainGrouper, '_ensure_all_chains_grouped'):
                grouper = ChainGrouper(
                    self.mock_parser, self.mock_coarse_grainer, self.hyperparams)

        # Should still create groups (likely singletons due to failed alignments)
        self.assertIsInstance(grouper.groups, list)

    def test_structure_grouping_with_superposition_errors(self):
        """Test structure grouping handles superposition errors gracefully."""
        # Create a minimal grouper instance for testing the method
        grouper = ChainGrouper.__new__(ChainGrouper)
        grouper.hyperparams = self.hyperparams

        coords1 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        coords2 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])

        # A non-converging SVD is treated as "not similar" and logged
        with patch.object(ChainGrouper, '_kabsch_rmsd',
                          side_effect=np.linalg.LinAlgError("SVD did not converge")):
            with self.assertLogs('ionerdss.model.pdb.chain_grouping', level='WARNING') as logs:
                result = grouper._are_structures_similar_coords(coords1, coords2)

        self.assertFalse(result)
        self.assertIn("Structure comparison failed", logs.output[0])

    def test_structure_grouping_does_not_swallow_unexpected_errors(self):
        """Unexpected errors must propagate instead of becoming 'not similar'."""
        grouper = ChainGrouper.__new__(ChainGrouper)
        grouper.hyperparams = self.hyperparams

        coords1 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
        coords2 = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])

        with patch.object(ChainGrouper, '_kabsch_rmsd',
                          side_effect=AttributeError("genuine bug")):
            with self.assertRaises(AttributeError):
                grouper._are_structures_similar_coords(coords1, coords2)

    def test_kabsch_rmsd_is_rotation_and_translation_invariant(self):
        """RMSD must be computed after optimal superposition, not on raw coords."""
        rng = np.random.default_rng(0)
        coords = rng.normal(size=(30, 3)) * 10.0

        # Build a proper rotation (determinant +1) and a translation
        rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(rotation) < 0:
            rotation[:, 0] *= -1
        moved = coords @ rotation.T + np.array([25.0, -13.0, 7.0])

        self.assertAlmostEqual(
            ChainGrouper._kabsch_rmsd(coords, moved), 0.0, places=6)
        self.assertAlmostEqual(
            ChainGrouper._kabsch_rmsd(coords, coords.copy()), 0.0, places=6)

    def test_kabsch_rmsd_rejects_reflections(self):
        """A mirror image is not a rotation and must not superimpose onto itself."""
        rng = np.random.default_rng(1)
        coords = rng.normal(size=(20, 3)) * 10.0
        mirrored = coords * np.array([1.0, 1.0, -1.0])

        self.assertGreater(
            ChainGrouper._kabsch_rmsd(coords, mirrored), 1.0)


class TestChainGrouperIntegration(unittest.TestCase):
    """Integration tests for ChainGrouper."""

    def test_full_grouping_workflow(self):
        """Test complete grouping workflow."""
        # Create mock objects
        mock_parser = Mock()
        mock_coarse_grainer = Mock()
        mock_coarse_grainer.get_chain_interfaces.return_value = []

        # Setup realistic data
        mock_parser.get_chain_ids.return_value = ["A", "B", "C", "D"]
        mock_parser.get_strand_ids.return_value = {
            "1": ["A", "B"],  # Entity 1: chains A, B
            "2": ["C"],       # Entity 2: chain C
            "3": ["D"]        # Entity 3: chain D
        }

        hyperparams = PDBModelHyperparameters(chain_grouping_matching_mode="default")

        # Create grouper
        grouper = ChainGrouper(mock_parser, mock_coarse_grainer, hyperparams)

        # Verify results
        self.assertEqual(len(grouper.groups), 3)  # 3 entities = 3 groups

        # Check specific groupings
        group_a = grouper.get_group_for_chain("A")
        self.assertEqual(set(group_a.members), {"A", "B"})

        group_c = grouper.get_group_for_chain("C")
        self.assertEqual(group_c.members, ["C"])

        group_d = grouper.get_group_for_chain("D")
        self.assertEqual(group_d.members, ["D"])

        # Check that all chains are mapped
        self.assertEqual(len(grouper.chain_to_group), 4)


class TestStructureGroupingOnRealStructure(unittest.TestCase):
    """Structure-mode grouping against a real multi-chain assembly."""

    @classmethod
    def setUpClass(cls):
        if not MULTI_CHAIN_CIF.exists():
            raise unittest.SkipTest(
                f"Test data file not found: {MULTI_CHAIN_CIF}")

        from ionerdss.model.pdb.parser import PDBParser

        cls.parser = PDBParser(MULTI_CHAIN_CIF)

    def _group(self, rmsd_threshold):
        coarse_grainer = Mock()
        coarse_grainer.get_chain_interfaces.return_value = []
        hyperparams = PDBModelHyperparameters(
            chain_grouping_matching_mode="structure",
            chain_grouping_rmsd_threshold=rmsd_threshold,
        )
        return ChainGrouper(self.parser, coarse_grainer, hyperparams)

    def test_identical_chains_collapse_into_one_group(self):
        """The 8 copies of the same chain must form a single structure group."""
        chain_ids = self.parser.get_chain_ids()
        self.assertGreater(len(chain_ids), 1, "fixture must be multi-chain")

        grouper = self._group(2.0)

        self.assertEqual(len(grouper.groups), 1)
        group = grouper.groups[0]
        self.assertEqual(group.grouping_method, "structure")
        self.assertEqual(group.members, sorted(chain_ids))

    def test_rmsd_threshold_is_applied(self):
        """A threshold below the real RMSD must split the chains apart."""
        # The chains superimpose to a small but non-zero RMSD, so a threshold
        # of 0 leaves every chain in its own group. If the superposition were
        # broken this test would pass for the wrong reason, which is why it is
        # paired with the grouping test above.
        grouper = self._group(0.0)

        self.assertEqual(len(grouper.groups),
                         len(self.parser.get_chain_ids()))
        for group in grouper.groups:
            self.assertEqual(len(group), 1)

    def test_pairwise_rmsd_between_copies_is_small(self):
        """Superimposing two copies of the same chain gives a near-zero RMSD."""
        chain_ids = self.parser.get_chain_ids()
        coords_a = self.parser.get_chain_data(chain_ids[0])["ca_coords"]
        coords_b = self.parser.get_chain_data(chain_ids[1])["ca_coords"]

        self.assertEqual(len(coords_a), len(coords_b))
        rmsd = ChainGrouper._kabsch_rmsd(coords_a, coords_b)
        self.assertGreater(rmsd, 0.0)
        self.assertLess(rmsd, 2.0)


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

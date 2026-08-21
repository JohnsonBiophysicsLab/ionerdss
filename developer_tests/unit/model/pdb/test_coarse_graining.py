"""
Unit tests for ionerdss.model.pdb.coarse_graining

Tests the CoarseGrainer class and its interface detection pipeline.
"""

import unittest
from unittest.mock import Mock, patch
import numpy as np

from ionerdss.model.pdb.coarse_graining import (
    InterfaceString, CoarseGrainedChain, CoarseGrainer
)
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters


class TestInterfaceString(unittest.TestCase):
    """Test cases for InterfaceString dataclass."""

    def test_interface_string_creation(self):
        """Test InterfaceString initialization."""
        coord_i = np.array([1.0, 2.0, 3.0])
        coord_j = np.array([4.0, 5.0, 6.0])
        residues_i = {1, 2, 3}
        residues_j = {4, 5, 6}

        interface = InterfaceString(
            chain_i="A",
            chain_j="B",
            coord_i=coord_i,
            coord_j=coord_j,
            residues_i=residues_i,
            residues_j=residues_j,
            residue_details_i=[],  # Empty list for test
            residue_details_j=[],  # Empty list for test
            energy=-2.5
        )

        self.assertEqual(interface.chain_i, "A")
        self.assertEqual(interface.chain_j, "B")
        np.testing.assert_array_equal(interface.coord_i, coord_i)
        np.testing.assert_array_equal(interface.coord_j, coord_j)
        self.assertEqual(interface.residues_i, residues_i)
        self.assertEqual(interface.residues_j, residues_j)
        self.assertEqual(interface.energy, -2.5)

    def test_interface_string_default_energy(self):
        """Test InterfaceString with default energy."""
        interface = InterfaceString(
            chain_i="A",
            chain_j="B",
            coord_i=np.array([0, 0, 0]),
            coord_j=np.array([1, 1, 1]),
            residues_i={1},
            residues_j={2},
            residue_details_i=[],  # Empty list for test
            residue_details_j=[]   # Empty list for test
        )

        self.assertEqual(interface.energy, -1.0)


class TestCoarseGrainedChain(unittest.TestCase):
    """Test cases for CoarseGrainedChain dataclass."""

    def test_coarse_grained_chain_creation(self):
        """Test CoarseGrainedChain initialization."""
        com = np.array([1.0, 2.0, 3.0])
        bbox_min = np.array([0.0, 0.0, 0.0])
        bbox_max = np.array([2.0, 4.0, 6.0])

        chain = CoarseGrainedChain(
            chain_id="A",
            com=com,
            radius=5.0,
            sequence="HCGK",
            bbox_min=bbox_min,
            bbox_max=bbox_max
        )

        self.assertEqual(chain.chain_id, "A")
        np.testing.assert_array_equal(chain.com, com)
        self.assertEqual(chain.radius, 5.0)
        self.assertEqual(chain.sequence, "HCGK")
        np.testing.assert_array_equal(chain.bbox_min, bbox_min)
        np.testing.assert_array_equal(chain.bbox_max, bbox_max)


class TestCoarseGrainer(unittest.TestCase):
    """Test cases for CoarseGrainer class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock parser
        self.mock_parser = Mock()

        # Create hyperparameters
        self.hyperparams = PDBModelHyperparameters(
            interface_detect_distance_cutoff=5.0,  # nm
            interface_detect_n_residue_cutoff=3
        )

        # Mock chain data
        self.chain_data = {
            "A": {
                "com": np.array([0.0, 0.0, 0.0]),
                "radius": 2.0,
                "sequence": "HCGK",
                "bbox_min": np.array([-1.0, -1.0, -1.0]),
                "bbox_max": np.array([1.0, 1.0, 1.0]),
                "ca_coords": np.array([
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [2.0, 0.0, 0.0],
                    [3.0, 0.0, 0.0]
                ]),
                "residues": [
                    {"id": 1, "name": "HIS"}, {"id": 2, "name": "CYS"}, {"id": 3, "name": "GLY"}, {"id": 4, "name": "LYS"}
                ]
            },
            "B": {
                "com": np.array([5.0, 0.0, 0.0]),
                "radius": 2.0,
                "sequence": "TGCA",
                "bbox_min": np.array([4.0, -1.0, -1.0]),
                "bbox_max": np.array([6.0, 1.0, 1.0]),
                "ca_coords": np.array([
                    [4.0, 0.0, 0.0],
                    [5.0, 0.0, 0.0],
                    [6.0, 0.0, 0.0],
                    [7.0, 0.0, 0.0]
                ]),
                "residues": [
                    {"id": 5, "name": "THR"}, {"id": 6, "name": "GLY"}, {"id": 7, "name": "CYS"}, {"id": 8, "name": "ALA"}
                ]
            },
            "C": {
                "com": np.array([20.0, 0.0, 0.0]),
                "radius": 2.0,
                "sequence": "AAAA",
                "bbox_min": np.array([19.0, -1.0, -1.0]),
                "bbox_max": np.array([21.0, 1.0, 1.0]),
                "ca_coords": np.array([
                    [19.0, 0.0, 0.0],
                    [20.0, 0.0, 0.0],
                    [21.0, 0.0, 0.0]
                ]),
                "residues": [
                    {"id": 9, "name": "ALA"}, {"id": 10, "name": "ALA"}, {"id": 11, "name": "ALA"}
                ]
            }
        }

        # Setup mock parser methods
        self.mock_parser.get_chain_ids.return_value = ["A", "B", "C"]
        self.mock_parser.get_chain_data.side_effect = lambda x: self.chain_data[x]
        self.mock_parser.convert_distance_to_angstrom.side_effect = lambda x: x * 10.0  # nm to Å

    def test_initialization(self):
        """Test CoarseGrainer initialization and pipeline execution."""
        with patch.object(CoarseGrainer, '_run_coarse_graining') as mock_run:
            grainer = CoarseGrainer.__new__(CoarseGrainer)
            grainer.parser = self.mock_parser
            grainer.hyperparams = self.hyperparams
            grainer.chains = {}
            grainer.interfaces = []
            grainer.partner_map = {}

            # Verify attributes are set
            self.assertEqual(grainer.parser, self.mock_parser)
            self.assertEqual(grainer.hyperparams, self.hyperparams)
            self.assertEqual(grainer.chains, {})
            self.assertEqual(grainer.interfaces, [])
            self.assertEqual(grainer.partner_map, {})

    def test_initialize_chains(self):
        """Test _initialize_chains method."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.parser = self.mock_parser
        grainer.hyperparams = self.hyperparams
        grainer.chains = {}

        grainer._initialize_chains()

        # Check that all chains were initialized
        self.assertEqual(len(grainer.chains), 3)
        self.assertIn("A", grainer.chains)
        self.assertIn("B", grainer.chains)
        self.assertIn("C", grainer.chains)

        # Check chain A details
        chain_a = grainer.chains["A"]
        self.assertEqual(chain_a.chain_id, "A")
        self.assertEqual(chain_a.sequence, "HCGK")
        self.assertEqual(chain_a.radius, 2.0)
        np.testing.assert_array_equal(chain_a.com, np.array([0.0, 0.0, 0.0]))
        # Interfaces are managed separately in grainer.interfaces, not in chain objects

    def test_can_chains_interact_true(self):
        """Test _can_chains_interact returns True for nearby chains."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.parser = self.mock_parser
        grainer.hyperparams = self.hyperparams
        grainer.chains = {}
        grainer._initialize_chains()

        # Chains A and B should be able to interact (close bounding boxes)
        result = grainer._can_chains_interact("A", "B")
        self.assertTrue(result)

    def test_can_chains_interact_false(self):
        """Test _can_chains_interact returns False for distant chains."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.parser = self.mock_parser
        grainer.hyperparams = self.hyperparams
        grainer.chains = {}
        grainer._initialize_chains()

        # Modify chain C to be much farther away
        # Current cutoff is 5.0 nm = 50.0 Å after conversion
        # Make sure the gap is larger than the cutoff
        grainer.chains["C"].bbox_min = np.array(
            [100.0, -1.0, -1.0])  # Much farther
        grainer.chains["C"].bbox_max = np.array([102.0, 1.0, 1.0])

        # Chains A and C should not interact (distant bounding boxes)
        result = grainer._can_chains_interact("A", "C")
        self.assertFalse(result)

    def test_detect_interface_valid(self):
        """Test _detect_interface with valid interface."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.parser = self.mock_parser
        grainer.hyperparams = self.hyperparams

        # Mock KDTree to return neighbors
        with patch('ionerdss.model.pdb.coarse_graining.KDTree') as mock_kdtree_class:
            mock_tree = Mock()
            # Simulate that first 3 residues of A have neighbors in B
            mock_tree.query_ball_point.return_value = [
                [0, 1],  # First residue has neighbors
                [1, 2],  # Second residue has neighbors
                [2, 3],  # Third residue has neighbors
                []       # Fourth residue has no neighbors
            ]
            mock_kdtree_class.return_value = mock_tree

            interfaces = grainer._detect_interface("A", "B")

            self.assertTrue(interfaces)
            interface = interfaces[0]
            self.assertEqual(interface.chain_i, "A")
            self.assertEqual(interface.chain_j, "B")
            self.assertEqual(len(interface.residues_i), 3)  # First 3 residues
            # All B residues contacted
            self.assertEqual(len(interface.residues_j), 4)

    def test_detect_interface_insufficient_contacts(self):
        """Test _detect_interface with insufficient contacts."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.parser = self.mock_parser
        grainer.hyperparams = self.hyperparams

        # Mock KDTree to return insufficient neighbors
        with patch('ionerdss.model.pdb.coarse_graining.KDTree') as mock_kdtree_class:
            mock_tree = Mock()
            # Only 2 residues have contacts (below residue_cutoff of 3)
            mock_tree.query_ball_point.return_value = [
                [0],  # First residue has neighbor
                [1],  # Second residue has neighbor
                [],   # Third residue has no neighbors
                []    # Fourth residue has no neighbors
            ]
            mock_kdtree_class.return_value = mock_tree

            interfaces = grainer._detect_interface("A", "B")

            self.assertEqual(interfaces, [])

    def test_detect_interface_empty_coordinates(self):
        """Test _detect_interface with empty coordinates."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.parser = self.mock_parser
        grainer.hyperparams = self.hyperparams

        # Mock empty coordinates
        empty_chain_data = self.chain_data["A"].copy()
        empty_chain_data["ca_coords"] = np.array([]).reshape(0, 3)
        empty_chain_data["residues"] = []

        with patch.object(self.mock_parser, 'get_chain_data') as mock_get_data:
            mock_get_data.side_effect = lambda x: empty_chain_data if x == "A" else self.chain_data[
                x]

            interfaces = grainer._detect_interface("A", "B")

            self.assertEqual(interfaces, [])

    def test_build_partner_mapping(self):
        """Test _build_partner_mapping method."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.chains = {"A": Mock(), "B": Mock(), "C": Mock()}
        grainer.partner_map = {}

        # Create mock interfaces
        interface1 = InterfaceString(
            chain_i="A", chain_j="B",
            coord_i=np.array([0, 0, 0]), coord_j=np.array([1, 1, 1]),
            residues_i={1}, residues_j={2},
            residue_details_i=[],  # Empty list for test
            residue_details_j=[]   # Empty list for test
        )
        interface2 = InterfaceString(
            chain_i="A", chain_j="C",
            coord_i=np.array([0, 0, 0]), coord_j=np.array([2, 2, 2]),
            residues_i={1}, residues_j={3},
            residue_details_i=[],  # Empty list for test
            residue_details_j=[]   # Empty list for test
        )
        grainer.interfaces = [interface1, interface2]

        grainer._build_partner_mapping()

        # Check bidirectional mapping
        self.assertEqual(grainer.partner_map[("A", 0)], ("B", 0))
        self.assertEqual(grainer.partner_map[("B", 0)], ("A", 0))
        self.assertEqual(grainer.partner_map[("A", 1)], ("C", 0))
        self.assertEqual(grainer.partner_map[("C", 0)], ("A", 1))

    def test_detect_all_interfaces(self):
        """Test _detect_all_interfaces method."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.parser = self.mock_parser
        grainer.hyperparams = self.hyperparams
        grainer.interfaces = []
        grainer.chains = {}
        grainer._initialize_chains()

        # Mock the interface detection methods
        with patch.object(grainer, '_can_chains_interact') as mock_can_interact:
            with patch.object(grainer, '_detect_interface') as mock_detect:
                # Set up interaction possibilities
                mock_can_interact.side_effect = lambda i, j: (
                    i, j) in [("A", "B"), ("B", "C")]

                # Mock interface detection results
                mock_interface_ab = Mock(spec=InterfaceString)
                mock_interface_ab.chain_i = "A"
                mock_interface_ab.chain_j = "B"

                def detect_side_effect(chain_i, chain_j):
                    if (chain_i, chain_j) == ("A", "B"):
                        return [mock_interface_ab]
                    return []

                mock_detect.side_effect = detect_side_effect

                grainer._detect_all_interfaces()

                # Check that interfaces were detected and added
                self.assertEqual(len(grainer.interfaces), 1)
                self.assertEqual(grainer.interfaces[0], mock_interface_ab)
                # Note: Interfaces are NO LONGER stored in chain.interfaces
                # They are stored centrally in grainer.interfaces

    def test_full_pipeline(self):
        """Test complete coarse-graining pipeline."""
        # This will test the actual __init__ method
        with patch.object(CoarseGrainer, '_detect_interface') as mock_detect:
            # Mock interface detection to return a valid interface for A-B
            mock_interface = InterfaceString(
                chain_i="A", chain_j="B",
                coord_i=np.array([2.0, 0.0, 0.0]),
                coord_j=np.array([4.0, 0.0, 0.0]),
                residues_i={3, 4}, residues_j={5, 6},
                residue_details_i=[],  # Empty list for test
                residue_details_j=[]   # Empty list for test
            )

            def detect_side_effect(chain_i, chain_j):
                if (chain_i, chain_j) == ("A", "B"):
                    return [mock_interface]
                return []

            mock_detect.side_effect = detect_side_effect

            grainer = CoarseGrainer(self.mock_parser, self.hyperparams)

            # Verify chains were initialized
            self.assertEqual(len(grainer.chains), 3)

            # Verify interfaces were detected
            self.assertEqual(len(grainer.interfaces), 1)
            self.assertEqual(grainer.interfaces[0], mock_interface)

            # Verify partner mapping was built
            self.assertIn(("A", 0), grainer.partner_map)
            self.assertIn(("B", 0), grainer.partner_map)

    def test_get_coarse_grained_chains(self):
        """Test get_coarse_grained_chains method."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.chains = {"A": Mock(), "B": Mock()}

        result = grainer.get_coarse_grained_chains()

        # Should return a copy
        self.assertEqual(len(result), 2)
        self.assertIn("A", result)
        self.assertIn("B", result)
        # Verify it's a copy, not the original
        self.assertIsNot(result, grainer.chains)

    def test_get_interfaces(self):
        """Test get_interfaces method."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        mock_interface = Mock()
        grainer.interfaces = [mock_interface]

        result = grainer.get_interfaces()

        # Should return a copy
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], mock_interface)
        # Verify it's a copy, not the original
        self.assertIsNot(result, grainer.interfaces)

    def test_get_partner_mapping(self):
        """Test get_partner_mapping method."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.partner_map = {("A", 0): ("B", 0)}

        result = grainer.get_partner_mapping()

        # Should return a copy
        self.assertEqual(len(result), 1)
        self.assertEqual(result[("A", 0)], ("B", 0))
        # Verify it's a copy, not the original
        self.assertIsNot(result, grainer.partner_map)

    def test_get_summary(self):
        """Test get_summary method."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.chains = {"A": Mock(), "B": Mock()}

        mock_interface = Mock()
        mock_interface.chain_i = "A"
        mock_interface.chain_j = "B"
        mock_interface.residue_details_i = [Mock(), Mock()]  # 2 residues
        mock_interface.residue_details_j = [Mock(), Mock(), Mock()]  # 3 residues
        grainer.interfaces = [mock_interface]

        summary = grainer.get_summary()

        # Check summary structure and content
        self.assertEqual(summary["num_chains"], 2)
        self.assertEqual(summary["num_interfaces"], 1)
        self.assertEqual(set(summary["chains"]), {"A", "B"})
        self.assertEqual(summary["interface_pairs"], [("A", "B")])

    def test_interface_coordinate_calculation(self):
        """Test that interface coordinates are calculated correctly."""
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.parser = self.mock_parser
        grainer.hyperparams = self.hyperparams

        # Mock KDTree to return specific neighbors
        with patch('ionerdss.model.pdb.coarse_graining.KDTree') as mock_kdtree_class:
            mock_tree = Mock()
            # First 3 residues of A contact first 3 residues of B
            mock_tree.query_ball_point.return_value = [
                [0, 1, 2],  # First A residue contacts first 3 B residues
                [0, 1, 2],  # Second A residue contacts first 3 B residues
                [0, 1, 2],  # Third A residue contacts first 3 B residues
                []          # Fourth A residue has no contacts
            ]
            mock_kdtree_class.return_value = mock_tree

            interfaces = grainer._detect_interface("A", "B")
            self.assertTrue(interfaces)
            interface = interfaces[0]

            # Check that interface coordinates are means of contacting residues
            # Mean of first 3 A coords
            expected_coord_i = np.array([1.0, 0.0, 0.0])
            # Mean of first 3 B coords
            expected_coord_j = np.array([5.0, 0.0, 0.0])

            np.testing.assert_array_almost_equal(
                interface.coord_i, expected_coord_i)
            np.testing.assert_array_almost_equal(
                interface.coord_j, expected_coord_j)

    def test_residue_cutoff_enforcement(self):
        """Test that residue cutoff is properly enforced."""
        # Test with residue_cutoff = 3
        grainer = CoarseGrainer.__new__(CoarseGrainer)
        grainer.parser = self.mock_parser
        grainer.hyperparams = self.hyperparams

        with patch('ionerdss.model.pdb.coarse_graining.KDTree') as mock_kdtree_class:
            mock_tree = Mock()

            # Test case 1: Exactly at cutoff (should pass)
            mock_tree.query_ball_point.return_value = [
                [0, 1, 2],  # 3 contacts
                [0, 1, 2],  # 3 contacts
                [0, 1, 2],  # 3 contacts
                []          # No contacts
            ]
            mock_kdtree_class.return_value = mock_tree

            interfaces = grainer._detect_interface("A", "B")
            self.assertTrue(interfaces)

            # Test case 2: Below cutoff (should fail)
            mock_tree.query_ball_point.return_value = [
                [0, 1],     # 2 contacts
                [0, 1],     # 2 contacts
                [],         # No contacts
                []          # No contacts
            ]

            interfaces = grainer._detect_interface("A", "B")
            self.assertEqual(interfaces, [])


class TestCoarseGrainerIntegration(unittest.TestCase):
    """Integration tests for CoarseGrainer."""

    def test_realistic_scenario(self):
        """Test with a realistic protein interaction scenario."""
        # Create mock parser with realistic data
        mock_parser = Mock()

        # Two chains that should interact
        chain_data = {
            "A": {
                "com": np.array([0.0, 0.0, 0.0]),
                "radius": 15.0,
                "sequence": "MKLAVQNC",
                "bbox_min": np.array([-10.0, -10.0, -10.0]),
                "bbox_max": np.array([10.0, 10.0, 10.0]),
                "ca_coords": np.array([
                    [-5.0, 0.0, 0.0],  # Close to B
                    [-3.0, 0.0, 0.0],  # Close to B
                    [-1.0, 0.0, 0.0],  # Close to B
                    [1.0, 0.0, 0.0],   # Close to B
                    [3.0, 0.0, 0.0],   # Close to B
                    [5.0, 5.0, 5.0],   # Far from B
                    [7.0, 7.0, 7.0],   # Far from B
                    [9.0, 9.0, 9.0]    # Far from B
                ]),
                "residues": [{"id": i, "name": "ALA"} for i in range(1, 9)]
            },
            "B": {
                "com": np.array([0.0, 0.0, 0.0]),
                "radius": 15.0,
                "sequence": "TGCARLMN",
                "bbox_min": np.array([-10.0, -10.0, -10.0]),
                "bbox_max": np.array([10.0, 10.0, 10.0]),
                "ca_coords": np.array([
                    [5.0, 0.0, 0.0],   # Close to A
                    [3.0, 0.0, 0.0],   # Close to A
                    [1.0, 0.0, 0.0],   # Close to A
                    [-1.0, 0.0, 0.0],  # Close to A
                    [-3.0, 0.0, 0.0],  # Close to A
                    [-5.0, -5.0, -5.0],  # Far from A
                    [-7.0, -7.0, -7.0],  # Far from A
                    [-9.0, -9.0, -9.0]  # Far from A
                ]),
                "residues": [{"id": i, "name": "GLY"} for i in range(9, 17)]
            }
        }

        mock_parser.get_chain_ids.return_value = ["A", "B"]
        mock_parser.get_chain_data.side_effect = lambda x: chain_data[x]
        mock_parser.convert_distance_to_angstrom.side_effect = lambda x: x * 10.0

        hyperparams = PDBModelHyperparameters(
            interface_detect_distance_cutoff=1.0,  # 10 Å cutoff
            interface_detect_n_residue_cutoff=3
        )

        # Create coarse grainer
        grainer = CoarseGrainer(mock_parser, hyperparams)

        # Verify results
        self.assertEqual(len(grainer.chains), 2)
        self.assertEqual(len(grainer.interfaces), 1)

        # Check interface details
        interface = grainer.interfaces[0]
        self.assertEqual(interface.chain_i, "A")
        self.assertEqual(interface.chain_j, "B")
        self.assertGreaterEqual(len(interface.residues_i), 3)
        self.assertGreaterEqual(len(interface.residues_j), 3)

        # Check partner mapping
        self.assertIn(("A", 0), grainer.partner_map)
        self.assertIn(("B", 0), grainer.partner_map)
        self.assertEqual(grainer.partner_map[("A", 0)], ("B", 0))
        self.assertEqual(grainer.partner_map[("B", 0)], ("A", 0))


if __name__ == '__main__':
    # Run with verbose output
    unittest.main(verbosity=2)

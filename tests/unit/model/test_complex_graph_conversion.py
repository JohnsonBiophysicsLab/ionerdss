"""
Unit tests for Complex <-> NetworkX graph conversion.

Tests the bidirectional conversion between Complex objects and NetworkX graphs,
as well as the graph-based naming scheme.
"""

import unittest
import networkx as nx
from ionerdss.model.complex import Complex
from ionerdss.model.complex_to_graph import (
    complex_to_networkx,
    generate_complex_name_from_graph,
    _classify_topology
)


class TestComplexGraphConversion(unittest.TestCase):
    """Test Complex to NetworkX conversion and naming."""
    
    def setUp(self):
        """Create mock molecules and reactions for testing."""
        # Create mock molecule template
        class MockTemplate:
            def __init__(self, name):
                self.name = name
                self.expression = f"{name}_binding"
        
        # Create mock molecules
        class MockMolecule:
            def __init__(self, name, template_name):
                self.name = name
                self.my_template = MockTemplate(template_name)
        
        # Create mock reaction
        class MockReaction:
            def __init__(self, expression):
                self.my_template = MockTemplate("reaction")
                self.my_template.expression = expression
        
        self.MockMolecule = MockMolecule
        self.MockReaction = MockReaction
    
    def test_single_molecule_complex(self):
        """Test conversion of single molecule complex."""
        complex_obj = Complex()
        mol = self.MockMolecule("A1", "A")
        complex_obj.add_interaction(mol, None, None)
        
        # Convert to graph
        G = complex_to_networkx(complex_obj)
        
        # Verify graph structure
        self.assertEqual(len(G.nodes), 1)
        self.assertEqual(len(G.edges), 0)
        self.assertEqual(G.nodes[0]['type'], 'A')
        
        # Test naming
        name = generate_complex_name_from_graph(G)
        self.assertEqual(name, 'A')
    
    def test_linear_dimer(self):
        """Test conversion of linear dimer A-B."""
        complex_obj = Complex()
        mol_a = self.MockMolecule("A1", "A")
        mol_b = self.MockMolecule("B1", "B")
        reaction = self.MockReaction("A_B_binding")
        
        complex_obj.add_interaction(mol_a, mol_b, reaction)
        complex_obj.add_interaction(mol_b, mol_a, reaction)
        
        # Convert to graph
        G = complex_to_networkx(complex_obj)
        
        # Verify graph structure
        self.assertEqual(len(G.nodes), 2)
        self.assertEqual(len(G.edges), 1)
        self.assertIn('type', G.nodes[0])
        self.assertIn('type', G.edges[0, 1])
    
    def test_linear_trimer(self):
        """Test linear trimer A-A-A."""
        complex_obj = Complex()
        mol1 = self.MockMolecule("A1", "A")
        mol2 = self.MockMolecule("A2", "A")
        mol3 = self.MockMolecule("A3", "A")
        reaction = self.MockReaction("A_A_binding")
        
        complex_obj.add_interaction(mol1, mol2, reaction)
        complex_obj.add_interaction(mol2, mol1, reaction)
        complex_obj.add_interaction(mol2, mol3, reaction)
        complex_obj.add_interaction(mol3, mol2, reaction)
        
        # Convert to graph
        G = complex_to_networkx(complex_obj)
        
        # Verify graph structure
        self.assertEqual(len(G.nodes), 3)
        self.assertEqual(len(G.edges), 2)
        
        # Test topology classification
        topology = _classify_topology(G)
        self.assertEqual(topology, 'linear')
        
        # Test naming includes composition
        name = generate_complex_name_from_graph(G, use_hash=False)
        self.assertIn('A3', name)
        self.assertIn('linear', name)
    
    def test_cyclic_trimer(self):
        """Test cyclic trimer (triangle)."""
        complex_obj = Complex()
        mol1 = self.MockMolecule("A1", "A")
        mol2 = self.MockMolecule("A2", "A")
        mol3 = self.MockMolecule("A3", "A")
        reaction = self.MockReaction("A_A_binding")
        
        # Create triangle
        complex_obj.add_interaction(mol1, mol2, reaction)
        complex_obj.add_interaction(mol2, mol1, reaction)
        complex_obj.add_interaction(mol2, mol3, reaction)
        complex_obj.add_interaction(mol3, mol2, reaction)
        complex_obj.add_interaction(mol3, mol1, reaction)
        complex_obj.add_interaction(mol1, mol3, reaction)
        
        # Convert to graph
        G = complex_to_networkx(complex_obj)
        
        # Verify it's complete (triangle)
        self.assertEqual(len(G.nodes), 3)
        self.assertEqual(len(G.edges), 3)
        
        topology = _classify_topology(G)
        self.assertEqual(topology, 'complete')
    
    def test_naming_determinism(self):
        """Test that isomorphic complexes get the same name."""
        # Create two isomorphic linear trimers with different node ordering
        complex1 = Complex()
        mol1_a = self.MockMolecule("A1", "A")
        mol1_b = self.MockMolecule("A2", "A")
        mol1_c = self.MockMolecule("A3", "A")
        reaction = self.MockReaction("binding")
        
        complex1.add_interaction(mol1_a, mol1_b, reaction)
        complex1.add_interaction(mol1_b, mol1_a, reaction)
        complex1.add_interaction(mol1_b, mol1_c, reaction)
        complex1.add_interaction(mol1_c, mol1_b, reaction)
        
        complex2 = Complex()
        mol2_a = self.MockMolecule("X1", "A")
        mol2_b = self.MockMolecule("X2", "A")
        mol2_c = self.MockMolecule("X3", "A")
        
        # Different construction order
        complex2.add_interaction(mol2_c, mol2_b, reaction)
        complex2.add_interaction(mol2_b, mol2_c, reaction)
        complex2.add_interaction(mol2_b, mol2_a, reaction)
        complex2.add_interaction(mol2_a, mol2_b, reaction)
        
        # Both should have same graph-based name
        G1 = complex_to_networkx(complex1)
        G2 = complex_to_networkx(complex2)
        
        name1 = generate_complex_name_from_graph(G1)
        name2 = generate_complex_name_from_graph(G2)
        
        # Names should be identical for isomorphic structures
        self.assertEqual(name1, name2)
    
    def test_heterogeneous_complex(self):
        """Test complex with different molecule types."""
        complex_obj = Complex()
        mol_a = self.MockMolecule("A1", "A")
        mol_b = self.MockMolecule("B1", "B")
        mol_c = self.MockMolecule("C1", "C")
        reaction_ab = self.MockReaction("A_B_binding")
        reaction_bc = self.MockReaction("B_C_binding")
        
        complex_obj.add_interaction(mol_a, mol_b, reaction_ab)
        complex_obj.add_interaction(mol_b, mol_a, reaction_ab)
        complex_obj.add_interaction(mol_b, mol_c, reaction_bc)
        complex_obj.add_interaction(mol_c, mol_b, reaction_bc)
        
        G = complex_to_networkx(complex_obj)
        name = generate_complex_name_from_graph(G, use_hash=False)
        
        # Should contain all molecule types
        self.assertIn('A1', name)
        self.assertIn('B1', name)
        self.assertIn('C1', name)
        self.assertIn('linear', name)


class TestTopologyClassification(unittest.TestCase):
    """Test topology classification function."""
    
    def test_linear_topology(self):
        """Test linear path detection."""
        G = nx.Graph()
        G.add_nodes_from([(0, {'type': 'A'}), (1, {'type': 'A'}), (2, {'type': 'A'})])
        G.add_edges_from([(0, 1, {'type': 'ab'}), (1, 2, {'type': 'ab'})])
        
        topology = _classify_topology(G)
        self.assertEqual(topology, 'linear')
    
    def test_cyclic_topology(self):
        """Test cycle detection."""
        G = nx.Graph()
        G.add_nodes_from([(i, {'type': 'A'}) for i in range(4)])
        G.add_edges_from([(0, 1, {'type': 'a'}), (1, 2, {'type': 'a'}), 
                          (2, 3, {'type': 'a'}), (3, 0, {'type': 'a'})])
        
        topology = _classify_topology(G)
        self.assertEqual(topology, 'cyclic')
    
    def test_complete_topology(self):
        """Test complete graph detection."""
        G = nx.complete_graph(4)
        nx.set_node_attributes(G, 'A', 'type')
        nx.set_edge_attributes(G, 'binding', 'type')
        
        topology = _classify_topology(G)
        self.assertEqual(topology, 'complete')
    
    def test_star_topology(self):
        """Test star graph detection."""
        G = nx.star_graph(3)
        nx.set_node_attributes(G, 'A', 'type')
        nx.set_edge_attributes(G, 'binding', 'type')
        
        topology = _classify_topology(G)
        self.assertEqual(topology, 'star')


if __name__ == '__main__':
    unittest.main()

"""
Unit tests for ionerdss.model.pdb.ring_regularizer

Tests the RingRegularizer class and its sphere projection capabilities.
"""

import unittest
from unittest.mock import Mock, MagicMock
import numpy as np

from ionerdss.model.pdb.ring_regularizer import RingRegularizer
from ionerdss.model.components.system import System
from ionerdss.model.components.instances import MoleculeInstance
from ionerdss.model.components.types import MoleculeType
from ionerdss.model.pdb.file_manager import WorkspaceManager


class TestRingRegularizer(unittest.TestCase):
    """Test cases for RingRegularizer (Sphere Projection) class."""

    def setUp(self):
        """Set up test fixtures."""
        self.system = System(workspace_path="/")
        self.workspace_manager = Mock(spec=WorkspaceManager)
        self.workspace_manager.logger = Mock()
        
        # Define molecule types
        self.type_a = Mock(spec=MoleculeType)
        self.type_a.name = "ProteinA"
        
        self.type_b = Mock(spec=MoleculeType)
        self.type_b.name = "ProteinB"
        
        # Helper to create molecules
        self.mol_count = 0
        
    def _create_mol(self, m_type, pos):
        self.mol_count += 1
        mol = MoleculeInstance(
            name=f"{m_type.name}_{self.mol_count}",
            molecule_type=m_type,
            com=np.array(pos),
            norm=np.array([0,0,1]),
            ref1=np.array([1,0,0]),
            ref2=np.array([0,1,0])
        )
        self.system.molecule_instances.add(mol)
        return mol

    def test_stoichiometry_calculation(self):
        """Test _get_stoichiometry identifies most abundant species."""
        # Add 3 A's and 2 B's
        for _ in range(3):
            self._create_mol(self.type_a, [0,0,0])
        for _ in range(2):
            self._create_mol(self.type_b, [0,0,0])
            
        regularizer = RingRegularizer(self.system, self.workspace_manager)
        name, counts = regularizer._get_stoichiometry()
        
        self.assertEqual(name, "ProteinA")
        self.assertEqual(counts["ProteinA"], 3)
        self.assertEqual(counts["ProteinB"], 2)

    def test_fit_sphere(self):
        """Test _fit_sphere with ideal points."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)
        
        # Points on a sphere of radius 10 at origin
        points = []
        for i in range(10):
            vec = np.random.randn(3)
            vec /= np.linalg.norm(vec)
            points.append(vec * 10.0)
        points = np.array(points)
        
        center, radius, error = regularizer._fit_sphere(points)
        
        # Use assertAlmostEqual with delta for numpy arrays
        np.testing.assert_allclose(center, [0,0,0], atol=1e-4)
        self.assertAlmostEqual(radius, 10.0, places=4)
        self.assertAlmostEqual(error, 0.0, places=4)

    def test_regularize_projection(self):
        """Test full regularization projects points to perfect spheres."""
        # ProteinA (Ref): Radius ~10 with noise
        # ProteinB: Radius ~20 with noise
        ref_center = np.array([1.0, 2.0, 3.0]) # Offset center
        
        mols_a = []
        for _ in range(10): # Reference species needs enough points for good fit
            vec = np.random.randn(3)
            vec /= np.linalg.norm(vec)
            pos = ref_center + vec * (10.0 + np.random.uniform(-0.5, 0.5))
            mols_a.append(self._create_mol(self.type_a, pos))
            
        mols_b = []
        for _ in range(5):
             vec = np.random.randn(3)
             vec /= np.linalg.norm(vec)
             pos = ref_center + vec * (20.0 + np.random.uniform(-0.5, 0.5))
             mols_b.append(self._create_mol(self.type_b, pos))
             
        regularizer = RingRegularizer(self.system, self.workspace_manager)
        success = regularizer.regularize()
        
        self.assertTrue(success)
        
        # Verify ProteinA (Ref) are on a sphere
        final_pos_a = np.array([m.com for m in mols_a])
        # Re-fit to find the exact center the regularizer settled on
        # (It should be close to ref_center)
        fit_center, fit_radius, fit_error = regularizer._fit_sphere(final_pos_a)
        
        self.assertLess(fit_error, 1e-6, "Reference species should have negligible radial error")
        np.testing.assert_allclose(fit_center, ref_center, atol=1.5)
        
        # Verify ProteinB are on a sphere centered at the predicted center
        final_pos_b = np.array([m.com for m in mols_b])
        radii_b = np.linalg.norm(final_pos_b - fit_center, axis=1)
        
        self.assertLess(np.std(radii_b), 1e-6, "Secondary species should be on a spherical shell")
        
        # Log message verification
        self.workspace_manager.logger.info.assert_any_call("Starting Sphere Projection Regularization")

    def test_regularize_empty_system(self):
        """Test regularization with empty system."""
        regularizer = RingRegularizer(self.system, self.workspace_manager)
        success = regularizer.regularize()
        self.assertFalse(success)

if __name__ == '__main__':
    unittest.main()

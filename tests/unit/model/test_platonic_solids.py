import unittest
import numpy as np
from ionerdss.model.PlatonicSolids import PlatonicSolidsModel
from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.components.instances import MoleculeInstance, InterfaceInstance
from ionerdss.model.components.reactions import ReactionRule
from ionerdss.model.components.system import System
import tempfile
import os
import shutil

class TestPlatonicSolidsModel(unittest.TestCase):
    def test_create_solid_cube(self):
        """Test creating a cube solid with standard components."""
        # Returns (System, List[ReactionRule])
        system, reactions = PlatonicSolidsModel.create_solid("cube", radius=10.0, sigma=1.0)
        
        # Verify System content
        self.assertIsInstance(system, System)
        self.assertEqual(len(system.molecule_types), 1)
        
        mol_type = system.molecule_types.get("cube")
        self.assertIsInstance(mol_type, MoleculeType)
        self.assertEqual(mol_type.name, "cube")
        self.assertEqual(mol_type.radius_nm, 10.0)
        
        # Verify InterfaceTypes in System
        # Cube has 4 sites -> 4 interfaces + 0 COM (COM is implicit origin)
        # Standard implementation logic: iterate 4 legs -> add 4 interfaces
        self.assertEqual(len(system.interface_types), 4)
        
        # Check an interface
        if1 = system.interface_types.get("cubecube1") # name format {this}{partner}{index}
        self.assertIsNotNone(if1)
        self.assertEqual(if1.interface_index, 1)
        self.assertTrue(isinstance(if1.local_coord, np.ndarray))

        # Verify Reactions
        # 4 sites combined with replacement: 4 self + 4*3/2 cross = 10 reactions
        # 4 sites combined with replacement: 4 self + 4*3/2 cross = 10 reactions
        self.assertEqual(len(reactions), 10)
        self.assertIsInstance(reactions[0], ReactionRule)
        self.assertEqual(reactions[0].geometry.sigma_nm, 1.0)

        # Verify Molecule Instances
        self.assertEqual(len(system.molecule_instances), 1)
        mol_inst = list(system.molecule_instances)[0]
        self.assertIsInstance(mol_inst, MoleculeInstance)
        self.assertEqual(mol_inst.molecule_type, mol_type)
        self.assertEqual(mol_inst.name, "cube_0")
        
        # Verify Interface Instances
        self.assertEqual(len(system.interface_instances), 4)
        for ii in system.interface_instances:
            self.assertIsInstance(ii, InterfaceInstance)
            self.assertEqual(ii.this_mol, mol_inst)
            # Check mapping
            self.assertIn(ii, mol_inst.interfaces_neighbors_map)
            self.assertIsNone(mol_inst.interfaces_neighbors_map[ii])

    def test_create_solid_dode(self):
        """Test creating a dodecahedron solid."""
        system, reactions = PlatonicSolidsModel.create_solid("dode", radius=10.0, sigma=1.0)
        
        # Dode has 5 sites
        self.assertEqual(len(system.interface_types), 5)
        # Reactions: 5 self + 5*4/2 = 15 total
        self.assertEqual(len(reactions), 15)

    def test_invalid_solid_type(self):
        """Test invalid solid type raises ValueError."""
        with self.assertRaises(ValueError):
            PlatonicSolidsModel.create_solid("invalid", radius=10.0, sigma=1.0)

    def test_missing_sigma_dode(self):
        """Test missing sigma raises ValueError."""
        with self.assertRaises(ValueError):
            PlatonicSolidsModel.create_solid("dode", radius=10.0, sigma=None)

    def test_reaction_attributes(self):
        """Verify generated reaction attributes using cube."""
        system, reactions = PlatonicSolidsModel.create_solid("cube", radius=10.0, sigma=2.0)
        reaction = reactions[0]
        
        self.assertIsInstance(reaction, ReactionRule)
        # Check geometry
        self.assertEqual(reaction.geometry.sigma_nm, 2.0)
        self.assertIsNotNone(reaction.geometry.theta1)
        self.assertTrue(len(reaction.geometry.norm1) == 3)
        self.assertTrue(isinstance(reaction.geometry.norm1, np.ndarray))
        
        # Check rate assignment (self vs cross)
        # First loop i=0, j=0 -> same site -> ka=120.0
        self.assertEqual(reaction.ka, 120.0)
        
        # Find a cross reaction (i != j)
        for r in reactions:
            if r.reactant_interfaces[0] != r.reactant_interfaces[1]:
                self.assertEqual(r.ka, 240.0)
                break

    def test_coordinates_validity(self):
        """Check coordinates are not all zero/None."""
        system, reactions = PlatonicSolidsModel.create_solid("cube", radius=10.0, sigma=1.0)
        
        for iface in system.interface_types:
            self.assertIsNotNone(iface.absolute_coord)
            self.assertIsNotNone(iface.local_coord)
            # Ensure they are numpy arrays
            self.assertIsNotNone(iface.local_coord)
            # Ensure they are numpy arrays
            self.assertTrue(isinstance(iface.absolute_coord, np.ndarray))

    def test_export_nerdss(self):
        """Test exporting to NERDSS format (integration test)."""
        system, reactions = PlatonicSolidsModel.create_solid("cube", radius=5.0, sigma=1.0)
        
        # Create temp directory
        with tempfile.TemporaryDirectory() as tmp_dir:
            PlatonicSolidsModel.export_nerdss(system, tmp_dir, reactions)
            
            # Check for NERDSS files
            # Expected structure: normal file organization handled by WorkspaceManager inside export_nerdss? 
            # PlatonicSolids.py: wm = WorkspaceManager(output_path, ...)
            # WorkspaceManager creates: structures/, outputs/, logs/... and exporter creates nerdss_files?
            # NerdssExporter: output_dir = workspace_manager.workspace_path / 'nerdss_files'
            
            nerdss_dir = os.path.join(tmp_dir, "nerdss_files")
            self.assertTrue(os.path.exists(nerdss_dir))
            
            # Check for .mol file
            mol_file = os.path.join(nerdss_dir, "cube.mol")
            self.assertTrue(os.path.exists(mol_file))
            
            # Check for parms.inp
            parms_file = os.path.join(nerdss_dir, "parms.inp")
            self.assertTrue(os.path.exists(parms_file))

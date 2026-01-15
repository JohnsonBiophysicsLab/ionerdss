"""Platonic Solids Model module for generating NERDSS molecule types and reactions.

This module uses a consolidated geometry generation approach to create geometric models 
of Platonic solids (cube, dodecahedron, etc.) using standard ionerdss components.
"""

from typing import List, Tuple, Dict
import numpy as np

# Standard component imports
from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.components.reactions import ReactionRule, ReactionGeometrySet
from ionerdss.model.components.system import System

# Import consolidated logic
from .platonic_solids.geometry import angle_cal
from .platonic_solids.solids import (
    CubeGenerator,
    DodecahedronGenerator,
    IcosahedronGenerator,
    OctahedronGenerator,
    TetrahedronGenerator,
    PlatonicSolidGenerator
)

class PlatonicSolidsModel:
    """A class for generating NERDSS molecule types and reactions for platonic solids."""
    
    # Registry of generator instances
    _GENERATORS: Dict[str, PlatonicSolidGenerator] = {
        "cube": CubeGenerator(),
        "dode": DodecahedronGenerator(),
        "icos": IcosahedronGenerator(),
        "octa": OctahedronGenerator(),
        "tetr": TetrahedronGenerator(),
    }

    @classmethod
    def create_solid(cls, solid_type: str, radius: float, sigma: float) -> Tuple[System, List[ReactionRule]]:
        """
        Create a System containing the Platonic solid definition and its reactions.

        Args:
            solid_type (str): The platonic solid type ["cube", "dode", "icos", "octa", "tetr"]
            radius (float): The radius of the circumscribed sphere (nm)
            sigma (float): Distance between two binding sites (nm)

        Returns:
            Tuple[System, List[ReactionRule]]: A tuple containing:
                - A System object populated with the MoleculeType and InterfaceTypes
                - A list of ReactionRule objects defining the binding interactions
        """
        if solid_type not in cls._GENERATORS:
            raise ValueError(f"Solid type must be one of {list(cls._GENERATORS.keys())}.")
        
        if sigma is None:
             raise ValueError(f"Sigma must be provided for solid type {solid_type}.")

        generator = cls._GENERATORS[solid_type]
        
        # 0. Initialize System
        system = System(workspace_path=".", pdb_id=f"{solid_type}_gen")

        # 1. Generate Coordinates (COM, Legs, Normal) for ALL faces
        # Returns List of [COM, leg1, leg2..., Normal]
        all_faces_coords = generator.generate_coordinates(radius, sigma)
        
        # Extract Representative Face Data (Face 0)
        # Structure: [COM, leg1, leg2, ..., Normal]
        # Normal is the LAST element. COM is the FIRST. legs are in between.
        face0_data = all_faces_coords[0]
        com = face0_data[0]
        normal = face0_data[-1]
        legs = face0_data[1:-1]
        
        # 2. Generate Angle Parameters
        # Use generator's angle indices to pick points from the generated faces
        idx1, idx2, idx3, idx4 = generator.angle_indices
        
        # Helper to extract point: (face_index, element_index)
        # element_index: 0=COM, 1=leg1...
        p1 = all_faces_coords[idx1[0]][idx1[1]]
        p2 = all_faces_coords[idx2[0]][idx2[1]]
        p3 = all_faces_coords[idx3[0]][idx3[1]]
        p4 = all_faces_coords[idx4[0]][idx4[1]]
        
        theta1, theta2, phi1, phi2, omega = angle_cal(p1, p2, p3, p4)

        # 3. Create MoleculeType
        mol_type = MoleculeType(name=solid_type, radius_nm=float(radius))
        mol_type.set_diffusion_constants_from_radius() # standard physics
        system.molecule_types.add(mol_type)

        # 4. Create InterfaceTypes and add to System
        # COM is treated as center of system/molecule? 
        # In this context, the entire solid is ONE particle in NERDSS if coarse-grained?
        # NO. Platonic solids simulation treats FACES as particles usually?
        # "dode_face_write" suggests we are simulating FACES as individual rigid bodies that assemble into the solid.
        # "create_Solid" implies creating a model OF THE SOLID.
        # But if we return one MoleculeType "cube" with 4 binding sites... that implies the Cube is ONE particle?
        # BUT `cube` MoleculeType has `radius` of the circumscribed sphere.
        # If the Cube is the particle, why do we need reaction angles between faces?
        # Standard NERDSS: "Patchy particles".
        # Yes, here 'cube' likely represents a single CUBE PARTICLE that binds to OTHER CUBE PARTICLES?
        # OR does 'cube' represent a SQUARE FACE that binds to form a cube? (Self-assembly of faces into solid).
        # "dode(lg1) + dode(lg1) <-> ..."
        # If it's self-assembly, then `MoleculeType` should be "Face".
        # But the name is `solid_type` ("cube").
        # If `num_sites=4` (legs of square face), then "cube" IS the face.
        # The naming is confusing: "cube" = "Square Face used to build a Cube".
        # "dode" = "Pentagon Face used to build a Dodecahedron".
        # This aligns with `num_sites` (4 for cube face, 5 for dode face).
        # So `com` calculated (Face COM) is the center of the particle.
        # And `legs` are the binding sites on the edges of the face.
        # `normal` is the orientation vector.
        
        # So for MoleculeType creation:
        # local_coord of interface = leg_coord - com.
        # Since `com` is the origin of the face particle essentially (or we define it so).
        # Actually `com` calculated by `generate_coordinates` is the position of the face in the assembled solid (relative to solid center).
        # But for the `MoleculeType` definition of a single free face, we want coordinates relative to the face center!
        # If `com` is [x,y,z], and `leg` is [lx, ly, lz].
        # Relative coord `leg - com` is correct for defining the reusable Face Template.
        
        interface_objects = []
        
        for i, leg_coord in enumerate(legs):
            index = i + 1
            local_coord = np.array(leg_coord) - np.array(com)
            
            # Absolute coord in the template definition is usually just the local coord if COM is origin.
            # `InterfaceType` constructor takes `absolute_coord` and `local_coord`. 
            # In `types.py`: local_coord = relative to COM. absolute_coord = global? 
            # But in a Type definition, global doesn't exist. Usually absolute=local for Type.
            
            interface = InterfaceType(
                this_mol_type_name=solid_type,
                partner_mol_type_name=solid_type, 
                interface_index=index,
                absolute_coord=local_coord, # For Type definition, absolute is usually same as local/relative to origin
                local_coord=local_coord,
                this_mol_type=mol_type,
                partner_mol_type=mol_type,
                energy=-1.0 
            )
            
            system.interface_types.add(interface)
            interface_objects.append(interface)

        # 5. Generate Reactions
        reactions = []
        
        for i in range(len(interface_objects)):
            for j in range(i, len(interface_objects)):
                site1 = interface_objects[i]
                site2 = interface_objects[j]
                
                geometry = ReactionGeometrySet(
                    theta1=theta1, theta2=theta2, 
                    phi1=phi1, phi2=phi2, 
                    omega=omega, 
                    sigma_nm=float(sigma),
                    norm1=normal, norm2=normal 
                )
                
                ka_val = 2.0 if i == j else 4.0
                
                reaction = ReactionRule(
                    expr="", 
                    reactant_interfaces=(site1, site2),
                    geometry=geometry,
                    ka=ka_val,
                    kb=1.0 
                )
                reactions.append(reaction)

        return system, reactions

    # Legacy alias
    create_Solid = create_solid

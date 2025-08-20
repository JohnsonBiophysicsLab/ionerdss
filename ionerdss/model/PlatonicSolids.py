"""Platonic Solids Model module for generating NERDSS molecule types and reactions, and corresponding files for specified platonic solid type.

This module defines the `PlantonicSolidsModel` class, inheriting from the `Model` class, which is used to generate NERDSS molecule types and reactions, and corresponding files for platonic solid.
"""

from .components import Model
from .components import Reaction
from .components import MoleculeType
from .components import MoleculeInterface
from ..utils.coords import Coords
from .platonic_solids.dode.dode_face_write import dode_face_write
from .platonic_solids.cube.cube_face_write import cube_face_write
from .platonic_solids.icos.icos_face_write import icos_face_write
from .platonic_solids.octa.octa_face_write import octa_face_write
from .platonic_solids.tetr.tetr_face_write import tetr_face_write
from dataclasses import dataclass, field
from typing import List, Dict, Tuple

class PlatonicSolid(Model):
    """"A class for generating NERDSS molecule types and reactions, and corresponding files for platonic solid.
    
    Attributes:
        pdb_file (str): The path to the PDB structure file.
        solid_type (str): The platonic solid type.
        binding_site_position (str): The binding site position.
    """
    name: str
    molecule_types: List[MoleculeType] = field(default_factory=list)
    reactions: List[Reaction] = field(default_factory=list)
  
    @classmethod
    def create_Solid(cls, solid_type: str ,radius:float,sigma: float=None,distances:float = None, mol_name=None) -> Model:
        """
        Parameters:
            cls:
            solid_type (str):  The platonic solid type within ["cube","dode","icos","octa","tetr"]
            radius (float): the radius of the circumscribed sphere around the platonic solid (nm) — that is, the distance from the center of the dodecahedron to any of its vertices.
            sigma (float): distance between two binding sites (nm)
            distances (float): not used
            mol_name: not used
        """
        types: list = ["cube","dode","icos","octa","tetr"]
        reactions_to_return = []   #this function returns all reactions generated in this module
        molecule_interfaces = []    #this returns the molecule information needed to generate a model class
        if solid_type not in types:
            raise ValueError(f"Solid type must be one of {types}.")
        if solid_type == 'dode':
            if sigma == None:
                raise ValueError(f"if Solid Type is {solid_type}. Sigma must be provided. Argument is currently sigma={sigma}")
            dode_reaction_parameters,dode_mol_information = dode_face_write(radius,sigma,create_Solid=True)
            dode_reactions:list = ['dode(lg1) + dode(lg1) <-> dode(lg1!1).dode(lg1!1)', 
                                    'dode(lg2) + dode(lg2) <-> dode(lg2!1).dode(lg2!1)', 
                                    'dode(lg3) + dode(lg3) <-> dode(lg3!1).dode(lg3!1)', 
                                    'dode(lg4) + dode(lg4) <-> dode(lg4!1).dode(lg4!1)', 
                                    'dode(lg5) + dode(lg5) <-> dode(lg5!1).dode(lg5!1)', 
                                    'dode(lg1) + dode(lg2) <-> dode(lg1!1).dode(lg2!1)', 
                                    'dode(lg1) + dode(lg3) <-> dode(lg1!1).dode(lg3!1)', 
                                    'dode(lg1) + dode(lg4) <-> dode(lg1!1).dode(lg4!1)', 
                                    'dode(lg1) + dode(lg5) <-> dode(lg1!1).dode(lg5!1)', 
                                    'dode(lg2) + dode(lg3) <-> dode(lg2!1).dode(lg3!1)', 
                                    'dode(lg2) + dode(lg4) <-> dode(lg2!1).dode(lg4!1)', 
                                    'dode(lg2) + dode(lg5) <-> dode(lg2!1).dode(lg5!1)', 
                                    'dode(lg3) + dode(lg4) <-> dode(lg3!1).dode(lg4!1)', 
                                    'dode(lg3) + dode(lg5) <-> dode(lg3!1).dode(lg5!1)', 
                                    'dode(lg4) + dode(lg5) <-> dode(lg4!1).dode(lg5!1)']
            
            norm = [float(dode_reaction_parameters['n'][0]),
                     float(dode_reaction_parameters['n'][1]),
                     float(dode_reaction_parameters['n'][2])]
            
            for i in dode_reactions:
                reactions_to_return.append(Reaction(
                    name=i,
                    binding_radius=float(sigma),
                    binding_angles=[dode_reaction_parameters['theta1'],
                                    dode_reaction_parameters['theta2'],
                                    dode_reaction_parameters['phi1'],
                                    dode_reaction_parameters['phi2'],
                                    dode_reaction_parameters['omega']],
                    norm1=norm,
                    norm2=norm,
                ))
            for i in dode_mol_information.keys():
                if i == "COM": continue
                #MoleculeInterface(name=vals[0], coord=Coords(x_coord, y_coord, z_coord)))
                x_coord = round(dode_mol_information[i][0],8)
                y_coord = round(dode_mol_information[i][1],8)
                z_coord = round(dode_mol_information[i][2],8)
                molecule_interfaces.append(MoleculeInterface(name=i,coord=Coords(x_coord,y_coord,z_coord)))

            molecule = [MoleculeType(name='dode',interfaces=molecule_interfaces)]
        if solid_type == 'cube':
            if sigma == None:
                raise ValueError(f"if Solid Type is {solid_type}. Sigma must be provided. Argument is currently sigma={sigma}")
            cube_reaction_parameters, cube_mol_information = cube_face_write(radius,sigma,create_Solid=True)
            
            cube_reactions:list = ['cube(lg1) + cube(lg1) <-> cube(lg1!1).cube(lg1!1)', 
                                    'cube(lg2) + cube(lg2) <-> cube(lg2!1).cube(lg2!1)', 
                                    'cube(lg3) + cube(lg3) <-> cube(lg3!1).cube(lg3!1)', 
                                    'cube(lg4) + cube(lg4) <-> cube(lg4!1).cube(lg4!1)', 
                                    'cube(lg1) + cube(lg2) <-> cube(lg1!1).cube(lg2!1)', 
                                    'cube(lg1) + cube(lg3) <-> cube(lg1!1).cube(lg3!1)', 
                                    'cube(lg1) + cube(lg4) <-> cube(lg1!1).cube(lg4!1)', 
                                    'cube(lg2) + cube(lg3) <-> cube(lg2!1).cube(lg3!1)', 
                                    'cube(lg2) + cube(lg4) <-> cube(lg2!1).cube(lg4!1)', 
                                    'cube(lg3) + cube(lg4) <-> cube(lg3!1).cube(lg4!1)', 
                                    ]
            
            norm = [float(cube_reaction_parameters['n'][0]),
                     float(cube_reaction_parameters['n'][1]),
                     float(cube_reaction_parameters['n'][2])]
            
            for i in cube_reactions:
                reactions_to_return.append(Reaction(
                    name=i,
                    binding_radius=float(sigma),
                    binding_angles=[cube_reaction_parameters['theta1'],
                                    cube_reaction_parameters['theta2'],
                                    cube_reaction_parameters['phi1'],
                                    cube_reaction_parameters['phi2'],
                                    cube_reaction_parameters['omega']],
                    norm1=norm,
                    norm2=norm,
                ))
            for i in cube_mol_information.keys():
                if i == "COM": continue
                
                x_coord = round(cube_mol_information[i][0],8)
                y_coord = round(cube_mol_information[i][1],8)
                z_coord = round(cube_mol_information[i][2],8)
                molecule_interfaces.append(MoleculeInterface(name=i,coord=Coords(x_coord,y_coord,z_coord)))

            molecule = [MoleculeType(name='cube',interfaces=molecule_interfaces)]
        if solid_type == 'icos':
            if sigma == None:
                raise ValueError(f"if Solid Type is {solid_type}. Sigma must be provided. Argument is currently sigma={sigma}")
            icos_reaction_parameters, icos_mol_information = icos_face_write(radius,sigma,create_Solid=True)
            
            icos_reactions:list = ['icos(lg1) + icos(lg1) <-> icos(lg1!1).icos(lg1!1)', 
                                    'icos(lg2) + icos(lg2) <-> icos(lg2!1).icos(lg2!1)', 
                                    'icos(lg3) + icos(lg3) <-> icos(lg3!1).icos(lg3!1)', 
                                    'icos(lg1) + icos(lg2) <-> icos(lg1!1).icos(lg2!1)', 
                                    'icos(lg1) + icos(lg3) <-> icos(lg1!1).icos(lg3!1)', 
                                    'icos(lg2) + icos(lg3) <-> icos(lg2!1).icos(lg3!1)', 
                                    ]
            
            norm = [float(icos_reaction_parameters['n'][0]),
                     float(icos_reaction_parameters['n'][1]),
                     float(icos_reaction_parameters['n'][2])]
            
            for i in icos_reactions:
                reactions_to_return.append(Reaction(
                    name=i,
                    binding_radius=float(sigma),
                    binding_angles=[icos_reaction_parameters['theta1'],
                                    icos_reaction_parameters['theta2'],
                                    icos_reaction_parameters['phi1'],
                                    icos_reaction_parameters['phi2'],
                                    icos_reaction_parameters['omega']],
                    norm1=norm,
                    norm2=norm,
                ))
            for i in icos_mol_information.keys():
                if i == "COM": continue
                
                x_coord = round(icos_mol_information[i][0],8)
                y_coord = round(icos_mol_information[i][1],8)
                z_coord = round(icos_mol_information[i][2],8)
                molecule_interfaces.append(MoleculeInterface(name=i,coord=Coords(x_coord,y_coord,z_coord)))

            molecule = [MoleculeType(name='icos',interfaces=molecule_interfaces)]
        if solid_type == 'octa':
            if sigma == None:
                raise ValueError(f"if Solid Type is {solid_type}. Sigma must be provided. Argument is currently sigma={sigma}")
            octa_reaction_parameters, octa_mol_information = octa_face_write(radius,sigma,create_Solid=True)
            
            octa_reactions:list = ['octa(lg1) + octa(lg1) <-> octa(lg1!1).octa(lg1!1)', 
                                    'octa(lg2) + octa(lg2) <-> octa(lg2!1).octa(lg2!1)', 
                                    'octa(lg3) + octa(lg3) <-> octa(lg3!1).octa(lg3!1)', 
                                    'octa(lg1) + octa(lg2) <-> octa(lg1!1).octa(lg2!1)', 
                                    'octa(lg1) + octa(lg3) <-> octa(lg1!1).octa(lg3!1)', 
                                    'octa(lg2) + octa(lg3) <-> octa(lg2!1).octa(lg3!1)', 
                                    ]
            
            norm = [float(octa_reaction_parameters['n'][0]),
                     float(octa_reaction_parameters['n'][1]),
                     float(octa_reaction_parameters['n'][2])]
            
            for i in octa_reactions:
                reactions_to_return.append(Reaction(
                    name=i,
                    binding_radius=float(sigma),
                    binding_angles=[octa_reaction_parameters['theta1'],
                                    octa_reaction_parameters['theta2'],
                                    octa_reaction_parameters['phi1'],
                                    octa_reaction_parameters['phi2'],
                                    octa_reaction_parameters['omega']],
                    norm1=norm,
                    norm2=norm,
                ))
            for i in octa_mol_information.keys():
                if i == "COM": continue
                
                x_coord = round(octa_mol_information[i][0],8)
                y_coord = round(octa_mol_information[i][1],8)
                z_coord = round(octa_mol_information[i][2],8)
                molecule_interfaces.append(MoleculeInterface(name=i,coord=Coords(x_coord,y_coord,z_coord)))

            molecule = [MoleculeType(name='octa',interfaces=molecule_interfaces)]
        if solid_type == 'tetr':
            if sigma == None:
                raise ValueError(f"if Solid Type is {solid_type}. Sigma must be provided. Argument is currently sigma={sigma}")
            tetr_reaction_parameters, tetr_mol_information = tetr_face_write(radius,sigma,create_Solid=True)
            
            tetr_reactions:list = ['tetr(lg1) + tetr(lg1) <-> tetr(lg1!1).tetr(lg1!1)', 
                                    'tetr(lg2) + tetr(lg2) <-> tetr(lg2!1).tetr(lg2!1)', 
                                    'tetr(lg3) + tetr(lg3) <-> tetr(lg3!1).tetr(lg3!1)', 
                                    'tetr(lg1) + tetr(lg2) <-> tetr(lg1!1).tetr(lg2!1)', 
                                    'tetr(lg1) + tetr(lg3) <-> tetr(lg1!1).tetr(lg3!1)', 
                                    'tetr(lg2) + tetr(lg3) <-> tetr(lg2!1).tetr(lg3!1)', 
                                    ]
            
            norm = [float(tetr_reaction_parameters['n'][0]),
                     float(tetr_reaction_parameters['n'][1]),
                     float(tetr_reaction_parameters['n'][2])]
            
            for i in tetr_reactions:
                reactions_to_return.append(Reaction(
                    name=i,
                    binding_radius=float(sigma),
                    binding_angles=[tetr_reaction_parameters['theta1'],
                                    tetr_reaction_parameters['theta2'],
                                    tetr_reaction_parameters['phi1'],
                                    tetr_reaction_parameters['phi2'],
                                    tetr_reaction_parameters['omega']],
                    norm1=norm,
                    norm2=norm,
                ))
            for i in tetr_mol_information.keys():
                if i == "COM": continue
                x_coord = round(tetr_mol_information[i][0],8)
                y_coord = round(tetr_mol_information[i][1],8)
                z_coord = round(tetr_mol_information[i][2],8)
                molecule_interfaces.append(MoleculeInterface(name=i,coord=Coords(x_coord,y_coord,z_coord)))

            molecule = [MoleculeType(name='tetr',interfaces=molecule_interfaces)]
        return cls(name=solid_type, molecule_types=molecule, reactions=reactions_to_return)
        

        

                    
"""
components.py

Author: sikaoguo@gmail.com, yying7@jh.edu

This module defines the core data structures and utilities for specifying and serializing
molecular models used in NERDSS (Non-Equilibrium Reaction-Diffusion Self-Assembly Simulator) simulations.

The primary components include:
- `MoleculeInterface`: Represents a binding or interaction site on a molecule.
- `MoleculeType`: Represents a coarse-grained molecular unit with defined interfaces and diffusion parameters.
- `ReactionType`: Represents a bimolecular binding interaction with spatial and angular constraints.
- `Model`: A container that aggregates molecules and reactions, and provides methods to save/load model definitions to/from JSON.

This file enables structured definition of molecular assemblies and reaction schemas,
with serialization support for use in automated pipeline generation and simulation input preparation.

Typical Usage:
    model = Model.load_model("my_model.json")
    model.save_model("backup_model.json")

Dependencies:
    - numpy
    - dataclasses
    - json
    - typing
    - local module: `coords.py` for Coords class
"""

import json

import numpy as np

from ionerdss.math.coords import Coords

#=============================
# Molecules
#=============================

class MoleculeInterface:
    """Represents an interface of a molecule type.
    
    Attributes:
        name (str): Name of the interface.
        coord (Tuple[float, float, float]): Coordinates of the interface.
    """
    def __init__(self, name, coord):
        self.name = name
        self.coord = coord

    def __repr__(self):
        return f"<MoleculeInterface {self.name} @ {self.coord}>"

class MoleculeType:
    """
    Represents a molecule type in the model.

    Attributes
    ----------
    name : str
        Name of the molecule type.
    interfaces : list of MoleculeInterface
        List of interfaces associated with the molecule.
    diffusion_translation : float
        Translational diffusion const (default 0.0).
    diffusion_rotation : float
        Rotational diffusion const (default 0.0).
    """
    def __init__(self, name, interfaces=None, translational_diffusion_constant=0.0, rotational_diffusion_constant=0.0):
        self.name = name
        self.interfaces = interfaces if interfaces is not None else []
        self.translational_diffusion_constant = translational_diffusion_constant
        self.rotational_diffusion_constant = rotational_diffusion_constant

    def __repr__(self):
        return f"<MoleculeType {self.name} with {len(self.interfaces)} interfaces>"

class MoleculeTemplate:
    """
    Represents a molecule type in NERDSS, including the molecule's center of mass (COM) 
    and a list of binding interfaces.

    Attributes:
        name (str): Identifier for the molecule type.
        interface_template_list (list): A list of BindingInterfaceTemplate objects that 
            describe the molecule’s binding sites.
        normal_point (list): Default normal vector direction.
    """
    def __init__(self, name: str = '',
                 com=None,
                 interfaces=None,
                 translational_diffusion_constant = 1.0,
                 rotational_diffusion_constant = 1.0,
                 radius: float = 1.0):
        """
        Initializes a MoleculeTemplate.

        Parameters
        ----------
        name : str, optional
            Name or identifier of the molecule type. Default is empty string.
        com : array-like of shape (3,), optional
            Center of mass coordinates. Defaults to [0.0, 0.0, 0.0].
        radius : float, optional
            Effective molecular radius. Default is 1.0.
        """
        self.name = name
        self.interface_template_list = []
        self.com = np.asarray(com if com is not None else [0.0, 0.0, 0.0], dtype=float)
        self.interfaces = interfaces if interfaces is not None else []
        self.normal_point = [0.0, 0.0, 1.0]
        self.translational_diffusion_constant = translational_diffusion_constant
        self.rotational_diffusion_constant = rotational_diffusion_constant
        self.radius = float(radius)

    def add_interface(self, interface):
        """
        Add an InterfaceTemplate to this molecule.

        Parameters
        ----------
        interface : InterfaceTemplate
            The interface object to add.
        """
        self.interfaces.append(interface)

    def __str__(self):
        interfaces = "\n  ".join(str(it) for it in self.interface_template_list)
        return f"Molecule Template: {self.name}\n  Interfaces:\n  {interfaces}"

    def __repr__(self):
        interfaces = "\n  ".join(str(it) for it in self.interface_template_list)
        return f"Molecule Template: {self.name}\n  Interfaces:\n  {interfaces}"

    def __eq__(self, other):
        if not isinstance(other, MoleculeTemplate):
            return False
        return self.name == other.name

class CoarseGrainedMolecule:
    """
    Represents a coarse-grained molecule in NERDSS, potentially derived from a PDB chain.

    Attributes:
        name (str): Identifier of the molecule.
        my_template (MoleculeTemplate): Reference to the associated molecule template.
        coord (Coords): Center-of-mass coordinates.
        interface_list (list): List of binding interfaces.
        normal_point (list): Normal vector direction.
    """

    def __init__(
            self,
            name: str,
            template=None,
            coord=None,
            interfaces=None,
            normal_point=None,
            translational_diffusion_constant=1.0,
            rotational_diffusion_constant=1.0,
            radius=1.0
        ):
        """
        Initializes a CoarseGrainedMolecule.

        Parameters
        ----------
        name : str
            Name or identifier of the molecule.
        template : Any, optional
            Reference structure or coordinate template.
        coord : Any, optional
            Center of mass coordinates or 3D origin of the molecule.
        interface_list : list, optional
            List of interaction site coordinates.
        normal_point : Any, optional
            Optional orientation vector or normal.
        diffusion_translation : list or np.ndarray, optional
            3D translational diffusion coefficients.
        diffusion_rotation : list or np.ndarray, optional
            3D rotational diffusion coefficients.
        radius : float, optional
            Effective molecular radius.
        """
        self.name = name
        self.template = template
        self.coord = coord
        self.interfaces = interfaces if interfaces is not None else []
        self.normal_point = normal_point
        self.translational_diffusion_constant = translational_diffusion_constant
        self.rotational_diffusion_constant = rotational_diffusion_constant
        self.radius = radius

    def __str__(self):
        interfaces = "\n  ".join(str(interface) for interface in self.interfaces)
        return (f"CoarseGrainedMolecule: {self.name}\n"
                f"  Template: {self.template}\n"
                f"  Coordinates: {self.coord}\n"
                f"  Interfaces:\n  {interfaces}")

    def __repr__(self):
        # Similar to __str__ but more formal for debugging
        return self.name

    def __eq__(self, other):
        if not isinstance(other, CoarseGrainedMolecule):
            return False
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

class Reaction:
    """
    Represents a reaction in the model.

    Attributes
    ----------
    name : str
        Reaction expression.
    binding_radius : float
        Binding radius of the reaction.
    binding_angles : tuple of float
        Binding angles (theta1, theta2, phi1, phi2, chi).
    norm1 : tuple of float
        First normal vector.
    norm2 : tuple of float
        Second normal vector.
    ka : float
        Forward rate constant.
    kb : float
        Reverse rate constant.
    """
    def __init__(self, name, binding_radius, binding_angles,
                 norm1, norm2, ka=0.0, kb=0.0):
        self.name = name
        self.binding_radius = binding_radius
        self.binding_angles = binding_angles
        self.norm1 = norm1
        self.norm2 = norm2
        self.ka = ka
        self.kb = kb

    def __repr__(self):
        return f"<Reaction {self.name}, ka={self.ka}, kb={self.kb}>"

class Model:
    """
    Parent class for all models to generate input files for NERDSS simulations.

    Attributes
    ----------
    name : str
        Name of the model.
    molecule_types : list of MoleculeType
        List of molecule types in the model.
    reactions : list of Reaction
        List of reactions in the model.
    """
    def __init__(self, name, molecule_types=None, reactions=None):
        self.name = name
        self.molecule_types = molecule_types if molecule_types is not None else []
        self.reactions = reactions if reactions is not None else []

    def __repr__(self):
        return f"<Model {self.name}: {len(self.molecule_types)} molecules, {len(self.reactions)} reactions>"

    def save_model(self, file_path: str) -> None:
        """Saves the model to a specified JSON file.
        
        Args:
            file_path (str): Path to the file where the model should be saved.
        """
        data = {
            "name": self.name,
            "molecule_types": [
                {
                    "name": mol.name,
                    "interfaces": [
                        {"name": iface.name, "coord": iface.coord} for iface in mol.interfaces
                    ],
                    "diffusion_translation": mol.diffusion_translation,
                    "diffusion_rotation": mol.diffusion_rotation,
                }
                for mol in self.molecule_types
            ],
            "reactions": [
                {
                    "name": rxn.name,
                    "binding_radius": rxn.binding_radius,
                    "binding_angles": rxn.binding_angles,
                    "norm1": rxn.norm1,
                    "norm2": rxn.norm2,
                    "ka": rxn.ka,
                    "kb": rxn.kb,
                }
                for rxn in self.reactions
            ],
        }
        with open(file_path, "w", encoding="utf-8") as file:
            json.dump(data, file, indent=4, cls=CustomJSONEncoder)  # Use custom encoder

    def load_model(self, cls, file_path: str) -> "Model":
        """Loads a model from a specified JSON file.
        
        Args:
            file_path (str): Path to the JSON file containing the model data.
        
        Returns:
            Model: An instance of the Model class with the loaded data.
        """
        with open(file_path, "r", encoding="utf-8") as file:
            data = json.load(file)

        molecule_types = [
            MoleculeType(
                name=mol["name"],
                interfaces=[MoleculeInterface(name=iface["name"], coord=Coords(**iface["coord"])) for iface in mol["interfaces"]],
                translational_diffusion_constant=mol["diffusion_translation"],
                rotational_diffusion_constant=mol["diffusion_rotation"],
            )
            for mol in data["molecule_types"]
        ]

        reactions = [
            Reaction(
                name=rxn["name"],
                binding_radius=rxn["binding_radius"],
                binding_angles=tuple(rxn["binding_angles"]),
                norm1=tuple(rxn["norm1"]),
                norm2=tuple(rxn["norm2"]),
                ka=rxn["ka"],
                kb=rxn["kb"],
            )
            for rxn in data["reactions"]
        ]

        return cls(name=data["name"], molecule_types=molecule_types, reactions=reactions)

class CustomJSONEncoder(json.JSONEncoder):
    """Custom JSON encoder to handle Coords serialization and NumPy types."""
    def default(self, o):
        if isinstance(o, Coords):
            return {"x": o.x, "y": o.y, "z": o.z}  # Convert Coords to dict
        elif isinstance(o, np.float32):  # Convert numpy float32 to standard float
            return float(o)
        elif isinstance(o, np.ndarray):  # Convert numpy array to list
            return o.tolist()
        return super().default(o)

class BindingInterfaceTemplate:
    """
    Represents a binding interface template between molecules.

    Attributes:
        name (str): Identifier of the interface template.
        coord (Coords): Relative coordinates of the interface.
        my_residues (list): Residues forming this interface.
        required_free_list (list): Other interface templates that must remain unbound 
            for this interface to bind.
        signature (dict): Stores interface geometry information.
    """

    def __init__(self, name: str = '', residues = None,
                 coords = None, energy = None):
        """
        Initializes a BindingInterfaceTemplate.

        Args:
            name (str): Identifier for the interface template.
        """
        self.name = name
        if coords is None:
            coords = [0.0, 0.0, 0.0]
        self.coords = np.asarray(coords, dtype=float)
        self.residues = residues if residues is not None else []
        # The list of interface templates that need
        # to be free to bind to this interface template
        self.required_free_list = []
        self.signature = {}
        self.energy = energy

    def __repr__(self):
        return f"<InterfaceTemplate {self.name} @ {self.coords.tolist()}>"

    def __str__(self):
        residues = ", ".join(self.residues)
        required_free = ", ".join(self.required_free_list)
        return (f"Interface Template: {self.name}\n"
                f"  Coordinates: {self.coords}\n"
                f"  Residues: {residues}\n"
                f"  Required Free: {required_free}")

    def __eq__(self, other):
        if not isinstance(other, BindingInterfaceTemplate):
            return False
        return self.name == other.name

class BindingInterface:
    """
    Represents a binding interface between molecules.

    Attributes:
        name (str): Identifier of the binding interface.
        coord (Coords): Position of the interface.
        my_template (BindingInterfaceTemplate): Reference to the associated interface template.
        my_residues (list): Residues included in the interface.
        signature (dict): Stores interface geometry information.
    """

    def __init__(self, name: str):
        """
        Initializes a BindingInterface.

        Args:
            name (str): Identifier for the binding interface.
        """
        self.name = name
        self.coord = None
        self.my_template = None
        self.my_residues = []
        self.signature = {}
        self.energy = None

    def __str__(self):
        return (f"BindingInterface: {self.name}\n"
                f"  Template: {self.my_template}\n"
                f"  Coordinates: {self.coord}\n"
                f"  Residue Count: {len(self.my_residues)}\n"
                f"  Residues: {self.my_residues}")

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, BindingInterface):
            return False
        return self.my_template == other.my_template

class ReactionTemplate:
    """
    Defines a reaction template between two MoleculeTemplates.

    Attributes:
        expression (str): Textual representation of the reaction.
        reactants (list): List of reactant molecule/interface templates.
        products (list): List of product molecule/interface templates.
        binding_angles (tuple): Tuple describing binding angles (theta1, theta2, phi1, phi2, omega).
        binding_radius (float): Distance between binding interfaces.
        norm1 (list): Normal vector of the first reactant.
        norm2 (list): Normal vector of the second reactant.
    """
    def __init__(self,
                 reactants = None,
                 products = None,
                 binding_angles = None,
                 binding_radius = 3.5,
                 norm1 = None,
                 norm2 = None,
                 ka = 1.0,
                 kb = 1.0):
        """
        Initializes a ReactionTemplate with default values.
        """
        self.expression = None
        self.reactants = reactants if reactants is not None else []
        self.products = products if products is not None else []
        self.binding_angles = binding_angles if binding_angles is not None else []
        self.binding_radius = binding_radius
        self.norm1 = norm1
        self.norm2 = norm2
        self.kd = None
        self.ka = ka
        self.kb = kb
        self.energy = None

    def __str__(self):
        return (f"Reaction Template: {self.expression}\n"
                f"  Reactants: {self.reactants}\n"
                f"  Products: {self.products}\n"
                f"  Binding Angles: {self.binding_angles}\n"
                f"  Binding Radius: {self.binding_radius / 10:.6f} nm\n"
                f"  norm1: {self.norm1}\n"
                f"  norm2: {self.norm2}")

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, ReactionTemplate):
            return False
        return self.expression == other.expression

"""
molecule.py

Defines the core data structures for representing coarse-grained molecules and their
interfaces within the NERDSS model generation pipeline.

Each `MoleculeTemplate` represents a protein chain (or homologous group) with:
- A center-of-mass (COM)
- A radius (used in steric clash detection or visualization)
- A list of `InterfaceTemplate` objects representing contact points

Each `InterfaceTemplate` includes:
- A 3D coordinate (relative or absolute)
- A list of residues contributing to the interface
- An optional energy (e.g., from contact scoring or docking)

These classes are instantiated after coarse-graining and regularization, and then passed
to modules like `reaction.py` for pairing and `export.py` for simulation generation.

Classes
-------
- MoleculeTemplate
- InterfaceTemplate
"""

import numpy as np

class InterfaceTemplate:
    """
    Represents a single protein interface (binding site) for coarse-grained modeling.

    Attributes
    ----------
    name : str
        Unique identifier (e.g., "I0_1").
    coords : ndarray of shape (3,)
        Interface center in 3D space (absolute or relative to COM).
    residues : list of int
        Residue indices that define this interface.
    energy : float
        Estimated binding energy or contact strength (optional).
    """
    def __init__(self, name, coords, residues=None, energy=None):
        self.name = name
        self.coords = np.asarray(coords, dtype=float)
        self.residues = residues if residues is not None else []
        self.energy = energy

    def __repr__(self):
        return f"<InterfaceTemplate {self.name} @ {self.coords.tolist()}>"


class MoleculeTemplate:
    """
    Represents a coarse-grained protein or subunit in the NERDSS model.

    Attributes
    ----------
    name : str
        Descriptive name (e.g., "MOL_0").
    com : ndarray of shape (3,)
        Center of mass of the molecule.
    radius : float
        Effective radius of the molecule.
    interfaces : list of InterfaceTemplate
        All contact sites exposed by this molecule.
    """
    def __init__(self, name, com, radius):
        self.name = name
        self.com = np.asarray(com, dtype=float)
        self.radius = float(radius)
        self.interfaces = []

    def add_interface(self, interface):
        """
        Add an InterfaceTemplate to this molecule.

        Parameters
        ----------
        interface : InterfaceTemplate
            The interface object to add.
        """
        self.interfaces.append(interface)

    def __repr__(self):
        return f"<MoleculeTemplate {self.name} with {len(self.interfaces)} interfaces>"

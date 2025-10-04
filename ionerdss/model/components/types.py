"""
ionerdss.model.components.types

Molecular type definitions and templates for simulation components.

This module provides the template classes that define the types of molecules
and interfaces used in molecular simulations. These classes serve as blueprints
or templates from which concrete instances are created during simulation runtime.

Classes:
    InterfaceType: Template definition for molecular binding interfaces, including
        spatial coordinates, energy parameters, and binding partner relationships.
    MoleculeType: Template definition for molecule types, including physical
        properties like size and diffusion constants, plus available interfaces.

Architecture Overview:
    The type classes follow a template pattern where:
    - InterfaceType defines the blueprint for binding sites
    - MoleculeType defines the blueprint for entire molecules
    - Each type can be instantiated multiple times during simulation
    - Types maintain references to their binding partners and relationships

Key Features:
    - Serialization support via to_dict() and from_dict() methods
    - Automatic diffusion constant calculation from molecular radius
    - Interface naming conventions: "{mol}_{partner}_{index}"
    - Coordinate handling for both absolute and relative positions

Dependencies:
    - numpy: For coordinate arrays and mathematical operations
    - ionerdss.utils.diffusion_constant: For computing physical properties

Example:
    ```python
    # Create a molecule type
    protein_A = MoleculeType(
        name="ProteinA",
        radius_nm=2.5
    )
    protein_A.set_diffusion_constants_from_radius()
    
    # Create an interface type
    binding_site = InterfaceType(
        this_mol_type_name="ProteinA",
        partner_mol_type_name="ProteinB",
        interface_index=1,
        absolute_coord=np.array([0.0, 0.0, 0.0]),
        local_coord=np.array([1.0, 0.0, 0.0]),
        energy=-5.0
    )
    ```

See Also:
    ionerdss.model.components.instances: Runtime instances created from these types
"""

from __future__ import annotations  # fix forward reference problem
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from ionerdss.utils.diffusion_constant import compute_diffusion_constants_nm_us

# ---------------------------------------------------- #


@dataclass
class InterfaceType:
    """Represents a type of molecular interfaces binding site between two molecules.

    An InterfaceType defines a specific binding interface type on a molecule that 
    can interact with a partner molecule. Each interface has spatial coordinates,
    energy information, and maintains references to its binding partner and
    parent molecule types.

    The interface identifier follows the format: "{this_mol}_{partner_mol}_{index}"
    where index allows multiple interfaces of the same type on one molecule.

    Attributes:
        this_mol_type_name: Name of the molecule containing this interface (e.g., "A").
        partner_mol_type_name: Name of the target binding partner molecule (e.g., "B").
        interface_index: Numeric index for distinguishing multiple interfaces 
            of the same type on one molecule.
        absolute_coord: Absolute Cartesian coordinates of the interface.
        local_coord: Coordinates relative to center of mass (template space).
        partner_interface_type: The corresponding InterfaceType on the partner
            molecule that this interface binds to.
        this_mol_type: Reference to the MoleculeType containing this interface.
        partner_mol_type: Reference to the partner's MoleculeType.
        required_free: List of interface names that must remain unbound to
            avoid steric clashes.
        energy: Binding energy of this interface interaction. Defaults to -1.0.
        signature: Dictionary containing interface-specific properties and
            metadata.
    """
    # Required fields first
    this_mol_type_name: str          # e.g., "A"
    partner_mol_type_name: str       # e.g., "B"
    interface_index: int        # e.g., 1 <- note, separate index for different `this_mol`
    absolute_coord: np.ndarray        # absolute cartesian coord
    local_coord: np.ndarray           # relative to COM (template space)

    # Optional fields with defaults
    partner_interface_type: Optional['InterfaceType'] = None  # e.g. B_A_1
    this_mol_type: Optional['MoleculeType'] = None           # e.g. A
    partner_mol_type: Optional['MoleculeType'] = None        # e.g. B
    # interfaces with steric clash
    required_free: List[str] = field(default_factory=list)
    energy: Optional[float] = -1.0
    signature: Dict = field(default_factory=dict)

    def get_name(self) -> str:
        """Return the formatted interface identifier string.

        Constructs the interface name using the format:
        "{this_mol_name}_{partner_mol_name}_{interface_index}"

        Returns:
            The interface identifier string (e.g., "A_B_1").
        """
        return self.this_mol_type_name + "_" +\
            self.partner_mol_type_name + "_" +\
            str(self.interface_index)

    def set_name(self, new_name: str) -> None:
        """Parse and set interface identifiers from a formatted name string.

        Parses a name string in the format "{mol}_{partner}_{index}" and
        updates the corresponding instance attributes.

        Args:
            new_name: Interface name string to parse (e.g., "A_B_1").

        Raises:
            IndexError: If the name string doesn't contain exactly 3 
                underscore-separated components.
            ValueError: If the interface index cannot be converted to integer.
        """
        substrings = new_name.split("_")
        self.this_mol_type_name = substrings[0]
        self.partner_mol_type_name = substrings[1]
        self.interface_index = int(substrings[2])

    def to_dict(self) -> dict:
        """Convert the InterfaceType to a dictionary representation.

        Serializes the interface data to a dictionary format suitable for
        JSON export or storage. Coordinate arrays are converted to lists.

        Returns:
            Dictionary containing the interface data with keys: 'name',
            'partner_interface_type', 'this_mol_type', 'absolute_coord', 
            'local_coord', 'required_free', 'energy', and 'signature'.
        """
        return {
            "name": self.get_name(),
            "partner_interface_type": self.partner_interface_type.get_name() if self.partner_interface_type else None,
            "this_mol_type": self.this_mol_type.name if self.this_mol_type else None,
            "absolute_coord": self.absolute_coord.tolist(),
            "local_coord": self.local_coord.tolist(),
            "required_free": list(self.required_free),
            "energy": self.energy,
            "signature": self.signature,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "InterfaceType":
        """Create an InterfaceType instance from a dictionary.

        Deserializes interface data from a dictionary format, typically
        loaded from JSON or other storage formats. Missing optional fields
        are filled with default values.

        Args:
            d: Dictionary containing interface data. Required keys are
                'this_mol_type_name', 'partner_mol_type_name', 'interface_index',
                'absolute_coord', and 'local_coord'. Optional keys include
                'partner_interface_type', 'this_mol_type', 'partner_mol_type',
                'required_free', 'energy', and 'signature'.

        Returns:
            New InterfaceType instance populated with the dictionary data.

        Raises:
            KeyError: If required keys are missing from the dictionary.
        """
        return cls(
            this_mol_type_name=d["this_mol_type_name"],
            partner_mol_type_name=d["partner_mol_type_name"],
            interface_index=d["interface_index"],
            absolute_coord=np.array(d["absolute_coord"]),
            local_coord=np.array(d["local_coord"]),
            partner_interface_type=d.get("partner_interface_type", None),
            this_mol_type=d.get("this_mol_type", None),
            partner_mol_type=d.get("partner_mol_type", None),
            required_free=list(d.get("required_free", [])),
            energy=d.get("energy", -1.0),
            signature=d.get("signature", {}),
        )


# ---------------------------------------------------- #


@dataclass
class MoleculeType:
    """Represents a type of molecule with its physical properties and binding interfaces.

    A MoleculeType defines the template for a specific kind of molecule, including
    its physical parameters (size, diffusion constants) and all possible binding
    interfaces it can form with other molecules. Multiple molecule instances can
    be created from the same MoleculeType.

    Attributes:
        name: Unique identifier for this molecule type (e.g., "ProteinA").
        interfaces_neighbors_map: Dictionary mapping interface names to neighbor
            molecule types that can bind through those interfaces.
        radius_nm: Molecular radius in nanometers. Defaults to 0.0.
        D_t_nm2_us: Translational diffusion constant in nm²/μs. Defaults to 0.0.
        D_r_rad2_us: Rotational diffusion constant in rad²/μs. Defaults to 0.0.
    """
    # Required fields
    name: str

    # Optional fields with defaults
    interfaces_neighbors_map: Dict = field(default_factory=dict)
    radius_nm: Optional[float] = 0.0
    D_t_nm2_us: Optional[float] = 0.0
    D_r_rad2_us: Optional[float] = 0.0

    def set_diffusion_constants_from_radius(self) -> None:
        """Calculate and set diffusion constants based on the molecular radius.

        Uses the Stokes-Einstein relation to compute translational and rotational
        diffusion constants from the molecular radius. Updates the D_t_nm2_us and
        D_r_rad2_us attributes in place.

        Note:
            Requires that radius_nm is set to a positive value before calling.
            The computation assumes spherical molecules in aqueous solution.
        """
        self.D_t_nm2_us, self.D_r_rad2_us =\
            compute_diffusion_constants_nm_us(radius_nm=self.radius_nm)

    def to_dict(self) -> dict:
        """Convert the MoleculeType to a dictionary representation.

        Serializes the molecule type data to a dictionary format suitable for
        JSON export or storage. Interface types are referenced by name to
        avoid circular references.

        Returns:
            Dictionary containing the molecule data with keys: 'name', 'radius',
            'diffusion_translation', 'diffusion_rotation', and 'interfaces_neighbors_map'.
        """
        return {
            "name": self.name,
            "radius": float(self.radius_nm),
            "diffusion_translation": float(self.D_t_nm2_us),
            "diffusion_rotation": float(self.D_r_rad2_us),
            "interfaces_neighbors_map": dict(self.interfaces_neighbors_map),
        }

    @classmethod
    def from_dict(cls, d: dict) -> "MoleculeType":
        """Create a MoleculeType instance from a dictionary.

        Deserializes molecule type data from a dictionary format, typically
        loaded from JSON or other storage formats. Missing optional fields
        are filled with default values.

        Args:
            d: Dictionary containing molecule data. Required key is 'name'.
                Optional keys include 'radius', 'diffusion_translation',
                'diffusion_rotation', and 'interfaces_neighbors_map'.

        Returns:
            New MoleculeType instance populated with the dictionary data.

        Raises:
            KeyError: If the required 'name' key is missing from the dictionary.
        """
        return cls(
            name=d["name"],
            radius_nm=float(d.get("radius", 0.0)),
            D_t_nm2_us=float(d.get("diffusion_translation", 0.0)),
            D_r_rad2_us=float(d.get("diffusion_rotation", 0.0)),
            interfaces_neighbors_map=d.get("interfaces_neighbors_map", {}),
        )

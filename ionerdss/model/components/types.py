"""
ionerdss.model.components.types
================================

Molecular type definitions and templates for simulation components.

This module provides template classes (``InterfaceType`` and ``MoleculeType``) that define
the blueprints for molecules and binding interfaces. These templates are instantiated
multiple times during simulation to create concrete molecule and interface instances.

**Unit Conventions:**
    - **Coordinates**: Nanometers (nm) - all spatial positions
    - **Radius**: Nanometers (nm) - molecular size
    - **Diffusion (translational)**: nm²/μs - translational mobility  
    - **Diffusion (rotational)**: rad²/μs - rotational mobility
    - **Energy**: kT (dimensionless) - binding free energy

**Coordinate Systems:**
    - **Global frame**: Laboratory/simulation box coordinates (nm)
    - **Local frame**: Molecule-centered coordinates relative to COM (nm)
    - **Reference vectors (ref1, ref2)**: Define molecule local orientation
        - ``ref1_local``: Primary reference axis (default: X-axis = [1, 0, 0])
        - ``ref2_local``: Secondary reference axis (default: Z-axis = [0, 0, 1])
        - Together define right-handed coordinate system for each molecule

**Architecture:**
    Templates follow a blueprint pattern:
    
    - ``MoleculeType``: Template for molecule species (e.g., "ProteinA")
        - Defines physical properties (size, diffusion)
        - Lists available binding interfaces
        - Stores local reference frame (ref1, ref2)
        
    - ``InterfaceType``: Template for binding site types (e.g., "ProteinA_ProteinB_1")
        - Defines binding geometry (position, orientation)
        - Specifies binding partner
        - Links to parent molecule type

**Key Classes:**
    - ``InterfaceType``: Binding site template with coordinates and energy
    - ``MoleculeType``: Molecule template with physical properties

**Reference Frame System:**
    Each molecule has a local coordinate system defined by two orthogonal reference vectors:
    
    .. code-block:: python
    
        # Default reference frame (right-handed)
        ref1_local = [1.0, 0.0, 0.0]  # X-axis (primary reference)
        ref2_local = [0.0, 0.0, 1.0]  # Z-axis (secondary reference)
        # Implied: ref3 = cross(ref1, ref2) = [0.0, 1.0, 0.0] (Y-axis)
    
    **Purpose**: Reference vectors allow NERDSS to:
    
    - Calculate binding angles relative to molecule orientation
    - Apply rotations when molecules bind
    - Maintain consistent geometry across multiple instances

**Examples:**
    Creating a molecule type with interface:
    
    >>> from ionerdss.model.components.types import MoleculeType, InterfaceType
    >>> import numpy as np
    >>> 
    >>> # Create molecule type
    >>> protein_a = MoleculeType(
    ...     name="ProteinA",
    ...     radius_nm=2.5  # nm
    ... )
    >>> protein_a.set_diffusion_constants_from_radius()
    >>> 
    >>> # Create binding interface
    >>> binding_site = InterfaceType(
    ...     this_mol_type_name="ProteinA",
    ...     partner_mol_type_name="ProteinB",
    ...     interface_index=1,
    ...     absolute_coord=np.array([10.0, 5.0, 3.0]),  # nm (global)
    ...     local_coord=np.array([1.5, 0.0, 0.0]),      # nm (relative to COM)
    ...     energy=-5.0  # kT
    ... )

See Also:
    - ``ionerdss.model.components.instances``: Runtime instances from these templates
    - ``ionerdss.utils.diffusion_constant``: Diffusion constant calculations

.. note::
    **Critical**: All coordinates must be in nanometers (nm) when creating types.
    The PDB parser converts Ångström→nm during data extraction.
"""

from __future__ import annotations # fix forward reference problem
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
    required_free: List[str] = field(default_factory=list) # interfaces with steric clash
    signature: Dict = field(default_factory=dict)
    
    partner_interface_type: Optional['InterfaceType'] = None # e.g. B_A_1
    this_mol_type: Optional['MoleculeType'] = None           # e.g. A
    partner_mol_type: Optional['MoleculeType'] = None        # e.g. B
    energy: Optional[float] = -1.0
    tag: Optional[str] = None   # None | 'f' | 'b'

    def get_name(self) -> str:
        """Return the formatted interface identifier string.
        
        Constructs the interface name using the format:
        "{this_mol_name}{partner_mol_name}{interface_index}"
        WITHOUT underscores to match the parser regex pattern.
        
        Returns:
            The interface identifier string (e.g., "AB1" or "AA1f").
        """
        core = self.this_mol_type_name + self.partner_mol_type_name + str(self.interface_index)
        return f"{core}{self.tag}" if self.tag else core

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
        
        # check if last character is a alphabet (assuming f or b)
        # A_A_1f -> substrings[2] gives 1f
        last_char = substrings[2][-1] # f

        if last_char.isalpha():
            # Split the string into the part before the last character and the last character
            string_without_last = substrings[2][:-1]
            self.interface_index = int(string_without_last)
            self.tag = last_char
        else:
            self.interface_index = int(substrings[2])
            self.tag = None

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
    radius_nm: Optional[float] = 0.0  # UNITS: Nanometers (nm)
    D_t_nm2_us: Optional[float] = 0.0  # UNITS: Translational diffusion in nm²/μs
    D_r_rad2_us: Optional[float] = 0.0  # UNITS: Rotational diffusion in rad²/μs
    
    # Optional ref coords (two vectors to robustly define molecule orientation)
    # COORDINATE SYSTEM: Local reference frame for this molecule type
    # ref1_local: Primary reference axis (default X-axis)
    # ref2_local: Secondary reference axis (default Z-axis)  
    # These define the molecule's intrinsic orientation for angle calculations
    ref1_local: Optional[np.ndarray] = field(default_factory=lambda: np.array([1.0, 0.0, 0.0]))  # X-axis (primary)
    ref2_local: Optional[np.ndarray] = field(default_factory=lambda: np.array([0.0, 0.0, 1.0]))  # Z-axis (secondary)

    
    def set_diffusion_constants_from_radius(self) -> None:
        """Calculate and set diffusion constants based on the molecular radius.
        
        Uses the Stokes-Einstein relation to compute translational and rotational
        diffusion constants from the molecular radius. Updates D_t_nm2_us and
        D_r_rad2_us attributes in place.
        
        **Physical Model:**
            Treats molecule as a hard sphere in aqueous solution:
            
            - D_t = kT / (6πηr)  - Translational diffusion
            - D_r = kT / (8πηr³) - Rotational diffusion
            
            where:
            - k: Boltzmann constant
            - T: Temperature (typically 293 K)
            - η: Viscosity of water
            - r: Molecular radius (nm)
        
        Note:
            Requires that radius_nm is set to a positive value before calling.
            The computation assumes spherical molecules in aqueous solution.
            
        Units:
            - Input: radius_nm (nanometers)
            - Output: D_t_nm2_us (nm²/μs), D_r_rad2_us (rad²/μs)
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
            "ref1_local": list(self.ref1_local),
            "ref2_local": list(self.ref2_local),
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
            ref1_local=np.array(d.get("ref1_local", [1.0, 0.0, 0.0])) if d.get("ref1_local") is not None else np.array([1.0, 0.0, 0.0]),
            ref2_local=np.array(d.get("ref2_local", [0.0, 0.0, 1.0])) if d.get("ref2_local") is not None else np.array([0.0, 0.0, 1.0]),
        )
"""
ionerdss.model.components.instances

Dependencies:
    - numpy: For coordinate arrays and mathematical operations
    - ionerdss.model.components.types: Template classes (MoleculeType, InterfaceType)

See Also:
    ionerdss.model.components.types: Template definitions for molecules and interfaces
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from ionerdss.model.components.types import MoleculeType, InterfaceType

#------------ helper ---------------
# make sure that mol name does not incude _, which is the separater
# for interface name, e.g. A_12 and B would have interface A-12_B_1
def _replace_underscore_with_dash(input_string : str) -> str:
    """Replace underscores with dashes in molecule names.
    
    Ensures molecule names don't contain underscores which are used as
    separators in interface naming conventions (e.g., "A_B_1").
    
    Args:
        input_string: The string to process.
        
    Returns:
        String with all underscores replaced by dashes.
        
    Example:
        "A_12" becomes "A-12" to avoid conflicts with interface naming.
    """
    return input_string.replace("_", "-")
    
# -------------- dataclasses ------------

@dataclass
class InterfaceInstance:
    """Represents a specific instance of a molecular interface in a simulation.
    
    An InterfaceInstance is a concrete realization of an InterfaceType, representing
    an actual binding site on a specific molecule instance with current spatial
    coordinates and binding state. It maintains references to its partner interface
    and parent molecule instance.
    
    Attributes:
        this_mol_name: Name of the molecule containing this interface (e.g., "A").
        partner_mol_name: Name of the target binding partner molecule (e.g., "B").
        interface_index: Numeric index for distinguishing multiple interfaces
            of the same type on one molecule.
        partner_interface: The InterfaceInstance on the partner molecule that
            this interface is bound to or can bind to.
        interface_type: Reference to the InterfaceType template that defines
            this interface's properties.
        this_mol: Reference to the MoleculeInstance containing this interface.
        absolute_coord: Current absolute Cartesian coordinates of the interface.
        residues: List of residue indices that participate in this interface.
        energy: Current binding energy of this interface. Defaults to -1.0.
        signature: Dictionary containing interface-specific runtime properties.
    """
    # Required fields (no defaults) must come first
    absolute_coord: np.ndarray # absolute cartesian coordinates

    # Optional fields with defaults
    partner_interface: Optional['InterfaceInstance'] = None
    interface_type: Optional[InterfaceType] = None
    this_mol: Optional['MoleculeInstance'] = None
    this_mol_name: str = "unnamed"
    partner_mol_name: str = "unnamed"
    interface_index: int = 0
    residues: List[int] = field(default_factory=list)
    energy: Optional[float] = -1.0
    signature: Dict = field(default_factory=dict)

    def _convert_numpy_types(self, obj):
        """Convert NumPy types to native Python types for JSON serialization."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, dict):
            return {key: self._convert_numpy_types(value) for key, value in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [self._convert_numpy_types(item) for item in obj]
        else:
            return obj

    def get_name(self) -> str:
        """Return the formatted interface identifier string.
        
        Constructs the interface name using the format:
        "{this_mol_name}_{partner_mol_name}_{interface_index}"
        
        Returns:
            The interface identifier string (e.g., "A_B_1").
        """
        return self.this_mol_name + "_" +\
            self.partner_mol_name + "_" +\
            str(self.interface_index)

    def set_name(self, new_name: str):
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
        self.this_mol_name = substrings[0]
        self.partner_mol_name = substrings[1]
        self.interface_index = int(substrings[2])

    def to_dict(self) -> dict:
        """Convert the InterfaceInstance to a dictionary representation."""
        data = {
            "name": self.get_name(),
            "type": self.interface_type.get_name() if self.interface_type else "unknown",
            "coord": self.absolute_coord,
            "residues": list(self.residues) if hasattr(self, 'residues') else [],
            "energy": self.energy,
        }
        return self._convert_numpy_types(data)

    @classmethod
    def from_dict(cls, d: dict) -> "InterfaceInstance":
        """Create an InterfaceInstance from a dictionary.
        
        Deserializes interface instance data from a dictionary format, typically
        loaded from JSON or other storage formats. Missing optional fields
        are filled with default values.
        
        Args:
            d: Dictionary containing interface data. Required key is
                'absolute_coords'. Optional keys include 'this_mol_name',
                'partner_mol_name', 'interface_index', 'this_mol',
                'partner_interface', 'interface_type', 'energy', and 'signature'.
        
        Returns:
            New InterfaceInstance populated with the dictionary data.
            
        Raises:
            KeyError: If required key 'absolute_coords' is missing.
        """
        return cls(
            absolute_coord=np.array(d["absolute_coords"]),
            partner_interface=d.get("partner_interface", None),
            interface_type=d.get("interface_type", None),
            this_mol=d.get("this_mol", None),
            this_mol_name=d.get("this_mol_name", "unnamed"),
            partner_mol_name=d.get("partner_mol_name", "unnamed"),
            interface_index=d.get("interface_index", 0),
            residues=d.get("residues", []),
            energy=d.get("energy", -1.0),
            signature=d.get("signature", {}),
        )
        
    # Add these methods to InterfaceInstance class
    def __hash__(self) -> int:
        """Make InterfaceInstance hashable for use as dictionary keys."""
        return hash((self.this_mol_name, self.partner_mol_name, self.interface_index))

    def __eq__(self, other) -> bool:
        """Define equality for InterfaceInstance objects."""
        if not isinstance(other, InterfaceInstance):
            return False
        return (self.this_mol_name == other.this_mol_name and 
                self.partner_mol_name == other.partner_mol_name and 
                self.interface_index == other.interface_index)

@dataclass
class MoleculeInstance:
    """Represents a specific instance of a molecule in a simulation.
    
    A MoleculeInstance is a concrete realization of a MoleculeType, representing
    an actual molecule in the simulation with current position, orientation, and
    binding state. It maintains a map of its interface instances and their
    current binding partners.
    
    Note : in implementation keys of interface_neighbors_map are
    changed from InterfaceInstances into string due to the instances
    containing mutable numpy arrays.
    
    Attributes:
        name: Unique identifier for this molecule instance.
        interfaces_neighbors_map: Dictionary mapping interface strings
            to their currently bound neighbor MoleculeInstance objects.
        molecule_type: Reference to the MoleculeType template that defines
            this molecule's properties.
        norm: Orientation vector (relative coordinates) defining molecular
            rotation state.
        com: Center of mass coordinates (absolute position) of the molecule.
    """
    # Required fields first (no defaults)
    name: str
    norm: np.ndarray
    ref1: np.ndarray
    ref2: np.ndarray
    com: np.ndarray # absolute coord

    # Optional fields with defaults
    molecule_type: Optional[MoleculeType] = None
    interfaces_neighbors_map: Dict = field(default_factory=dict)
    
    def get_name(self):
        """Adding a get_name method to conform to the InterfaceInstance class"""
        return self.name

    def _convert_numpy_types(self, obj):
        """Convert NumPy types to native Python types for JSON serialization."""
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        elif isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, dict):
            return {key: self._convert_numpy_types(value) for key, value in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [self._convert_numpy_types(item) for item in obj]
        else:
            return obj

    # In MoleculeInstance class
    def to_dict(self) -> dict:
        """Convert the MoleculeInstance to a dictionary representation."""
        # Convert interfaces_neighbors_map to serializable format
        interfaces_list = []
        for interface_instance, partner_mol in self.interfaces_neighbors_map.items():
            interfaces_list.append({
                "interface_name": interface_instance.get_name(),
                "partner_molecule": partner_mol.name if partner_mol else "unknown"
            })
        
        data = {
            "name": self.name,
            "type": self.molecule_type.name if self.molecule_type else "unknown",
            "com": self.com,
            "norm": self.norm,
            "ref1": self.ref1,
            "ref2": self.ref2,
            "interfaces": interfaces_list
        }
        return self._convert_numpy_types(data)


    # In MoleculeInstance class
    @classmethod
    def from_dict(cls, d: dict) -> "MoleculeInstance":
        """Create a MoleculeInstance from a dictionary.
        
        Note: The interfaces_neighbors_map will be empty and needs to be 
        rebuilt by the system's _rebuild_cross_references() method.
        """
        return cls(
            name=_replace_underscore_with_dash(d.get("name", "unnamed")),
            norm=np.array(d.get("norm", [0.0, 0.0, 1.0])),
            ref1=np.array(d.get("ref1", [1.0, 0.0, 0.0])),
            ref2=np.array(d.get("ref2", [0.0, 0.0, -1.0])),
            com=np.array(d.get("coord", [0.0, 0.0, 0.0])),
            molecule_type=d.get("molecule_type", None),
            interfaces_neighbors_map={}  # Will be rebuilt by system
        )
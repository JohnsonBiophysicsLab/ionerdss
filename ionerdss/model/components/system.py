"""
ionerdss.model.components.system

System-level container for complete molecular simulation models.

This module provides the System class which serves as the top-level container
for all molecular simulation components including molecule types, instances,
interface definitions, and reaction rules. The System manages the complete
state of a molecular model and provides serialization capabilities for
persistence and data exchange.

Classes:
    System: Complete molecular simulation system containing all registries,
        configuration data, and serialization methods.

Key Features:
    - Centralized management of all molecular components
    - Complete serialization/deserialization to/from JSON
    - Workspace path management for file organization
    - PDB structure integration support
    - Unit system management
    - Registry-based component organization

Architecture Overview:
    The System class acts as the root container that holds:
    - MoleculeTypeRegistry: Templates for all molecule types
    - MoleculeInstanceRegistry: Runtime instances of molecules
    - InterfaceTypeRegistry: Templates for all interface types  
    - InterfaceInstanceRegistry: Runtime instances of interfaces
    - Configuration data: workspace path, PDB ID, units
    
    All data can be serialized to JSON and fully reconstructed, maintaining
    all relationships and references between components.

Example Usage:
    ```python
    # Create a new system
    system = System(
        workspace_path="/path/to/workspace",
        pdb_id="1ABC",
        units=Units()
    )
    
    # Add components to registries
    system.molecule_types.add(protein_a)
    system.interface_types.add(binding_site)
    
    # Save to file
    system.to_json("model.json")
    
    # Load from file
    restored_system = System.from_json("model.json")
    ```

Dependencies:
    - json: For serialization/deserialization
    - pathlib: For path handling
    - All component modules: types, instances, registry, units
"""
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Any, Union

from ionerdss.model.components.units import Units
from ionerdss.model.components.registry import (
    MoleculeTypeRegistry, MoleculeInstanceRegistry,
    InterfaceTypeRegistry, InterfaceInstanceRegistry
)


class System:
    """Complete molecular system containing all components and registries.

    Attributes:
        workspace_path: Path to workspace directory.
        pdb_id: PDB identifier for this system.
        units: Unit system used throughout the system.
        molecule_types: Registry of molecule type definitions.
        molecule_instances: Registry of molecule instances.
        interface_types: Registry of interface type definitions.
        interface_instances: Registry of interface instances.
    """

    def __init__(self, workspace_path: str, pdb_id: Optional[str] = None,
                 units: Optional[Units] = None):
        """Initialize system with empty registries.

        Args:
            workspace_path: Path to workspace directory.
            pdb_id: PDB identifier (optional).
            units: Unit system (defaults to standard units).
        """
        self.workspace_path = workspace_path
        self.pdb_id = pdb_id
        self.units = units or Units()

        # Initialize registries
        self.molecule_types = MoleculeTypeRegistry()
        self.molecule_instances = MoleculeInstanceRegistry()
        self.interface_types = InterfaceTypeRegistry()
        self.interface_instances = InterfaceInstanceRegistry()

    def _convert_numpy_types(self, obj: Any) -> Any:
        """Recursively convert NumPy types to native Python types for JSON serialization.

        Args:
            obj: Object that may contain NumPy types.

        Returns:
            Object with NumPy types converted to native Python types.
        """
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

    def to_dict(self) -> Dict[str, Any]:
        """Convert system to dictionary representation.

        Returns:
            Dictionary containing complete system data.
        """
        system_data = {
            "metadata": {
                "workspace_path": self.workspace_path,
                "pdb_id": self.pdb_id,
                "units": self.units.to_dict()
            },
            "registries": {
                "molecule_types": self.molecule_types.to_list(),
                "molecule_instances": self.molecule_instances.to_list(),
                "interface_types": self.interface_types.to_list(),
                "interface_instances": self.interface_instances.to_list()
            }
        }

        # Convert NumPy types to native Python types
        return self._convert_numpy_types(system_data)

    def to_json(self, filepath: Union[str, Path]) -> None:
        """Save system to JSON file.

        Args:
            filepath: Path to output JSON file.
        """
        filepath = Path(filepath)

        # Get system data with NumPy types converted
        system_data = self.to_dict()

        class SystemJSONEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                elif hasattr(obj, 'to_dict'):
                    return obj.to_dict()
                elif isinstance(obj, set):
                    return list(obj)
                else:
                    return super().default(obj)
        
        system_data = self.to_dict()
        
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(system_data, f, indent=2, ensure_ascii=False, cls=SystemJSONEncoder)


    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "System":
        """Create system from dictionary representation.

        Args:
            data: Dictionary containing system data.

        Returns:
            New System instance.
        """
        # Extract metadata
        metadata = data.get("metadata", {})
        workspace_path = metadata.get("workspace_path", ".")
        pdb_id = metadata.get("pdb_id")
        units = Units.from_dict(metadata.get("units", {}))

        # Create system
        system = cls(workspace_path=workspace_path, pdb_id=pdb_id, units=units)

        # Load registries
        registries = data.get("registries", {})

        system.molecule_types = MoleculeTypeRegistry.from_list(
            registries.get("molecule_types", [])
        )
        system.molecule_instances = MoleculeInstanceRegistry.from_list(
            registries.get("molecule_instances", [])
        )
        system.interface_types = InterfaceTypeRegistry.from_list(
            registries.get("interface_types", [])
        )
        system.interface_instances = InterfaceInstanceRegistry.from_list(
            registries.get("interface_instances", [])
        )

        # Rebuild cross-references
        system._rebuild_cross_references()

        return system

    @classmethod
    def from_json(cls, filepath: Union[str, Path]) -> "System":
        """Load system from JSON file.

        Args:
            filepath: Path to JSON file.

        Returns:
            New System instance.
        """
        filepath = Path(filepath)

        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)

        return cls.from_dict(data)

    def _rebuild_cross_references(self) -> None:
        """Rebuild cross-references between components after loading.

        Internal method that restores object references between components
        that were lost during JSON serialization. This includes linking:
        - InterfaceTypes to their parent MoleculeTypes
        - InterfaceTypes to their partner InterfaceTypes
        - MoleculeInstances to their MoleculeTypes
        - InterfaceInstances to their InterfaceTypes and MoleculeInstances

        This method is called automatically by from_json() and should not
        be called directly by user code.
        """
        # Rebuild molecule type references in interface types
        for interface_type in self.interface_types:
            # Link to parent molecule type
            if interface_type.this_mol_type_name in self.molecule_types:
                interface_type.this_mol_type = self.molecule_types.get(
                    interface_type.this_mol_type_name
                )

            # Link to partner molecule type
            if interface_type.partner_mol_type_name in self.molecule_types:
                interface_type.partner_mol_type = self.molecule_types.get(
                    interface_type.partner_mol_type_name
                )

        # Rebuild partner interface references
        for interface_type in self.interface_types:
            partner_name = f"{interface_type.partner_mol_type_name}_{interface_type.this_mol_type_name}_{interface_type.interface_index}"
            if partner_name in self.interface_types:
                interface_type.partner_interface_type = self.interface_types.get(
                    partner_name)

        # Rebuild molecule instance references to their types
        for mol_instance in self.molecule_instances:
            if hasattr(mol_instance, 'molecule_type') and mol_instance.molecule_type:
                type_name = mol_instance.molecule_type.name if hasattr(
                    mol_instance.molecule_type, 'name') else str(mol_instance.molecule_type)
                if type_name in self.molecule_types:
                    mol_instance.molecule_type = self.molecule_types.get(
                        type_name)

        # Rebuild interface instance references
        for interface_instance in self.interface_instances:
            # Link to interface type
            interface_name = interface_instance.get_name()
            if interface_name in self.interface_types:
                interface_instance.interface_type = self.interface_types.get(
                    interface_name)

            # Link to parent molecule instance
            if hasattr(interface_instance, 'this_mol_name'):
                if interface_instance.this_mol_name in self.molecule_instances:
                    interface_instance.this_mol = self.molecule_instances.get(
                        interface_instance.this_mol_name
                    )

    def validate_system(self) -> Dict[str, List[str]]:
        """Validate system integrity.

        Returns:
            Dictionary with 'errors' and 'warnings' lists.
        """
        errors = []
        warnings = []

        # Check that we have some components
        if len(self.molecule_types.molecule_types) == 0:
            errors.append("No molecule types defined")

        if len(self.molecule_instances.molecule_instances) == 0:
            warnings.append("No molecule instances defined")

        # Check interface consistency
        for interface_type in self.interface_types.interface_types.values():
            # Check that referenced molecule types exist
            if interface_type.this_mol_type_name not in self.molecule_types.molecule_types:
                errors.append(
                    f"Interface type {interface_type.get_name()} references unknown molecule type: {interface_type.this_mol_type_name}")

            if interface_type.partner_mol_type_name not in self.molecule_types.molecule_types:
                errors.append(
                    f"Interface type {interface_type.get_name()} references unknown partner molecule type: {interface_type.partner_mol_type_name}")

        # Check instance consistency
        for mol_instance in self.molecule_instances.molecule_instances.values():
            if mol_instance.molecule_type and mol_instance.molecule_type.name not in self.molecule_types.molecule_types:
                errors.append(
                    f"Molecule instance {mol_instance.name} references unknown molecule type: {mol_instance.molecule_type.name}")

        return {
            "errors": errors,
            "warnings": warnings
        }

    def get_summary(self) -> Dict[str, Any]:
        """Get summary statistics of the system.

        Returns:
            Dictionary with system statistics.
        """
        return {
            "workspace_path": self.workspace_path,
            "pdb_id": self.pdb_id,
            "units": self.units.to_dict(),
            "molecule_types_count": len(self.molecule_types.molecule_types),
            "molecule_instances_count": len(self.molecule_instances.molecule_instances),
            "interface_types_count": len(self.interface_types.interface_types),
            "interface_instances_count": len(self.interface_instances.interface_instances),
            "molecule_type_names": list(self.molecule_types.molecule_types.keys()),
            "interface_type_names": list(self.interface_types.interface_types.keys())
        }

    def __str__(self) -> str:
        """String representation of the system."""
        return (f"System(pdb_id={self.pdb_id}, "
                f"molecule_types={len(self.molecule_types.molecule_types)}, "
                f"molecule_instances={len(self.molecule_instances.molecule_instances)}, "
                f"interface_types={len(self.interface_types.interface_types)}, "
                f"interface_instances={len(self.interface_instances.interface_instances)})")

    def __repr__(self) -> str:
        """Detailed representation of the system."""
        return self.__str__()

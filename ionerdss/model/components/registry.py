"""
ionerdss.model.components.registry

Registry classes for managing collections of molecular components.

This module provides registry classes that act as centralized collections
for managing molecule types, molecule instances, interface types, and
interface instances. Each registry provides methods for adding, retrieving,
and serializing collections of related components.

Classes:
    MoleculeTypeRegistry: Registry for MoleculeType objects
    MoleculeInstanceRegistry: Registry for MoleculeInstance objects
    InterfaceTypeRegistry: Registry for InterfaceType objects
    InterfaceInstanceRegistry: Registry for InterfaceInstance objects

Key Features:
    - Centralized component management
    - Name-based lookup and storage
    - Iteration support for easy traversal
    - Serialization/deserialization capabilities
    - Membership testing with 'in' operator
    - List conversion for bulk operations

Example:
    ```python
    # Create and populate a molecule type registry
    mol_registry = MoleculeTypeRegistry()
    
    protein_a = MoleculeType(name="ProteinA", radius_nm=2.5)
    protein_b = MoleculeType(name="ProteinB", radius_nm=3.0)
    
    mol_registry.add(protein_a)
    mol_registry.add(protein_b)
    
    # Iterate through all molecule types
    for mol_type in mol_registry:
        print(f"Molecule: {mol_type.name}, Radius: {mol_type.radius_nm}")
    
    # Check membership and retrieve
    if "ProteinA" in mol_registry:
        retrieved = mol_registry.get("ProteinA")
    
    # Serialize for storage
    serialized = mol_registry.to_list()
    restored_registry = MoleculeTypeRegistry.from_list(serialized)
    ```

See Also:
    ionerdss.model.components.types: Type definitions for molecules and interfaces
    ionerdss.model.components.instances: Instance classes for simulation runtime
"""

from __future__ import annotations
from typing import Dict, List, Iterator

from ionerdss.model.components.types import (
    MoleculeType,
    InterfaceType,
)

from ionerdss.model.components.instances import (
    MoleculeInstance,
    InterfaceInstance,
)


class MoleculeTypeRegistry:
    """Registry for managing collections of MoleculeType objects.

    Provides centralized storage and retrieval of molecule type definitions
    with name-based indexing and iteration capabilities.

    Attributes:
        molecule_types: Dictionary mapping molecule names to MoleculeType objects.
    """

    def __init__(self) -> None:
        """Initialize an empty molecule type registry."""
        self.molecule_types: Dict[str, MoleculeType] = {}

    def add(self, mt: MoleculeType) -> None:
        """Add a molecule type to the registry.

        Args:
            mt: MoleculeType object to add to the registry.

        Raises:
            ValueError: If a molecule type with the same name already exists.
        """
        if mt.name in self.molecule_types:
            raise ValueError(
                f"MoleculeType '{mt.name}' already exists in registry")
        self.molecule_types[mt.name] = mt

    def get(self, name: str) -> MoleculeType:
        """Retrieve a molecule type by name.

        Args:
            name: Name of the molecule type to retrieve.

        Returns:
            The MoleculeType object with the specified name.

        Raises:
            KeyError: If no molecule type with the given name exists.
        """
        if name not in self.molecule_types:
            raise KeyError(f"MoleculeType '{name}' not found in registry")
        return self.molecule_types[name]

    def remove(self, name: str) -> MoleculeType:
        """Remove and return a molecule type from the registry.

        Args:
            name: Name of the molecule type to remove.

        Returns:
            The removed MoleculeType object.

        Raises:
            KeyError: If no molecule type with the given name exists.
        """
        if name not in self.molecule_types:
            raise KeyError(f"MoleculeType '{name}' not found in registry")
        return self.molecule_types.pop(name)

    def clear(self) -> None:
        """Remove all molecule types from the registry."""
        self.molecule_types.clear()

    def __contains__(self, name: str) -> bool:
        """Check if a molecule type exists in the registry.

        Args:
            name: Name of the molecule type to check.

        Returns:
            True if the molecule type exists, False otherwise.
        """
        return name in self.molecule_types

    def __iter__(self) -> Iterator[MoleculeType]:
        """Iterate over all molecule types in the registry.

        Returns:
            Iterator over MoleculeType objects.
        """
        return iter(self.molecule_types.values())

    def __len__(self) -> int:
        """Return the number of molecule types in the registry.

        Returns:
            Number of molecule types stored.
        """
        return len(self.molecule_types)

    def __repr__(self) -> str:
        """Return string representation of the registry.

        Returns:
            String describing the registry contents.
        """
        names = list(self.molecule_types.keys())
        return f"MoleculeTypeRegistry({len(names)} types: {names})"

    def to_list(self) -> List[dict]:
        """Convert all molecule types to a list of dictionaries.

        Returns:
            List of dictionaries representing each molecule type.
        """
        return [mt.to_dict() for mt in self.molecule_types.values()]

    @classmethod
    def from_list(cls, items: List[dict]) -> "MoleculeTypeRegistry":
        """Create a registry from a list of molecule type dictionaries.

        Args:
            items: List of dictionaries representing molecule types.
                   Can be None or empty.

        Returns:
            New MoleculeTypeRegistry populated with the provided data.
        """
        reg = cls()
        for d in items or []:
            reg.add(MoleculeType.from_dict(d))
        return reg


class MoleculeInstanceRegistry:
    """Registry for managing collections of MoleculeInstance objects.

    Provides centralized storage and retrieval of molecule instances
    with name-based indexing and iteration capabilities.

    Attributes:
        molecule_instances: Dictionary mapping instance names to MoleculeInstance objects.
    """

    def __init__(self) -> None:
        """Initialize an empty molecule instance registry."""
        self.molecule_instances: Dict[str, MoleculeInstance] = {}

    def add(self, mi: MoleculeInstance) -> None:
        """Add a molecule instance to the registry.

        Args:
            mi: MoleculeInstance object to add to the registry.

        Raises:
            ValueError: If a molecule instance with the same name already exists.
        """
        if mi.name in self.molecule_instances:
            raise ValueError(
                f"MoleculeInstance '{mi.name}' already exists in registry")
        self.molecule_instances[mi.name] = mi

    def get(self, name: str) -> MoleculeInstance:
        """Retrieve a molecule instance by name.

        Args:
            name: Name of the molecule instance to retrieve.

        Returns:
            The MoleculeInstance object with the specified name.

        Raises:
            KeyError: If no molecule instance with the given name exists.
        """
        if name not in self.molecule_instances:
            raise KeyError(f"MoleculeInstance '{name}' not found in registry")
        return self.molecule_instances[name]

    def remove(self, name: str) -> MoleculeInstance:
        """Remove and return a molecule instance from the registry.

        Args:
            name: Name of the molecule instance to remove.

        Returns:
            The removed MoleculeInstance object.

        Raises:
            KeyError: If no molecule instance with the given name exists.
        """
        if name not in self.molecule_instances:
            raise KeyError(f"MoleculeInstance '{name}' not found in registry")
        return self.molecule_instances.pop(name)

    def clear(self) -> None:
        """Remove all molecule instances from the registry."""
        self.molecule_instances.clear()

    def __contains__(self, name: str) -> bool:
        """Check if a molecule instance exists in the registry.

        Args:
            name: Name of the molecule instance to check.

        Returns:
            True if the molecule instance exists, False otherwise.
        """
        return name in self.molecule_instances

    def __iter__(self) -> Iterator[MoleculeInstance]:
        """Iterate over all molecule instances in the registry.

        Returns:
            Iterator over MoleculeInstance objects.
        """
        return iter(self.molecule_instances.values())

    def __len__(self) -> int:
        """Return the number of molecule instances in the registry.

        Returns:
            Number of molecule instances stored.
        """
        return len(self.molecule_instances)

    def __repr__(self) -> str:
        """Return string representation of the registry.

        Returns:
            String describing the registry contents.
        """
        names = list(self.molecule_instances.keys())
        return f"MoleculeInstanceRegistry({len(names)} instances: {names})"

    def to_list(self) -> List[dict]:
        """Convert all molecule instances to a list of dictionaries.

        Returns:
            List of dictionaries representing each molecule instance.
        """
        return [mi.to_dict() for mi in self.molecule_instances.values()]

    @classmethod
    def from_list(cls, items: List[dict]) -> "MoleculeInstanceRegistry":
        """Create a registry from a list of molecule instance dictionaries.

        Args:
            items: List of dictionaries representing molecule instances.
                   Can be None or empty.

        Returns:
            New MoleculeInstanceRegistry populated with the provided data.
        """
        reg = cls()
        for d in items or []:
            reg.add(MoleculeInstance.from_dict(d))
        return reg


class InterfaceTypeRegistry:
    """Registry for managing collections of InterfaceType objects.

    Provides centralized storage and retrieval of interface type definitions
    with name-based indexing and iteration capabilities.

    Attributes:
        interface_types: Dictionary mapping interface names to InterfaceType objects.
    """

    def __init__(self) -> None:
        """Initialize an empty interface type registry."""
        self.interface_types: Dict[str, InterfaceType] = {}

    def add(self, it: InterfaceType) -> None:
        """Add an interface type to the registry.

        Args:
            it: InterfaceType object to add to the registry.

        Raises:
            ValueError: If an interface type with the same name already exists.
        """
        name = it.get_name()
        if name in self.interface_types:
            raise ValueError(
                f"InterfaceType '{name}' already exists in registry")
        self.interface_types[name] = it

    def get(self, name: str) -> InterfaceType:
        """Retrieve an interface type by name.

        Args:
            name: Name of the interface type to retrieve.

        Returns:
            The InterfaceType object with the specified name.

        Raises:
            KeyError: If no interface type with the given name exists.
        """
        if name not in self.interface_types:
            raise KeyError(f"InterfaceType '{name}' not found in registry")
        return self.interface_types[name]

    def remove(self, name: str) -> InterfaceType:
        """Remove and return an interface type from the registry.

        Args:
            name: Name of the interface type to remove.

        Returns:
            The removed InterfaceType object.

        Raises:
            KeyError: If no interface type with the given name exists.
        """
        if name not in self.interface_types:
            raise KeyError(f"InterfaceType '{name}' not found in registry")
        return self.interface_types.pop(name)

    def clear(self) -> None:
        """Remove all interface types from the registry."""
        self.interface_types.clear()

    def __contains__(self, name: str) -> bool:
        """Check if an interface type exists in the registry.

        Args:
            name: Name of the interface type to check.

        Returns:
            True if the interface type exists, False otherwise.
        """
        return name in self.interface_types

    def __iter__(self) -> Iterator[InterfaceType]:
        """Iterate over all interface types in the registry.

        Returns:
            Iterator over InterfaceType objects.
        """
        return iter(self.interface_types.values())

    def __len__(self) -> int:
        """Return the number of interface types in the registry.

        Returns:
            Number of interface types stored.
        """
        return len(self.interface_types)

    def __repr__(self) -> str:
        """Return string representation of the registry.

        Returns:
            String describing the registry contents.
        """
        names = list(self.interface_types.keys())
        return f"InterfaceTypeRegistry({len(names)} types: {names})"

    def to_list(self) -> List[dict]:
        """Convert all interface types to a list of dictionaries.

        Returns:
            List of dictionaries representing each interface type.
        """
        return [it.to_dict() for it in self.interface_types.values()]

    @classmethod
    def from_list(cls, items: List[dict]) -> "InterfaceTypeRegistry":
        """Create a registry from a list of interface type dictionaries.

        Args:
            items: List of dictionaries representing interface types.
                   Can be None or empty.

        Returns:
            New InterfaceTypeRegistry populated with the provided data.
        """
        reg = cls()
        for d in items or []:
            reg.add(InterfaceType.from_dict(d))
        return reg


class InterfaceInstanceRegistry:
    """Registry for managing collections of InterfaceInstance objects.

    Provides centralized storage and retrieval of interface instances
    with name-based indexing and iteration capabilities.

    Attributes:
        interface_instances: Dictionary mapping instance names to InterfaceInstance objects.
    """

    def __init__(self) -> None:
        """Initialize an empty interface instance registry."""
        self.interface_instances: Dict[str, InterfaceInstance] = {}

    def add(self, ii: InterfaceInstance) -> None:
        """Add an interface instance to the registry.

        Args:
            ii: InterfaceInstance object to add to the registry.

        Raises:
            ValueError: If an interface instance with the same name already exists.
        """
        name = ii.get_name()
        if name in self.interface_instances:
            raise ValueError(
                f"InterfaceInstance '{name}' already exists in registry")
        self.interface_instances[name] = ii

    def get(self, name: str) -> InterfaceInstance:
        """Retrieve an interface instance by name.

        Args:
            name: Name of the interface instance to retrieve.

        Returns:
            The InterfaceInstance object with the specified name.

        Raises:
            KeyError: If no interface instance with the given name exists.
        """
        if name not in self.interface_instances:
            raise KeyError(f"InterfaceInstance '{name}' not found in registry")
        return self.interface_instances[name]

    def remove(self, name: str) -> InterfaceInstance:
        """Remove and return an interface instance from the registry.

        Args:
            name: Name of the interface instance to remove.

        Returns:
            The removed InterfaceInstance object.

        Raises:
            KeyError: If no interface instance with the given name exists.
        """
        if name not in self.interface_instances:
            raise KeyError(f"InterfaceInstance '{name}' not found in registry")
        return self.interface_instances.pop(name)

    def clear(self) -> None:
        """Remove all interface instances from the registry."""
        self.interface_instances.clear()

    def __contains__(self, name: str) -> bool:
        """Check if an interface instance exists in the registry.

        Args:
            name: Name of the interface instance to check.

        Returns:
            True if the interface instance exists, False otherwise.
        """
        return name in self.interface_instances

    def __iter__(self) -> Iterator[InterfaceInstance]:
        """Iterate over all interface instances in the registry.

        Returns:
            Iterator over InterfaceInstance objects.
        """
        return iter(self.interface_instances.values())

    def __len__(self) -> int:
        """Return the number of interface instances in the registry.

        Returns:
            Number of interface instances stored.
        """
        return len(self.interface_instances)

    def __repr__(self) -> str:
        """Return string representation of the registry.

        Returns:
            String describing the registry contents.
        """
        names = list(self.interface_instances.keys())
        return f"InterfaceInstanceRegistry({len(names)} instances: {names})"

    def to_list(self) -> List[dict]:
        """Convert all interface instances to a list of dictionaries.

        Returns:
            List of dictionaries representing each interface instance.
        """
        return [ii.to_dict() for ii in self.interface_instances.values()]

    @classmethod
    def from_list(cls, items: List[dict]) -> "InterfaceInstanceRegistry":
        """Create a registry from a list of interface instance dictionaries.

        Args:
            items: List of dictionaries representing interface instances.
                   Can be None or empty.

        Returns:
            New InterfaceInstanceRegistry populated with the provided data.
        """
        reg = cls()
        for d in items or []:
            reg.add(InterfaceInstance.from_dict(d))
        return reg

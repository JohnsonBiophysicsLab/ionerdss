"""
Unit tests for ionerdss.model.components.registry module.

Tests all registry classes including their iteration capabilities,
serialization, and error handling.
"""

import unittest
import numpy as np

# Import the modules under test
from ionerdss.model.components.registry import (
    MoleculeTypeRegistry,
    MoleculeInstanceRegistry,
    InterfaceTypeRegistry,
    InterfaceInstanceRegistry
)
from ionerdss.model.components.types import MoleculeType, InterfaceType
from ionerdss.model.components.instances import MoleculeInstance, InterfaceInstance


class TestMoleculeTypeRegistry(unittest.TestCase):
    """Test the MoleculeTypeRegistry class."""

    def setUp(self):
        """Set up test fixtures."""
        self.registry = MoleculeTypeRegistry()
        self.mol_type_a = MoleculeType(name="ProteinA", radius_nm=2.5)
        self.mol_type_b = MoleculeType(name="ProteinB", radius_nm=3.0)

    def test_init(self):
        """Test registry initialization."""
        registry = MoleculeTypeRegistry()
        self.assertEqual(len(registry), 0)
        self.assertEqual(len(registry.molecule_types), 0)

    def test_add_single(self):
        """Test adding a single molecule type."""
        self.registry.add(self.mol_type_a)

        self.assertEqual(len(self.registry), 1)
        self.assertIn("ProteinA", self.registry)
        self.assertEqual(self.registry.get("ProteinA"), self.mol_type_a)

    def test_add_multiple(self):
        """Test adding multiple molecule types."""
        self.registry.add(self.mol_type_a)
        self.registry.add(self.mol_type_b)

        self.assertEqual(len(self.registry), 2)
        self.assertIn("ProteinA", self.registry)
        self.assertIn("ProteinB", self.registry)

    def test_add_duplicate_raises_error(self):
        """Test that adding duplicate names raises ValueError."""
        self.registry.add(self.mol_type_a)

        duplicate = MoleculeType(name="ProteinA", radius_nm=1.0)
        with self.assertRaises(ValueError) as context:
            self.registry.add(duplicate)

        self.assertIn("already exists", str(context.exception))

    def test_get_existing(self):
        """Test retrieving existing molecule type."""
        self.registry.add(self.mol_type_a)
        retrieved = self.registry.get("ProteinA")

        self.assertEqual(retrieved, self.mol_type_a)
        self.assertEqual(retrieved.name, "ProteinA")

    def test_get_nonexistent_raises_error(self):
        """Test that getting nonexistent molecule type raises KeyError."""
        with self.assertRaises(KeyError) as context:
            self.registry.get("NonExistent")

        self.assertIn("not found", str(context.exception))

    def test_remove_existing(self):
        """Test removing existing molecule type."""
        self.registry.add(self.mol_type_a)
        self.registry.add(self.mol_type_b)

        removed = self.registry.remove("ProteinA")

        self.assertEqual(removed, self.mol_type_a)
        self.assertEqual(len(self.registry), 1)
        self.assertNotIn("ProteinA", self.registry)
        self.assertIn("ProteinB", self.registry)

    def test_remove_nonexistent_raises_error(self):
        """Test that removing nonexistent molecule type raises KeyError."""
        with self.assertRaises(KeyError) as context:
            self.registry.remove("NonExistent")

        self.assertIn("not found", str(context.exception))

    def test_clear(self):
        """Test clearing all molecule types."""
        self.registry.add(self.mol_type_a)
        self.registry.add(self.mol_type_b)

        self.registry.clear()

        self.assertEqual(len(self.registry), 0)
        self.assertNotIn("ProteinA", self.registry)
        self.assertNotIn("ProteinB", self.registry)

    def test_contains(self):
        """Test membership testing with 'in' operator."""
        self.assertNotIn("ProteinA", self.registry)

        self.registry.add(self.mol_type_a)
        self.assertIn("ProteinA", self.registry)
        self.assertNotIn("ProteinB", self.registry)

    def test_iteration(self):
        """Test iterating over molecule types."""
        self.registry.add(self.mol_type_a)
        self.registry.add(self.mol_type_b)

        # Test iteration
        types = list(self.registry)
        self.assertEqual(len(types), 2)
        self.assertIn(self.mol_type_a, types)
        self.assertIn(self.mol_type_b, types)

        # Test that we can iterate multiple times
        count = sum(1 for _ in self.registry)
        self.assertEqual(count, 2)

    def test_iteration_empty(self):
        """Test iterating over empty registry."""
        types = list(self.registry)
        self.assertEqual(len(types), 0)

    def test_len(self):
        """Test length reporting."""
        self.assertEqual(len(self.registry), 0)

        self.registry.add(self.mol_type_a)
        self.assertEqual(len(self.registry), 1)

        self.registry.add(self.mol_type_b)
        self.assertEqual(len(self.registry), 2)

        self.registry.remove("ProteinA")
        self.assertEqual(len(self.registry), 1)

    def test_repr(self):
        """Test string representation."""
        # Empty registry
        repr_str = repr(self.registry)
        self.assertIn("MoleculeTypeRegistry", repr_str)
        self.assertIn("0 types", repr_str)

        # Registry with items
        self.registry.add(self.mol_type_a)
        self.registry.add(self.mol_type_b)
        repr_str = repr(self.registry)
        self.assertIn("MoleculeTypeRegistry", repr_str)
        self.assertIn("2 types", repr_str)
        self.assertIn("ProteinA", repr_str)
        self.assertIn("ProteinB", repr_str)

    def test_to_list(self):
        """Test conversion to list of dictionaries."""
        self.registry.add(self.mol_type_a)
        self.registry.add(self.mol_type_b)

        result = self.registry.to_list()

        self.assertEqual(len(result), 2)
        self.assertIsInstance(result, list)
        self.assertIsInstance(result[0], dict)

        # Check that all names are present
        names = [item["name"] for item in result]
        self.assertIn("ProteinA", names)
        self.assertIn("ProteinB", names)

    def test_from_list_empty(self):
        """Test creating registry from empty list."""
        registry = MoleculeTypeRegistry.from_list([])
        self.assertEqual(len(registry), 0)

        registry_none = MoleculeTypeRegistry.from_list(None)
        self.assertEqual(len(registry_none), 0)

    def test_from_list_populated(self):
        """Test creating registry from populated list."""
        data = [
            {"name": "TestA", "radius": 1.0, "diffusion_translation": 0.5,
             "diffusion_rotation": 0.25, "interfaces_neighbors_map": {}},
            {"name": "TestB", "radius": 2.0, "diffusion_translation": 0.3,
             "diffusion_rotation": 0.15, "interfaces_neighbors_map": {}}
        ]

        registry = MoleculeTypeRegistry.from_list(data)

        self.assertEqual(len(registry), 2)
        self.assertIn("TestA", registry)
        self.assertIn("TestB", registry)

        mol_a = registry.get("TestA")
        self.assertEqual(mol_a.name, "TestA")
        self.assertEqual(mol_a.radius_nm, 1.0)

    def test_serialization_roundtrip(self):
        """Test serialization and deserialization roundtrip."""
        self.registry.add(self.mol_type_a)
        self.registry.add(self.mol_type_b)

        # Serialize
        data = self.registry.to_list()

        # Deserialize
        new_registry = MoleculeTypeRegistry.from_list(data)

        # Verify
        self.assertEqual(len(new_registry), len(self.registry))
        for mol_type in self.registry:
            self.assertIn(mol_type.name, new_registry)
            retrieved = new_registry.get(mol_type.name)
            self.assertEqual(retrieved.name, mol_type.name)


class TestMoleculeInstanceRegistry(unittest.TestCase):
    """Test the MoleculeInstanceRegistry class."""

    def setUp(self):
        """Set up test fixtures."""
        self.registry = MoleculeInstanceRegistry()
        self.mol_instance_a = MoleculeInstance(
            name="MolA_001",
            norm=np.array([0.0, 0.0, 1.0]),
            ref1=np.array([1.0, 0.0, 0.0]),
            ref2=np.array([0.0, 0.0, 1.0]),
            com=np.array([0.0, 0.0, 0.0])
        )
        self.mol_instance_b = MoleculeInstance(
            name="MolB_001",
            norm=np.array([1.0, 0.0, 0.0]),
            ref1=np.array([0.0, 1.0, 0.0]),
            ref2=np.array([0.0, 0.0, 1.0]),
            com=np.array([1.0, 1.0, 1.0])
        )

    def test_basic_operations(self):
        """Test basic registry operations."""
        # Add
        self.registry.add(self.mol_instance_a)
        self.assertEqual(len(self.registry), 1)
        self.assertIn("MolA_001", self.registry)

        # Get
        retrieved = self.registry.get("MolA_001")
        self.assertEqual(retrieved, self.mol_instance_a)

        # Remove
        removed = self.registry.remove("MolA_001")
        self.assertEqual(removed, self.mol_instance_a)
        self.assertEqual(len(self.registry), 0)

    def test_iteration(self):
        """Test iteration over molecule instances."""
        self.registry.add(self.mol_instance_a)
        self.registry.add(self.mol_instance_b)

        instances = list(self.registry)
        self.assertEqual(len(instances), 2)
        self.assertIn(self.mol_instance_a, instances)
        self.assertIn(self.mol_instance_b, instances)

    def test_repr(self):
        """Test string representation."""
        self.registry.add(self.mol_instance_a)
        repr_str = repr(self.registry)

        self.assertIn("MoleculeInstanceRegistry", repr_str)
        self.assertIn("1 instances", repr_str)
        self.assertIn("MolA_001", repr_str)


class TestInterfaceTypeRegistry(unittest.TestCase):
    """Test the InterfaceTypeRegistry class."""

    def setUp(self):
        """Set up test fixtures."""
        self.registry = InterfaceTypeRegistry()
        self.interface_type_a = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            local_coord=np.array([0.5, 0.0, 0.0])
        )
        self.interface_type_b = InterfaceType(
            this_mol_type_name="B",
            partner_mol_type_name="A",
            interface_index=1,
            absolute_coord=np.array([-1.0, 0.0, 0.0]),
            local_coord=np.array([-0.5, 0.0, 0.0])
        )

    def test_add_and_get_by_name(self):
        """Test adding and retrieving by generated name."""
        self.registry.add(self.interface_type_a)

        name = self.interface_type_a.get_name()  # Should be "AB1"
        self.assertIn(name, self.registry)

        retrieved = self.registry.get(name)
        self.assertEqual(retrieved, self.interface_type_a)

    def test_iteration(self):
        """Test iteration over interface types."""
        self.registry.add(self.interface_type_a)
        self.registry.add(self.interface_type_b)

        types = list(self.registry)
        self.assertEqual(len(types), 2)
        self.assertIn(self.interface_type_a, types)
        self.assertIn(self.interface_type_b, types)

    def test_repr(self):
        """Test string representation."""
        self.registry.add(self.interface_type_a)
        repr_str = repr(self.registry)

        self.assertIn("InterfaceTypeRegistry", repr_str)
        self.assertIn("1 types", repr_str)
        self.assertIn("AB1", repr_str)

    def test_duplicate_interface_names(self):
        """Test handling of duplicate interface names."""
        self.registry.add(self.interface_type_a)

        # Create another interface with same name components
        duplicate = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([2.0, 0.0, 0.0]),
            local_coord=np.array([1.0, 0.0, 0.0])
        )

        with self.assertRaises(ValueError) as context:
            self.registry.add(duplicate)

        self.assertIn("already exists", str(context.exception))


class TestInterfaceInstanceRegistry(unittest.TestCase):
    """Test the InterfaceInstanceRegistry class."""

    def setUp(self):
        """Set up test fixtures."""
        self.registry = InterfaceInstanceRegistry()
        self.interface_instance_a = InterfaceInstance(
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            this_mol_name="A",
            partner_mol_name="B",
            interface_index=1
        )
        self.interface_instance_b = InterfaceInstance(
            absolute_coord=np.array([-1.0, 0.0, 0.0]),
            this_mol_name="B",
            partner_mol_name="A",
            interface_index=1
        )

    def test_basic_operations(self):
        """Test basic registry operations."""
        # Add
        self.registry.add(self.interface_instance_a)
        self.assertEqual(len(self.registry), 1)

        name = self.interface_instance_a.get_name()  # "A_B_1"
        self.assertIn(name, self.registry)

        # Get
        retrieved = self.registry.get(name)
        self.assertEqual(retrieved.get_name(),
                         self.interface_instance_a.get_name())

        # Remove
        removed = self.registry.remove(name)
        self.assertEqual(removed.get_name(),
                         self.interface_instance_a.get_name())
        self.assertEqual(len(self.registry), 0)

    def test_iteration(self):
        """Test iteration over interface instances."""
        self.registry.add(self.interface_instance_a)
        self.registry.add(self.interface_instance_b)

        instances = list(self.registry)
        self.assertEqual(len(instances), 2)

        # Compare by names instead of direct object comparison
        instance_names = [inst.get_name() for inst in instances]
        self.assertIn("A_B_1", instance_names)
        self.assertIn("B_A_1", instance_names)

    def test_repr(self):
        """Test string representation."""
        self.registry.add(self.interface_instance_a)
        repr_str = repr(self.registry)

        self.assertIn("InterfaceInstanceRegistry", repr_str)
        self.assertIn("1 instances", repr_str)
        self.assertIn("A_B_1", repr_str)


class TestRegistryIteration(unittest.TestCase):
    """Test iteration capabilities across all registry types."""

    def test_molecule_type_registry_iteration_patterns(self):
        """Test various iteration patterns for MoleculeTypeRegistry."""
        registry = MoleculeTypeRegistry()

        # Add test data
        mol_types = [
            MoleculeType(name=f"Protein{i}", radius_nm=float(i))
            for i in range(5)
        ]

        for mol_type in mol_types:
            registry.add(mol_type)

        # Test list comprehension
        names = [mt.name for mt in registry]
        self.assertEqual(len(names), 5)
        self.assertIn("Protein0", names)
        self.assertIn("Protein4", names)

        # Test filtering during iteration
        large_proteins = [mt for mt in registry if mt.radius_nm > 2.0]
        self.assertEqual(len(large_proteins), 2)  # Protein3 and Protein4

        # Test enumerate
        for i, mol_type in enumerate(registry):
            self.assertIsInstance(mol_type, MoleculeType)
            if i >= 4:  # Should only have 5 items
                break
        else:
            # This should execute (no break occurred)
            pass

        # Test sum with generator expression
        total_radius = sum(mt.radius_nm for mt in registry)
        expected_total = sum(float(i) for i in range(5))  # 0+1+2+3+4 = 10
        self.assertEqual(total_radius, expected_total)

    def test_molecule_instance_registry_iteration_patterns(self):
        """Test various iteration patterns for MoleculeInstanceRegistry."""
        registry = MoleculeInstanceRegistry()

        # Add test data
        for i in range(3):
            instance = MoleculeInstance(
                name=f"Mol_{i:03d}",
                norm=np.array([0.0, 0.0, 1.0]),
                ref1=np.array([1.0, 0.0, 0.0]),
                ref2=np.array([0.0, 0.0, 1.0]),
                com=np.array([float(i), 0.0, 0.0])
            )
            registry.add(instance)

        # Test iteration with coordinate access
        x_positions = [mi.com[0] for mi in registry]
        self.assertEqual(len(x_positions), 3)
        self.assertIn(0.0, x_positions)
        self.assertIn(1.0, x_positions)
        self.assertIn(2.0, x_positions)

        # Test any/all operations
        all_have_names = all(mi.name.startswith("Mol_") for mi in registry)
        self.assertTrue(all_have_names)

        any_at_origin = any(np.allclose(
            mi.com, [0.0, 0.0, 0.0]) for mi in registry)
        self.assertTrue(any_at_origin)

    def test_interface_type_registry_iteration_patterns(self):
        """Test various iteration patterns for InterfaceTypeRegistry."""
        registry = InterfaceTypeRegistry()

        # Add test data
        interface_configs = [
            ("A", "B", 1), ("B", "A", 1), ("A", "C", 1), ("C", "A", 1)
        ]

        for this_mol, partner_mol, index in interface_configs:
            interface = InterfaceType(
                this_mol_type_name=this_mol,
                partner_mol_type_name=partner_mol,
                interface_index=index,
                absolute_coord=np.array([1.0, 0.0, 0.0]),
                local_coord=np.array([0.5, 0.0, 0.0])
            )
            registry.add(interface)

        # Test filtering by molecule type
        a_interfaces = [it for it in registry if it.this_mol_type_name == "A"]
        self.assertEqual(len(a_interfaces), 2)  # AB1 and AC1

        # Test grouping by partner
        partners = set(it.partner_mol_type_name for it in registry)
        self.assertEqual(partners, {"A", "B", "C"})

    def test_interface_instance_registry_iteration_patterns(self):
        """Test various iteration patterns for InterfaceInstanceRegistry."""
        registry = InterfaceInstanceRegistry()

        # Add test data
        for i, (this_mol, partner_mol) in enumerate([("A", "B"), ("B", "A"), ("A", "C")]):
            instance = InterfaceInstance(
                absolute_coord=np.array([float(i), 0.0, 0.0]),
                this_mol_name=this_mol,
                partner_mol_name=partner_mol,
                interface_index=1,
                energy=-float(i + 1)  # -1, -2, -3
            )
            registry.add(instance)

        # Test energy-based filtering
        strong_binders = [ii for ii in registry if ii.energy < -1.5]
        self.assertEqual(len(strong_binders), 2)  # energy -2 and -3

        # Test coordinate-based operations
        max_x = max(ii.absolute_coord[0] for ii in registry)
        self.assertEqual(max_x, 2.0)

        # Test name-based operations instead of direct object comparison
        interface_names = [ii.get_name() for ii in registry]
        self.assertIn("A_B_1", interface_names)
        self.assertIn("B_A_1", interface_names)
        self.assertIn("A_C_1", interface_names)


class TestRegistryEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions for all registries."""

    def test_empty_registry_operations(self):
        """Test operations on empty registries."""
        registries = [
            MoleculeTypeRegistry(),
            MoleculeInstanceRegistry(),
            InterfaceTypeRegistry(),
            InterfaceInstanceRegistry()
        ]

        for registry in registries:
            # Test iteration on empty registry
            items = list(registry)
            self.assertEqual(len(items), 0)

            # Test length
            self.assertEqual(len(registry), 0)

            # Test contains
            self.assertNotIn("anything", registry)

            # Test to_list
            result = registry.to_list()
            self.assertEqual(result, [])

            # Test clear (should not raise error)
            registry.clear()

    def test_registry_with_none_values(self):
        """Test registry behavior with None or invalid inputs."""
        _ = MoleculeTypeRegistry()

        # Test from_list with None
        registry_from_none = MoleculeTypeRegistry.from_list(None)
        self.assertEqual(len(registry_from_none), 0)

        # Test from_list with empty list
        registry_from_empty = MoleculeTypeRegistry.from_list([])
        self.assertEqual(len(registry_from_empty), 0)

    def test_large_registry_performance(self):
        """Test registry performance with large numbers of items."""
        registry = MoleculeTypeRegistry()

        # Add many items
        n_items = 1000
        for i in range(n_items):
            mol_type = MoleculeType(name=f"Mol_{i:04d}", radius_nm=1.0)
            registry.add(mol_type)

        # Test that operations still work efficiently
        self.assertEqual(len(registry), n_items)

        # Test membership
        self.assertIn("Mol_0500", registry)
        self.assertNotIn("Mol_9999", registry)

        # Test iteration (should complete without issues)
        count = sum(1 for _ in registry)
        self.assertEqual(count, n_items)

        # Test retrieval
        retrieved = registry.get("Mol_0750")
        self.assertEqual(retrieved.name, "Mol_0750")

    def test_registry_name_edge_cases(self):
        """Test registries with edge case names."""
        registry = MoleculeTypeRegistry()

        # Test with special characters
        special_names = [
            "Protein-A", "Protein_B", "Protein.C", "Protein@D",
            "123Protein", "", "   ", "Protein with spaces"
        ]

        for name in special_names:
            if name.strip():  # Skip empty/whitespace names for MoleculeType
                mol_type = MoleculeType(name=name, radius_nm=1.0)
                registry.add(mol_type)

                # Test that we can retrieve it
                retrieved = registry.get(name)
                self.assertEqual(retrieved.name, name)


class TestRegistryIntegration(unittest.TestCase):
    """Integration tests for registries working together."""

    def test_cross_registry_references(self):
        """Test that objects in different registries can reference each other."""
        mol_type_registry = MoleculeTypeRegistry()
        mol_instance_registry = MoleculeInstanceRegistry()

        # Create molecule type
        mol_type = MoleculeType(name="TestProtein", radius_nm=2.0)
        mol_type_registry.add(mol_type)

        # Create molecule instance referencing the type
        mol_instance = MoleculeInstance(
            name="TestProtein_001",
            norm=np.array([0.0, 0.0, 1.0]),
            ref1=np.array([1.0, 0.0, 0.0]),
            ref2=np.array([0.0, 0.0, 1.0]),
            com=np.array([0.0, 0.0, 0.0]),
            molecule_type=mol_type
        )
        mol_instance_registry.add(mol_instance)

        # Verify cross-references work
        retrieved_instance = mol_instance_registry.get("TestProtein_001")
        self.assertEqual(retrieved_instance.molecule_type, mol_type)
        self.assertEqual(retrieved_instance.molecule_type.name, "TestProtein")

    def test_registry_serialization_integration(self):
        """Test serialization of registries with cross-references."""
        mol_type_registry = MoleculeTypeRegistry()
        interface_type_registry = InterfaceTypeRegistry()

        # Create molecule types
        mol_type_a = MoleculeType(name="A", radius_nm=1.0)
        mol_type_b = MoleculeType(name="B", radius_nm=2.0)
        mol_type_registry.add(mol_type_a)
        mol_type_registry.add(mol_type_b)

        # Create interface types with references
        interface_ab = InterfaceType(
            this_mol_type_name="A",
            partner_mol_type_name="B",
            interface_index=1,
            absolute_coord=np.array([1.0, 0.0, 0.0]),
            local_coord=np.array([0.5, 0.0, 0.0]),
            this_mol_type=mol_type_a,
            partner_mol_type=mol_type_b
        )
        interface_type_registry.add(interface_ab)

        # Test serialization
        mol_data = mol_type_registry.to_list()
        _ = interface_type_registry.to_list()

        # Test deserialization
        new_mol_registry = MoleculeTypeRegistry.from_list(mol_data)
        # Note: Interface registry deserialization would need special handling
        # for cross-references, which is beyond basic serialization

        self.assertEqual(len(new_mol_registry), 2)
        self.assertIn("A", new_mol_registry)
        self.assertIn("B", new_mol_registry)


class TestRegistryDocumentation(unittest.TestCase):
    """Test that registries behave as documented in examples."""

    def test_documentation_example(self):
        """Test the example from the module docstring."""
        # Create and populate a molecule type registry
        mol_registry = MoleculeTypeRegistry()

        protein_a = MoleculeType(name="ProteinA", radius_nm=2.5)
        protein_b = MoleculeType(name="ProteinB", radius_nm=3.0)

        mol_registry.add(protein_a)
        mol_registry.add(protein_b)

        # Iterate through all molecule types
        names_and_radii = []
        for mol_type in mol_registry:
            names_and_radii.append((mol_type.name, mol_type.radius_nm))

        self.assertEqual(len(names_and_radii), 2)
        self.assertIn(("ProteinA", 2.5), names_and_radii)
        self.assertIn(("ProteinB", 3.0), names_and_radii)

        # Check membership and retrieve
        if "ProteinA" in mol_registry:
            retrieved = mol_registry.get("ProteinA")
            self.assertEqual(retrieved.name, "ProteinA")

        # Serialize for storage
        serialized = mol_registry.to_list()
        restored_registry = MoleculeTypeRegistry.from_list(serialized)

        # Verify restoration
        self.assertEqual(len(restored_registry), len(mol_registry))
        for original_mol in mol_registry:
            self.assertIn(original_mol.name, restored_registry)

    def test_iteration_use_cases(self):
        """Test common iteration use cases mentioned in documentation."""
        registry = MoleculeTypeRegistry()

        # Add diverse molecule types
        molecules = [
            MoleculeType(name="SmallProtein", radius_nm=1.0),
            MoleculeType(name="MediumProtein", radius_nm=2.5),
            MoleculeType(name="LargeProtein", radius_nm=4.0),
            MoleculeType(name="HugeProtein", radius_nm=6.0)
        ]

        for mol in molecules:
            registry.add(mol)

        # Use case 1: Find all large proteins
        large_proteins = [mol for mol in registry if mol.radius_nm > 3.0]
        self.assertEqual(len(large_proteins), 2)

        # Use case 2: Calculate average radius
        avg_radius = sum(mol.radius_nm for mol in registry) / len(registry)
        expected_avg = (1.0 + 2.5 + 4.0 + 6.0) / 4
        self.assertEqual(avg_radius, expected_avg)

        # Use case 3: Check if any protein meets criteria
        has_huge_protein = any(mol.radius_nm > 5.0 for mol in registry)
        self.assertTrue(has_huge_protein)

        # Use case 4: Get sorted list by size
        sorted_by_size = sorted(registry, key=lambda mol: mol.radius_nm)
        self.assertEqual(sorted_by_size[0].name, "SmallProtein")
        self.assertEqual(sorted_by_size[-1].name, "HugeProtein")


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)

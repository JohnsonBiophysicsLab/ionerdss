"""Tests for the pre-simulation over-assembly warning.

A design whose coarse-grained molecules leave declared interfaces unbound can
assemble past the deposited stoichiometry once more than one copy is supplied.
The capacity is not visible on an instance: its ``interfaces_neighbors_map``
holds only the contacts actually observed, so a free interface never appears as
``None`` and has to be derived from what the molecule *type* declares.
"""

import numpy as np

from ionerdss.model.components.instances import InterfaceInstance, MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import InterfaceType, MoleculeType
from ionerdss.model.pdb.structure_validation import (
    get_free_interface_capacity,
    get_free_interface_message,
)


def _make_instance(name: str, molecule_type: MoleculeType, com,
                   interface_count: int) -> MoleculeInstance:
    return MoleculeInstance(
        name=name,
        molecule_type=molecule_type,
        com=np.asarray(com, dtype=float),
        norm=np.array([0.0, 0.0, 1.0]),
        ref1=np.array([1.0, 0.0, 0.0]),
        ref2=np.array([0.0, 1.0, 0.0]),
        interfaces_neighbors_map={f"iface_{i}": None for i in range(interface_count)},
    )


def _free_interface_system(*, saturate: bool):
    """Two A copies bound through AA1. Type A also declares AA2, left free unless saturated."""
    system = System(workspace_path=".")
    type_a = MoleculeType(name="A")
    type_a.interfaces_neighbors_map = {"AA1": "A", "AA2": "A"}
    system.molecule_types.add(type_a)

    interface_type_1 = InterfaceType(
        this_mol_type_name="A", partner_mol_type_name="A", interface_index=1,
        absolute_coord=np.zeros(3), local_coord=np.zeros(3),
    )
    interface_type_2 = InterfaceType(
        this_mol_type_name="A", partner_mol_type_name="A", interface_index=2,
        absolute_coord=np.zeros(3), local_coord=np.zeros(3),
    )

    first = _make_instance("a0", type_a, [0.0, 0.0, 0.0], interface_count=0)
    second = _make_instance("a1", type_a, [2.0, 0.0, 0.0], interface_count=0)

    def _interface(owner, owner_name, partner_name, interface_type):
        return InterfaceInstance(
            absolute_coord=np.zeros(3), this_mol=owner, this_mol_name=owner_name,
            partner_mol_name=partner_name, interface_index=interface_type.interface_index,
            interface_type=interface_type,
        )

    first.interfaces_neighbors_map = {
        _interface(first, "a0", "a1", interface_type_1): second
    }
    second.interfaces_neighbors_map = {
        _interface(second, "a1", "a0", interface_type_1): first
    }
    if saturate:
        first.interfaces_neighbors_map[_interface(first, "a0", "a1", interface_type_2)] = second
        second.interfaces_neighbors_map[_interface(second, "a1", "a0", interface_type_2)] = first

    system.molecule_instances.add(first)
    system.molecule_instances.add(second)
    return system


def test_free_interface_capacity_counts_unbound_declared_interfaces():
    # AA1 is used by both copies; AA2 is declared by the type but never realized.
    capacity = get_free_interface_capacity(_free_interface_system(saturate=False))
    assert capacity == {"A": {"AA2": 2}}


def test_free_interface_capacity_empty_when_assembly_is_saturated():
    # Every declared interface is bound, so more copies cannot extend the assembly.
    assert get_free_interface_capacity(_free_interface_system(saturate=True)) == {}
    assert get_free_interface_message(
        _free_interface_system(saturate=True), prefix="Validation preflight warning"
    ) is None


def test_free_interface_message_reports_over_assembly_risk():
    message = get_free_interface_message(
        _free_interface_system(saturate=False), prefix="Validation preflight warning"
    )
    assert message is not None
    assert "leaves 2 interface slots unbound" in message
    assert "over-assembly" in message
    assert "A (AA2 on 2 copies)" in message


def test_free_interface_capacity_ignores_instances_without_a_type():
    system = _free_interface_system(saturate=True)
    orphan = _make_instance("x0", None, [9.0, 9.0, 9.0], interface_count=0)
    system.molecule_instances.add(orphan)
    assert get_free_interface_capacity(system) == {}

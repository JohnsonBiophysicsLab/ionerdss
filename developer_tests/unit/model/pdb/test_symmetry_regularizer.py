"""Tests for automatic point-group detection and Cn regularization."""

import numpy as np
import pytest

from ionerdss.model.components.instances import InterfaceInstance, MoleculeInstance
from ionerdss.model.components.system import System
from ionerdss.model.components.types import MoleculeType
from ionerdss.model.pdb.symmetry_regularizer import SymmetryRegularizer


def _frame(instance):
    e1 = np.asarray(instance.ref1, float)
    e1 = e1 / np.linalg.norm(e1)
    raw = np.asarray(instance.ref2, float)
    e2 = raw - (raw @ e1) * e1
    e2 = e2 / np.linalg.norm(e2)
    return np.column_stack([e1, e2, np.cross(e1, e2)])


def _ordered(system):
    """Instances in ring order a0, a1, ... (the registry is not subscriptable)."""
    return sorted(system.molecule_instances, key=lambda i: int(i.name[1:]))


def _twist(a, b):
    rot = _frame(b) @ _frame(a).T
    return float(np.degrees(np.arccos(np.clip((np.trace(rot) - 1) / 2, -1, 1))))


def _ring_system(n, *, radius=10.0, angle_jitter=None, share_one_orientation=False,
                 z_jitter=None):
    """Cyclic homomer of n subunits bonded head-to-tail around the z axis."""
    system = System(workspace_path=".")
    mol_type = MoleculeType(name="A")
    mol_type.interfaces_neighbors_map = {"AA1f": "A", "AA1b": "A"}
    system.molecule_types.add(mol_type)

    step = 2.0 * np.pi / n
    instances = []
    for k in range(n):
        theta = k * step + (angle_jitter[k] if angle_jitter is not None else 0.0)
        z = z_jitter[k] if z_jitter is not None else 0.0
        com = np.array([radius * np.cos(theta), radius * np.sin(theta), z])
        if share_one_orientation:
            ref1, ref2 = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
        else:
            ref1 = np.array([np.cos(theta), np.sin(theta), 0.0])
            ref2 = np.array([-np.sin(theta), np.cos(theta), 0.0])
        inst = MoleculeInstance(name=f"a{k}", molecule_type=mol_type, com=com,
                                norm=np.array([0.0, 0.0, 1.0]), ref1=ref1, ref2=ref2,
                                interfaces_neighbors_map={})
        instances.append(inst)

    for k, inst in enumerate(instances):
        for tag, partner in (("AA1f", instances[(k + 1) % n]), ("AA1b", instances[(k - 1) % n])):
            interface = InterfaceInstance(
                absolute_coord=np.asarray(inst.com, float) + np.array([0.1, 0.0, 0.0]),
                this_mol=inst, this_mol_name=inst.name, partner_mol_name=partner.name,
                interface_index=0,
            )
            inst.interfaces_neighbors_map[interface] = partner
        system.molecule_instances.add(inst)
    return system


def test_detects_perfect_cyclic_ring():
    detection = SymmetryRegularizer(_ring_system(6)).detect()
    assert detection.group == "C6"
    assert detection.order == 6
    assert detection.regularizable


def test_regularizes_jittered_ring_to_exact_fold_angle():
    rng = np.random.default_rng(0)
    jitter = rng.normal(scale=0.03, size=5)
    system = _ring_system(5, angle_jitter=jitter)
    ring = _ordered(system)
    before = [_twist(ring[k], ring[(k + 1) % 5]) for k in range(5)]
    assert max(abs(t - 72.0) for t in before) > 0.5   # genuinely off before

    assert SymmetryRegularizer(system).apply("auto") is True

    ring = _ordered(system)
    after = [_twist(ring[k], ring[(k + 1) % 5]) for k in range(5)]
    assert max(abs(t - 72.0) for t in after) < 1e-6   # exactly 5-fold after


def test_synthesises_orientations_when_alignment_left_them_identical():
    # The dominant real-world failure: every copy keeps the representative's frame, so
    # the generator transform is a pure translation and the ring can never close.
    system = _ring_system(4, share_one_orientation=True)
    ring = _ordered(system)
    before = [_twist(ring[k], ring[(k + 1) % 4]) for k in range(4)]
    assert max(before) < 1e-9                          # all copies share one frame

    assert SymmetryRegularizer(system).apply("auto") is True

    ring = _ordered(system)
    after = [_twist(ring[k], ring[(k + 1) % 4]) for k in range(4)]
    assert max(abs(t - 90.0) for t in after) < 1e-6


def test_interfaces_follow_their_molecule():
    system = _ring_system(4, angle_jitter=np.array([0.0, 0.05, -0.04, 0.03]))
    offsets_before = {
        inst.name: [np.linalg.norm(np.asarray(i.absolute_coord, float) - np.asarray(inst.com, float))
                    for i in inst.interfaces_neighbors_map]
        for inst in system.molecule_instances
    }
    SymmetryRegularizer(system).apply("auto")
    for inst in system.molecule_instances:
        after = [np.linalg.norm(np.asarray(i.absolute_coord, float) - np.asarray(inst.com, float))
                 for i in inst.interfaces_neighbors_map]
        assert np.allclose(sorted(after), sorted(offsets_before[inst.name]), atol=1e-9)


def test_refuses_when_a_subunit_would_move_too_far():
    # A ring so distorted that snapping it would relocate subunits wholesale.
    system = _ring_system(4, z_jitter=np.array([0.0, 40.0, 0.0, -40.0]))
    regularizer = SymmetryRegularizer(system, com_shift_cap=1.0)
    coms_before = [np.array(inst.com, float) for inst in _ordered(system)]
    regularizer.apply("auto")
    for before, inst in zip(coms_before, _ordered(system)):
        assert np.allclose(before, np.asarray(inst.com, float))


def test_leaves_non_cyclic_assemblies_alone():
    system = _ring_system(4)
    # Break the cycle: one subunit loses a partner, leaving an open chain.
    first = _ordered(system)[0]
    key = next(iter(first.interfaces_neighbors_map))
    first.interfaces_neighbors_map[key] = None

    detection = SymmetryRegularizer(system).detect()
    assert not detection.regularizable

    coms_before = [np.array(inst.com, float) for inst in _ordered(system)]
    assert SymmetryRegularizer(system).apply("auto") is False
    for before, inst in zip(coms_before, _ordered(system)):
        assert np.allclose(before, np.asarray(inst.com, float))


def test_off_mode_is_a_no_op():
    system = _ring_system(5, angle_jitter=np.array([0.0, 0.1, -0.1, 0.05, -0.05]))
    coms_before = [np.array(inst.com, float) for inst in _ordered(system)]
    assert SymmetryRegularizer(system).apply("off") is False
    for before, inst in zip(coms_before, _ordered(system)):
        assert np.allclose(before, np.asarray(inst.com, float))


def test_dihedral_symbol_still_yields_a_usable_cyclic_order():
    # Dn contains Cn, and the n-fold rotation is the only thing regularization
    # applies. Rejecting 'D6' here would push dihedral assemblies out of the ring
    # path their cyclic subgroup qualifies them for.
    order_of = SymmetryRegularizer._cyclic_order_from_symbol
    assert order_of("C6") == 6
    assert order_of("D6") == 6
    assert order_of("T") is None
    assert order_of("Cinfv") is None
    assert order_of(None) is None


def test_family_label_is_not_mistaken_for_a_verified_ring():
    # 'C-family' must not satisfy regularizable on its leading letter, and must be
    # visibly distinct from a verified 'C6' wherever the two are tabulated together.
    from ionerdss.model.pdb.symmetry_regularizer import SymmetryDetection

    family = SymmetryDetection(group=SymmetryRegularizer._group_family("C4"), order=8)
    assert family.group == "C-family"
    assert not family.regularizable

    verified = SymmetryDetection(group="C6", order=6, axis=np.array([0.0, 0.0, 1.0]))
    assert verified.regularizable


def test_dimer_is_too_small_to_regularize():
    detection = SymmetryRegularizer(_ring_system(2)).detect()
    assert not detection.regularizable
    assert "nothing to regularize" in detection.reason


# --------------------------------------------------------------------------
# Dihedral resolution: PointGroup sees points, an assembly has oriented bodies
# --------------------------------------------------------------------------


def _stacked_ring_system(n, *, chiral):
    """Two rings stacked along z, bonded within each ring.

    With chiral=False the lower ring is the mirror partner of the upper one, so a
    C2 in the plane maps the assembly onto itself and the group really is Dn.
    With chiral=True every subunit carries the same handed frame, which breaks
    those axes and leaves only Cn.
    """
    system = System(workspace_path=".")
    mol_type = MoleculeType(name="A")
    mol_type.interfaces_neighbors_map = {"AA1f": "A", "AA1b": "A"}
    system.molecule_types.add(mol_type)

    instances = []
    for tier, z in enumerate((6.0, -6.0)):
        for k in range(n):
            theta = 2.0 * np.pi * k / n
            com = np.array([10 * np.cos(theta), 10 * np.sin(theta), z])
            ref1 = np.array([np.cos(theta), np.sin(theta), 0.0])
            if chiral or tier == 0:
                ref2 = np.array([0.0, 0.0, 1.0])
            else:
                ref2 = np.array([0.0, 0.0, -1.0])   # lower tier flipped
            inst = MoleculeInstance(
                name=f"a{tier * n + k}", molecule_type=mol_type, com=com,
                norm=np.array([0.0, 0.0, 1.0]), ref1=ref1, ref2=ref2,
                interfaces_neighbors_map={})
            instances.append(inst)
            system.molecule_instances.add(inst)
    return system


def test_point_group_reports_dihedral_when_frames_agree():
    system = _stacked_ring_system(4, chiral=False)
    regularizer = SymmetryRegularizer(system)
    instances = _ordered(system)
    assert regularizer._point_group_symbol(instances) == "D4"


def test_dihedral_is_downgraded_when_frames_are_all_one_hand():
    # The centres of mass are identical to the case above, so PointGroup alone
    # still says D4; only the frames distinguish them.
    system = _stacked_ring_system(4, chiral=True)
    regularizer = SymmetryRegularizer(system)
    instances = _ordered(system)
    assert regularizer._point_group_symbol(instances) == "C4"


def test_flat_protein_ring_is_not_reported_as_dihedral():
    # The dominant real case: n chiral subunits on a circle. The point set is Dn,
    # the assembly is Cn, and reporting Dn would relabel most cyclic rings.
    system = _ring_system(6)
    assert SymmetryRegularizer(system)._point_group_symbol(_ordered(system)) == "C6"


def test_flat_ring_still_detects_as_a_cyclic_ring():
    detection = SymmetryRegularizer(_ring_system(6)).detect()
    assert detection.group == "C6"
    assert detection.point_group_symbol == "C6"


# --------------------------------------------------------------------------
# Orientation alignment: pairing Cα atoms by residue number rather than
# positionally, so copies that resolve different loops still get a real frame.
# --------------------------------------------------------------------------

from ionerdss.model.pdb.system_builder import _paired_ca_coordinates


def _chain(residue_ids, names=None, coords=None):
    names = names or {}
    residues = []
    for k, rid in enumerate(residue_ids):
        residues.append({
            "id": rid,
            "name": names.get(rid, "ALA"),
            "ca_coord": np.array(coords[k] if coords is not None else [float(rid), 0.0, 0.0]),
        })
    return {"residues": residues,
            "ca_coords": np.array([r["ca_coord"] for r in residues])}


def test_pairs_chains_of_unequal_length_on_shared_residues():
    # The representative resolves a loop (residues 5-6) that the copy does not.
    rep = _chain([1, 2, 3, 4, 5, 6, 7])
    curr = _chain([1, 2, 3, 4, 7])
    P, Q = _paired_ca_coordinates(rep, curr)
    assert P is not None and len(P) == 5
    assert np.allclose(P, Q)          # same residues, same synthetic coordinates


def test_pairing_ignores_residues_missing_from_either_chain():
    rep = _chain([10, 11, 12, 13])
    curr = _chain([11, 12, 13, 99])
    P, _Q = _paired_ca_coordinates(rep, curr)
    assert len(P) == 3                # 11, 12, 13 only


def test_pairing_rejects_mismatched_residue_names():
    rep = _chain([1, 2, 3, 4], names={3: "GLY"})
    curr = _chain([1, 2, 3, 4], names={3: "TRP"})
    P, _Q = _paired_ca_coordinates(rep, curr)
    assert len(P) == 3                # residue 3 dropped as a numbering conflict


def test_pairing_drops_ambiguous_duplicate_residue_numbers():
    # Insertion codes are not preserved in chain_data, so a repeated number is
    # ambiguous and must be discarded rather than guessed.
    rep = _chain([1, 2, 2, 3, 4])
    curr = _chain([1, 2, 3, 4])
    P, _Q = _paired_ca_coordinates(rep, curr)
    assert len(P) == 3                # 1, 3, 4


def test_pairing_reports_failure_when_overlap_is_too_small():
    rep = _chain([1, 2, 3])
    curr = _chain([7, 8, 9])
    P, Q = _paired_ca_coordinates(rep, curr)
    assert P is None and Q is None


def test_pairing_recovers_a_real_rotation_between_copies():
    from ionerdss.model.pdb.system_builder import _kabsch_rotation
    from ionerdss.utils.rotations import rotation_matrix

    rng = np.random.default_rng(3)
    base = rng.normal(size=(9, 3))
    turn = rotation_matrix(np.array([0.0, 0.0, 1.0]), np.deg2rad(120.0))

    rep = _chain(list(range(1, 10)), coords=base)
    # The copy is rotated by 120 degrees AND is missing two residues.
    rotated = base @ turn.T
    keep = [0, 1, 2, 4, 5, 6, 8]
    curr = _chain([k + 1 for k in keep], coords=rotated[keep])

    P, Q = _paired_ca_coordinates(rep, curr)
    assert len(P) == len(keep)
    assert np.allclose(_kabsch_rotation(P, Q), turn, atol=1e-8)

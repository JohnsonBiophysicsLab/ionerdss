"""
ionerdss.model.pdb.symmetry_regularizer

Automatic point-group detection and geometric regularization for coarse-grained
assemblies.

WHY THIS EXISTS
---------------
A cyclic homomer built from head-to-tail self-binding interfaces (``AA1f`` binds
``AA1b``) is a polymer subunit that happens to close into a ring. The ring closes in
simulation only if the generator transform ``T`` taking subunit *k* to subunit *k+1*
satisfies ``T**n == I``. Two things break that:

* the per-subunit orientation was never resolved, so every copy shares one frame and
  ``T`` is a pure translation -- ``T**n`` is then never identity for any *n*, and the
  design is a straight polymer that grows without bound; or
* the deposited geometry is only approximately n-fold, so the closure error
  accumulates over *n* bonds until the last one falls outside NERDSS's binding
  tolerance.

Either way the chain elongates instead of closing, which surfaces as over-assembly
once more than one copy of the deposited stoichiometry is supplied.

This module snaps such an assembly onto its exact point group. Orientations are
*synthesised from the group element* rather than recovered by structural alignment,
which is what makes it robust to the first failure above.

SCOPE
-----
Cyclic (``Cn``) assemblies are detected and regularized. Helical/open assemblies are
detected and deliberately left alone -- a filament such as actin genuinely extends
past the deposited asymmetric unit, so forcing closure there would be wrong. Dihedral
(``Dn``) assemblies are detected, and are regularized only when their subunits form a
single n-fold orbit, in which case the operation applied is the ``Cn`` rotation about
the principal axis; the perpendicular C2 axes are never touched. Cubic groups
(``T``/``O``/``I``) are reported but not regularized. Anything unrecognised is left
untouched.

A caution specific to ``Dn``: PointGroup classifies centres of mass, and n points
evenly spaced on a circle carry the perpendicular C2 axes of ``Dn`` whether or not the
bodies sitting on them do. Protein subunits are chiral, so a flat ``Cn`` ring of them
has dihedral centres of mass and a cyclic assembly. Every ``Dn`` symbol is therefore
re-checked against the subunit frames before it is reported, and downgraded to ``Cn``
when the C2 operation fails to map those frames onto each other.

Detection reuses :class:`ionerdss.model.graph_based.symmetry.pointgroup.PointGroup`
and the rotation helper in :mod:`ionerdss.utils.rotations` rather than reimplementing
either.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ionerdss.model.components.system import System
from ionerdss.model.components.instances import MoleculeInstance
from ionerdss.utils.rotations import rotation_matrix
from .file_manager import WorkspaceManager

__all__ = ["SymmetryDetection", "SymmetryRegularizer"]


# A ring is accepted as n-fold only if rotating the centres of mass by 2*pi/n maps the
# set onto itself to within this fraction of the ring radius.
_RING_MATCH_RELATIVE_TOLERANCE = 0.15
# Refuse to move any subunit further than this, in the system's own length units.
_DEFAULT_COM_SHIFT_CAP = 6.0
# Step of the sweep for a C2 axis perpendicular to the principal axis, in degrees.
_C2_SCAN_STEP_DEG = 1.0
# How far a subunit frame may be turned from its image and still count as matched.
_FRAME_MATCH_TOLERANCE_DEG = 20.0


@dataclass
class SymmetryDetection:
    """What the detector concluded about an assembly."""

    group: str = "unknown"                      # 'C<n>' verified, '<X>-family', 'helical', 'none'
    order: int = 1                              # n for Cn
    axis: Optional[np.ndarray] = None           # unit vector, cyclic groups only
    centre: Optional[np.ndarray] = None
    ring: Tuple[str, ...] = ()                  # instance names in angular order
    point_group_symbol: Optional[str] = None    # what PointGroup said, for cross-check
    reason: str = ""                            # why regularization is or is not possible

    @property
    def regularizable(self) -> bool:
        # Only a verified 'C<n>' qualifies; 'C-family' means the ring was never
        # confirmed, and must not match on the leading letter alone.
        verified_cyclic = self.group.startswith("C") and self.group[1:].isdigit()
        return verified_cyclic and self.order >= 3 and self.axis is not None


class SymmetryRegularizer:
    """Detect an assembly's point group and snap its geometry onto it.

    Attributes:
        system: ionerdss System to process, modified in place.
        workspace_manager: Workspace manager, used only for logging.
        com_shift_cap: Refuse to regularize if any subunit would move further than this.
    """

    def __init__(
        self,
        system: System,
        workspace_manager: Optional[WorkspaceManager] = None,
        *,
        com_shift_cap: float = _DEFAULT_COM_SHIFT_CAP,
        fold_tolerance: float = _RING_MATCH_RELATIVE_TOLERANCE,
    ):
        self.system = system
        self.workspace_manager = workspace_manager
        self.com_shift_cap = float(com_shift_cap)
        self.fold_tolerance = float(fold_tolerance)

    # ------------------------------------------------------------------ logging
    def _log(self, level: str, message: str, *args) -> None:
        if self.workspace_manager is not None:
            getattr(self.workspace_manager.logger, level)(message, *args)

    # ------------------------------------------------------------------ public
    def apply(self, mode: str = "auto") -> bool:
        """Detect and regularize. Returns True if the geometry was changed."""
        if mode in ("off", None):
            return False

        detection = self.detect()
        self.system.symmetry_detection = detection  # recorded either way

        if mode == "auto" and not detection.regularizable:
            self._log("info", "Symmetry regularization skipped: %s", detection.reason)
            return False
        if mode not in ("auto", "cn"):
            self._log("warning", "Unsupported geometric_regularization mode %r; skipping", mode)
            return False
        if not detection.regularizable:
            self._log("info", "Symmetry regularization skipped: %s", detection.reason)
            return False

        return self._regularize_cyclic(detection)

    def detect(self) -> SymmetryDetection:
        """Classify the assembly from its contact graph, then confirm geometrically."""
        instances = list(self.system.molecule_instances)
        n = len(instances)
        if n < 3:
            return SymmetryDetection(group="none", reason=f"only {n} chain(s); nothing to regularize")

        adjacency = self._contact_graph(instances)
        degrees = {name: len(partners) for name, partners in adjacency.items()}
        bonds = sum(degrees.values()) // 2

        # Open head-to-tail chain: a filament, which is meant to extend.
        if bonds == n - 1 and max(degrees.values(), default=0) <= 2:
            return SymmetryDetection(
                group="helical", order=n,
                reason="open head-to-tail chain (filament); extension is expected, not corrected")

        symbol = self._point_group_symbol(instances)

        # A clean degree-2 cycle is the easy case. Otherwise fall back on the
        # graph_based point-group detector: a Cn assembly whose subunits also touch
        # non-neighbours has a denser contact graph but is still an n-fold orbit, and
        # the geometric fold test below is what decides either way.
        ring = self._cycle_order(instances, adjacency) if all(d == 2 for d in degrees.values()) else None
        if ring is None:
            order = self._cyclic_order_from_symbol(symbol)
            if order is None or order < 3 or n % order != 0:
                return SymmetryDetection(
                    group=self._group_family(symbol), order=n, point_group_symbol=symbol,
                    reason=(f"contact graph is not a simple cycle (degrees "
                            f"{sorted(set(degrees.values()))}) and point group "
                            f"{symbol or 'unknown'} is not a usable Cn orbit"))
            ring = instances

        centre, axis = self._ring_frame([inst.com for inst in ring])
        if axis is None:
            return SymmetryDetection(group="none", order=n, reason="ring is degenerate; no axis")

        ordered, phase, radius = self._angular_order(ring, centre, axis)

        residual = self._fold_residual([inst.com for inst in ordered], centre, axis, len(ordered))
        if radius <= 0 or residual > self.fold_tolerance:
            return SymmetryDetection(
                group="none", order=n, point_group_symbol=symbol,
                reason=(f"cycle of {n} is not {n}-fold symmetric "
                        f"(residual {residual:.3f} > {self.fold_tolerance})"))

        return SymmetryDetection(
            group=f"C{len(ordered)}", order=len(ordered), axis=axis, centre=centre,
            ring=tuple(inst.name for inst in ordered), point_group_symbol=symbol,
            reason=f"cyclic C{len(ordered)} ring, fold residual {residual:.3f}")

    # ------------------------------------------------------------- regularize
    def _regularize_cyclic(self, detection: SymmetryDetection) -> bool:
        """Place an n-fold ring on exact Cn geometry, orientations included."""
        by_name = {inst.name: inst for inst in self.system.molecule_instances}
        ring = [by_name[name] for name in detection.ring]
        n = len(ring)
        centre = np.asarray(detection.centre, float)
        axis = np.asarray(detection.axis, float)

        basis_u, basis_v = self._plane_basis(axis)
        angles, radii, heights = [], [], []
        for inst in ring:
            offset = np.asarray(inst.com, float) - centre
            height = float(offset @ axis)
            planar = offset - height * axis
            angles.append(float(np.arctan2(planar @ basis_v, planar @ basis_u)))
            radii.append(float(np.linalg.norm(planar)))
            heights.append(height)

        step = 2.0 * np.pi / n
        # Circular mean of the per-subunit phase residual, so no single subunit is
        # treated as ground truth and the whole ring moves as little as possible.
        residual_phases = np.array([angles[k] - k * step for k in range(n)])
        phase = float(np.arctan2(np.sin(residual_phases).mean(), np.cos(residual_phases).mean()))
        radius = float(np.mean(radii))
        height = float(np.mean(heights))

        # The seed keeps its own (real) orientation; every other frame is generated
        # from it by the group element, never by structural alignment.
        seed = int(np.argmin(np.abs(self._wrap(residual_phases - phase))))
        seed_frame_ok = self._has_frame(ring[seed])

        proposed: List[Tuple[MoleculeInstance, np.ndarray, np.ndarray]] = []
        for k, inst in enumerate(ring):
            theta = phase + k * step
            new_com = centre + radius * (np.cos(theta) * basis_u + np.sin(theta) * basis_v) \
                + height * axis
            delta = rotation_matrix(axis, (k - seed) * step)
            proposed.append((inst, new_com, delta))

        shifts = [float(np.linalg.norm(new_com - np.asarray(inst.com, float)))
                  for inst, new_com, _ in proposed]
        worst = max(shifts) if shifts else 0.0
        if worst > self.com_shift_cap:
            self._log("warning",
                      "Symmetry regularization refused for %s: subunit would move %.2f "
                      "(cap %.2f)", detection.group, worst, self.com_shift_cap)
            return False

        seed_rotation = self._frame(ring[seed]) if seed_frame_ok else None
        for k, (inst, new_com, delta) in enumerate(proposed):
            old_com = np.asarray(inst.com, float)
            if seed_rotation is not None:
                new_frame = delta @ seed_rotation
                rotation = new_frame @ self._frame(inst).T if self._has_frame(inst) else delta
                inst.ref1 = new_frame[:, 0]
                inst.ref2 = new_frame[:, 1]
            else:
                rotation = delta
            inst.com = new_com
            if getattr(inst, "norm", None) is not None:
                inst.norm = rotation @ np.asarray(inst.norm, float)
            self._move_interfaces(inst, old_com, new_com, rotation)

        self._log("info",
                  "Regularized %s ring onto exact geometry (radius %.3f, worst COM shift "
                  "%.3f, orientations synthesised from the group element)",
                  detection.group, radius, worst)
        return True

    # ------------------------------------------------------------------ helpers
    @staticmethod
    def _wrap(angles: np.ndarray) -> np.ndarray:
        return (np.asarray(angles) + np.pi) % (2 * np.pi) - np.pi

    @staticmethod
    def _contact_graph(instances: Sequence[MoleculeInstance]) -> Dict[str, List[str]]:
        adjacency: Dict[str, List[str]] = {inst.name: [] for inst in instances}
        for inst in instances:
            for _interface, partner in inst.interfaces_neighbors_map.items():
                if partner is not None and partner.name in adjacency:
                    if partner.name not in adjacency[inst.name]:
                        adjacency[inst.name].append(partner.name)
        return adjacency

    @staticmethod
    def _cycle_order(instances: Sequence[MoleculeInstance],
                     adjacency: Dict[str, List[str]]) -> Optional[List[MoleculeInstance]]:
        by_name = {inst.name: inst for inst in instances}
        start = instances[0].name
        order, previous, current = [start], None, start
        while True:
            candidates = [x for x in adjacency[current] if x != previous]
            if not candidates:
                return None
            nxt = candidates[0]
            if nxt == start:
                break
            if nxt in order:
                return None
            order.append(nxt)
            previous, current = current, nxt
        return [by_name[name] for name in order] if len(order) == len(instances) else None

    @staticmethod
    def _ring_frame(coms: Sequence[np.ndarray]) -> Tuple[np.ndarray, Optional[np.ndarray]]:
        points = np.asarray([np.asarray(c, float) for c in coms])
        centre = points.mean(axis=0)
        centred = points - centre
        if np.allclose(centred, 0):
            return centre, None
        # Least-variance direction of the COM cloud is the ring normal.
        _u, _s, vt = np.linalg.svd(centred, full_matrices=False)
        axis = vt[-1]
        norm = np.linalg.norm(axis)
        return centre, (axis / norm if norm > 1e-9 else None)

    @staticmethod
    def _plane_basis(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        seed = np.array([1.0, 0.0, 0.0])
        if abs(float(seed @ axis)) > 0.9:
            seed = np.array([0.0, 1.0, 0.0])
        u = seed - (seed @ axis) * axis
        u /= np.linalg.norm(u)
        return u, np.cross(axis, u)

    def _angular_order(self, ring: Sequence[MoleculeInstance], centre: np.ndarray,
                       axis: np.ndarray) -> Tuple[List[MoleculeInstance], float, float]:
        u, v = self._plane_basis(axis)
        entries = []
        for inst in ring:
            offset = np.asarray(inst.com, float) - centre
            planar = offset - (offset @ axis) * axis
            entries.append((float(np.arctan2(planar @ v, planar @ u)),
                            float(np.linalg.norm(planar)), inst))
        entries.sort(key=lambda e: e[0])
        radius = float(np.mean([e[1] for e in entries])) if entries else 0.0
        return [e[2] for e in entries], (entries[0][0] if entries else 0.0), radius

    def _fold_residual(self, coms: Sequence[np.ndarray], centre: np.ndarray,
                       axis: np.ndarray, n: int) -> float:
        """How far a 2*pi/n rotation is from mapping the COM set onto itself."""
        points = np.asarray([np.asarray(c, float) for c in coms]) - centre
        rotated = points @ rotation_matrix(axis, 2.0 * np.pi / n).T
        radius = float(np.linalg.norm(points, axis=1).mean())
        if radius <= 0:
            return np.inf
        worst = 0.0
        for point in rotated:
            worst = max(worst, float(np.linalg.norm(points - point, axis=1).min()))
        return worst / radius

    def _point_group_symbol(self, instances: Sequence[MoleculeInstance]) -> Optional[str]:
        """Ask the graph_based PointGroup implementation what it sees."""
        try:
            from ionerdss.model.graph_based.symmetry.pointgroup import PointGroup

            coms = np.asarray([np.asarray(inst.com, float) for inst in instances])
            symbols = [inst.molecule_type.name if inst.molecule_type else "X"
                       for inst in instances]
            symbol = PointGroup(positions=coms - coms.mean(axis=0),
                                symbols=symbols).get_point_group()
            return self._resolve_dihedral(symbol, instances)
        except Exception as exc:  # detection must never break the build
            self._log("debug", "PointGroup detection unavailable: %s", exc)
            return None

    def _resolve_dihedral(self, symbol: Optional[str],
                          instances: Sequence[MoleculeInstance]) -> Optional[str]:
        """Downgrade Dn to Cn unless the subunit frames back it up.

        PointGroup classifies a cloud of points, and n points evenly spaced on a
        circle carry the n perpendicular C2 axes of Dn whether or not the bodies
        sitting on them do. Protein subunits are chiral, so a flat Cn ring of
        them has Dn centres of mass and a Cn assembly -- taking the point-set
        answer at face value would relabel most cyclic rings as dihedral.

        A C2 perpendicular to the principal axis is real only if it maps subunit
        *frames* onto each other as well as their positions. Test that; if it
        does not hold, or there are no frames to test, report the cyclic
        subgroup that is certain rather than the dihedral group that is not.
        """
        if not symbol or symbol[0].upper() != "D" or not symbol[1:].isdigit():
            return symbol
        if self._perpendicular_c2_maps_frames(instances):
            return symbol
        self._log("debug", "Dihedral %s not confirmed by subunit frames; reporting C%s",
                  symbol, symbol[1:])
        return "C" + symbol[1:]

    def _perpendicular_c2_maps_frames(self, instances: Sequence[MoleculeInstance]) -> bool:
        """Is there a C2 axis, perpendicular to the principal axis, that maps the
        oriented assembly onto itself?"""
        if len(instances) < 2 or not all(self._has_frame(i) for i in instances):
            return False

        coms = np.asarray([np.asarray(i.com, float) for i in instances])
        _centre, axis = self._ring_frame(coms)
        if axis is None:
            return False

        local = coms - coms.mean(axis=0)
        radius = float(np.linalg.norm(local, axis=1).max())
        if radius <= 0:
            return False
        frames = [self._frame(i) for i in instances]
        names = [i.molecule_type.name if i.molecule_type else "X" for i in instances]

        base, _v = self._plane_basis(axis)
        # The n C2 axes of Dn are spaced by pi/n, so half a turn meets one.
        for angle in np.arange(0.0, np.pi, np.deg2rad(_C2_SCAN_STEP_DEG)):
            u = rotation_matrix(axis, float(angle)) @ base
            if self._rigid_op_is_a_symmetry(rotation_matrix(u, np.pi),
                                            local, frames, names, radius):
                return True
        return False

    def _rigid_op_is_a_symmetry(self, rotation: np.ndarray, local: np.ndarray,
                                frames: Sequence[np.ndarray], names: Sequence[str],
                                radius: float) -> bool:
        """Does `rotation` send every subunit onto a same-type subunit, frame included?"""
        position_tolerance = self.fold_tolerance * radius
        moved = local @ rotation.T
        for k, point in enumerate(moved):
            distances = np.linalg.norm(local - point, axis=1)
            turned = rotation @ frames[k]
            matched = False
            for j in np.argsort(distances):
                if distances[j] > position_tolerance:
                    break                       # sorted, so nothing closer remains
                if names[j] != names[k]:
                    continue
                cosine = (np.trace(turned.T @ frames[j]) - 1.0) / 2.0
                if np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))) <= _FRAME_MATCH_TOLERANCE_DEG:
                    matched = True
                    break
            if not matched:
                return False
        return True

    @staticmethod
    def _cyclic_order_from_symbol(symbol: Optional[str]) -> Optional[int]:
        """n from a Schoenflies symbol such as 'C6' or 'D6'; None if there is none.

        Dn counts as well as Cn. Its principal axis carries the same n-fold
        rotation, and that rotation is the only operation regularization ever
        applies -- the perpendicular C2 axes are what distinguish the two groups
        and are left alone. Rejecting Dn here would drop dihedral assemblies out
        of the ring path that their cyclic subgroup qualifies them for.
        """
        if not symbol or symbol[0].upper() not in ("C", "D"):
            return None
        tail = symbol[1:]
        return int(tail) if tail.isdigit() else None

    @staticmethod
    def _group_family(symbol: Optional[str]) -> str:
        """Coarse label for an assembly whose n-fold ring could not be verified.

        Deliberately not spelled 'Cn'/'Dn': those sit beside verified labels like
        'C6' in every table and figure, where they read as the family that
        contains C6 when they mean the opposite -- the point group looked
        cyclic but no ring passed the fold test.
        """
        if not symbol:
            return "unknown"
        head = symbol[0].upper()
        return {"C": "C-family", "D": "D-family", "T": "T-family",
                "O": "O-family", "I": "I-family"}.get(head, "unknown")

    @staticmethod
    def _has_frame(inst: MoleculeInstance) -> bool:
        ref1 = np.asarray(getattr(inst, "ref1", None), float) if getattr(inst, "ref1", None) is not None else None
        ref2 = np.asarray(getattr(inst, "ref2", None), float) if getattr(inst, "ref2", None) is not None else None
        if ref1 is None or ref2 is None:
            return False
        if np.linalg.norm(ref1) < 1e-9 or np.linalg.norm(ref2) < 1e-9:
            return False
        return np.linalg.norm(np.cross(ref1, ref2)) > 1e-9

    @staticmethod
    def _frame(inst: MoleculeInstance) -> np.ndarray:
        """Orthonormal frame from the instance's reference vectors."""
        e1 = np.asarray(inst.ref1, float)
        e1 = e1 / np.linalg.norm(e1)
        raw = np.asarray(inst.ref2, float)
        e2 = raw - (raw @ e1) * e1
        e2 = e2 / np.linalg.norm(e2)
        return np.column_stack([e1, e2, np.cross(e1, e2)])

    @staticmethod
    def _move_interfaces(inst: MoleculeInstance, old_com: np.ndarray,
                         new_com: np.ndarray, rotation: np.ndarray) -> None:
        """Carry the molecule's interfaces along with its rigid motion."""
        for interface in inst.interfaces_neighbors_map:
            coord = getattr(interface, "absolute_coord", None)
            if coord is None:
                continue
            interface.absolute_coord = new_com + rotation @ (np.asarray(coord, float) - old_com)

"""Measure how faithfully a NERDSS run reproduces the assembly its model was built from.

Everything here is read from the bonded-complex JSON snapshots NERDSS writes when
``bondedComplexWrite`` is set, the exported ``.mol`` and ``parms.inp``, and the
``outputs/systems/*_system.json`` the PDB pipeline writes next to them.

Three things are measured.

**Geometric closure.** NERDSS places a bond between two complexes exactly, but a bond
*inside* one complex -- closing a ring -- only forms when the two sites already sit
within ``bindRadSameCom`` (1.1) times the binding radius, and nothing moves to help it.
A model whose bond geometries disagree with each other can therefore never close its
rings: subunits keep free sites, and another subunit binds them in the wrong place.
The transforms are taken from NERDSS's own placement of isolated dimers in the run.

**Mis-assembly.** Whether each simulated complex is a piece of the designed assembly.
When the design is a linear polymer -- the deposited subunits lie along one axis and
every copy of an interface type joins subunits a fixed number of positions apart -- each
complex is embedded back into that lattice. A bond that contradicts the embedding, or
two subunits landing in one lattice position, is a branch or a doubled strand. Steric
overlap, two subunit centres closer than the clash distance, is checked for every design.

**Stretched bonds.** Bonds whose sites sit well outside the loop-closure window. NERDSS's
rigid-body placement cannot produce one: they appear when two associations between the
same pair of complexes fire in the same timestep and the second is carried out as a loop
closure with no distance check. They glue otherwise clean complexes into one reported
complex, so they are counted separately from mis-assembly rather than mixed into it.
"""

from __future__ import annotations

import json
import re
from collections import defaultdict, deque
from dataclasses import dataclass, field
from itertools import combinations, product
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np

# NERDSS bonds two sites inside one complex when they are within bindRadSameCom * sigma.
# The default is 1.1 and ioNERDSS does not export another value.
LOOP_CLOSURE_FACTOR = 1.1
# A little above it, so a bond that merely closed a strained loop is not called stretched.
STRETCHED_BOND_FACTOR = 1.15
# Enumerating orientations of self-binding bonds is exponential; beyond this many in one
# complex, the greedy assignment stands.
MAX_ENUMERATED_AMBIGUOUS_BONDS = 12

SitePair = Tuple[Tuple[str, str], Tuple[str, str]]


def quaternion_to_matrix(q: Sequence[float]) -> np.ndarray:
    """Rotation matrix for the (w, x, y, z) quaternion NERDSS writes per molecule.

    It turns the ``.mol`` template into the molecule's current orientation, so a site
    sits at ``com + R @ local_coord``.
    """
    w, x, y, z = (float(v) for v in q)
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def read_mol_sites(mol_path: Path) -> Dict[str, np.ndarray]:
    """Site label -> local coordinate (nm) from one exported ``.mol`` file."""
    sites: Dict[str, np.ndarray] = {}
    in_coords = False
    for line in Path(mol_path).read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("# Coordinates"):
            in_coords = True
            continue
        if in_coords:
            if stripped.startswith("#"):
                break
            parts = stripped.split()
            if len(parts) == 4 and parts[0] != "COM":
                sites[parts[0]] = np.array([float(v) for v in parts[1:]])
    return sites


def read_binding_radii(parms_path: Path) -> Dict[SitePair, float]:
    """Binding radius of each reaction, keyed by its unordered pair of (molecule, site)."""
    text = Path(parms_path).read_text()
    reactions = re.finditer(
        r"^\s*(\w+)\(([^)]+)\)\s*\+\s*(\w+)\(([^)]+)\).*?^\s*sigma\s*=\s*([\d.eE+-]+)",
        text,
        re.M | re.S,
    )
    radii: Dict[SitePair, float] = {}
    for match in reactions:
        mol1, site1, mol2, site2, sigma = match.groups()
        # A site list may carry required-free interfaces; the bound one comes first.
        key = _site_pair((mol1, site1.split(",")[0].strip()), (mol2, site2.split(",")[0].strip()))
        radii[key] = float(sigma)
    return radii


def _site_pair(a: Tuple[str, str], b: Tuple[str, str]) -> SitePair:
    return tuple(sorted([a, b]))  # type: ignore[return-value]


@dataclass
class Model:
    """The exported NERDSS model: site geometry per molecule type plus binding radii."""

    sites: Dict[str, Dict[str, np.ndarray]]
    radii: Dict[SitePair, float]

    @classmethod
    def from_nerdss_dir(cls, nerdss_dir: Path) -> "Model":
        nerdss_dir = Path(nerdss_dir)
        sites = {path.stem: read_mol_sites(path) for path in sorted(nerdss_dir.glob("*.mol"))}
        return cls(sites=sites, radii=read_binding_radii(nerdss_dir / "parms.inp"))

    def radius(self, mol1: str, site1: str, mol2: str, site2: str) -> Optional[float]:
        return self.radii.get(_site_pair((mol1, site1), (mol2, site2)))

    def site_gap(self, complex_json: dict, bond: dict) -> Optional[float]:
        """Site separation of one bond, in units of that reaction's binding radius."""
        i, j = bond["molindex1"], bond["molindex2"]
        mol_i, mol_j = complex_json["names"][i], complex_json["names"][j]
        radius = self.radius(mol_i, bond["site1"], mol_j, bond["site2"])
        if not radius:
            return None
        coords = complex_json["coords"]
        pos_i = np.asarray(coords[i], float) + quaternion_to_matrix(complex_json["rotations"][i]) @ self.sites[mol_i][bond["site1"]]
        pos_j = np.asarray(coords[j], float) + quaternion_to_matrix(complex_json["rotations"][j]) @ self.sites[mol_j][bond["site2"]]
        return float(np.linalg.norm(pos_i - pos_j)) / radius


@dataclass
class LatticeDesign:
    """Signed lattice step of each interface type in a linear-polymer design.

    ``step`` holds the sites whose two ends are distinguishable, so the direction is
    known. ``magnitude`` holds self-binding sites, where the design fixes how far apart
    the partners are but not which way round, and the embedding has to try both.
    """

    step: Dict[str, int] = field(default_factory=dict)
    magnitude: Dict[str, int] = field(default_factory=dict)

    def steps_for(self, site: str) -> List[int]:
        if site in self.step:
            return [self.step[site]]
        if site in self.magnitude:
            return [self.magnitude[site], -self.magnitude[site]]
        return []


def designed_lattice(system_json: Path) -> Optional[LatticeDesign]:
    """Derive the lattice of a linear-polymer design, or None when it is not one.

    Subunits are ranked along the deposited assembly's principal axis -- no dependence on
    how the chains happen to be named -- and every copy of an interface type has to join
    subunits a fixed number of ranks apart. The derived lattice is then checked by
    re-embedding the design itself, so anything that is not a linear polymer (a ring, a
    cage, a branched design) is rejected here rather than mismeasured later.
    """
    registries = json.loads(Path(system_json).read_text())["registries"]
    instances = {inst["name"]: inst for inst in registries["molecule_instances"]}
    if len(instances) < 3:
        return None

    coms = np.array([inst["com"] for inst in instances.values()], float)
    centred = coms - coms.mean(axis=0)
    axis = np.linalg.eigh(centred.T @ centred)[1][:, -1]
    order = sorted(instances, key=lambda name: float(np.dot(np.asarray(instances[name]["com"], float) - coms.mean(axis=0), axis)))
    rank = {name: i for i, name in enumerate(order)}

    interface_type = {inst["name"]: inst["type"] for inst in registries["interface_instances"]}
    observed: Dict[str, set] = defaultdict(set)
    edges: List[Tuple[str, str, str]] = []
    for name, inst in instances.items():
        for interface in inst["interfaces"]:
            partner = interface["partner_molecule"]
            if partner not in rank:
                continue
            site = interface_type[interface["interface_name"]].lower()
            observed[site].add(rank[partner] - rank[name])
            edges.append((name, partner, site))

    lattice = LatticeDesign()
    for site, steps in observed.items():
        if len(steps) == 1:
            lattice.step[site] = steps.pop()
        elif len(steps) == 2 and sum(steps) == 0:
            lattice.magnitude[site] = abs(max(steps))
        else:
            return None  # this interface type does not sit at a fixed spacing

    # Confirm the design itself embeds: every edge has to agree with the ranks it came from.
    for owner, partner, site in edges:
        if (rank[partner] - rank[owner]) not in lattice.steps_for(site):
            return None
    return lattice


def _embed(n_molecules: int, bonds: Sequence[dict], lattice: LatticeDesign) -> Tuple[Dict[int, int], int]:
    """Assign a lattice position to every molecule; return the positions and conflict count.

    Self-binding bonds have two possible directions. They are resolved greedily, and when
    that leaves the complex inconsistent the orientations are enumerated (up to
    ``MAX_ENUMERATED_AMBIGUOUS_BONDS``) and the best assignment kept.
    """
    fixed, ambiguous = [], []
    for bond in bonds:
        options = lattice.steps_for(bond["site1"])
        if len(options) == 1:
            fixed.append((bond["molindex1"], bond["molindex2"], options[0]))
        elif options:
            ambiguous.append((bond["molindex1"], bond["molindex2"], options[0]))
        else:
            return {}, len(bonds)  # a site the design never used

    def run(signs: Sequence[int]) -> Tuple[Dict[int, int], int]:
        adjacency = defaultdict(list)
        for i, j, step in fixed:
            adjacency[i].append((j, step))
            adjacency[j].append((i, -step))
        for (i, j, step), sign in zip(ambiguous, signs):
            adjacency[i].append((j, sign * step))
            adjacency[j].append((i, -sign * step))
        positions: Dict[int, int] = {}
        conflicts = 0
        for start in range(n_molecules):
            if start in positions:
                continue
            positions[start] = 0
            queue = deque([start])
            while queue:
                current = queue.popleft()
                for neighbour, step in adjacency[current]:
                    expected = positions[current] + step
                    if neighbour not in positions:
                        positions[neighbour] = expected
                        queue.append(neighbour)
                    elif positions[neighbour] != expected:
                        conflicts += 1
        return positions, conflicts // 2  # each bad bond is seen from both ends

    best = run([1] * len(ambiguous))
    if best[1] == 0 or len(ambiguous) > MAX_ENUMERATED_AMBIGUOUS_BONDS:
        return best
    for signs in product((1, -1), repeat=len(ambiguous)):
        candidate = run(signs)
        if candidate[1] < best[1]:
            best = candidate
        if best[1] == 0:
            break
    return best


def analyse_complex(complex_json: dict, model: Model, lattice: Optional[LatticeDesign], clash_nm: float) -> dict:
    """Classify one bonded complex from a snapshot."""
    coords = np.asarray(complex_json["coords"], float)
    n_molecules = len(coords)

    placed, stretched = [], 0
    for bond in complex_json["bonds"]:
        gap = model.site_gap(complex_json, bond)
        if gap is not None and gap > STRETCHED_BOND_FACTOR:
            stretched += 1
        else:
            placed.append(bond)

    self_binding = sum(1 for bond in placed if bond["site1"] == bond["site2"])
    min_separation = min(
        (float(np.linalg.norm(coords[a] - coords[b])) for a, b in combinations(range(n_molecules), 2)),
        default=float("inf"),
    )

    result = {
        "n_molecules": n_molecules,
        "n_bonds": len(complex_json["bonds"]),
        "stretched_bonds": stretched,
        "self_binding_bonds": self_binding,
        "min_com_separation_nm": min_separation,
        "overlapping": min_separation < clash_nm,
        "lattice_conflicts": 0,
        "lattice_collisions": 0,
        "lattice_bonds_formed": None,
        "lattice_bonds_possible": None,
    }

    # Stretched bonds are a simulator artefact, so judge the pieces they hold together
    # separately: each connected component of the remaining bonds is one real structure.
    for component in _components(n_molecules, [(b["molindex1"], b["molindex2"]) for b in placed]):
        members = set(component)
        bonds = [b for b in placed if b["molindex1"] in members]
        if lattice is None:
            continue
        positions, conflicts = _embed(n_molecules, bonds, lattice)
        occupied = [positions[m] for m in component if m in positions]
        result["lattice_conflicts"] += conflicts
        result["lattice_collisions"] += len(occupied) - len(set(occupied))
        if conflicts == 0 and len(set(occupied)) == len(occupied):
            taken = set(occupied)
            possible = sum(
                1 for position in taken for step in _positive_steps(lattice) if position + step in taken
            )
            result["lattice_bonds_formed"] = (result["lattice_bonds_formed"] or 0) + len(bonds)
            result["lattice_bonds_possible"] = (result["lattice_bonds_possible"] or 0) + possible

    result["mis_assembled"] = bool(
        result["overlapping"] or result["lattice_conflicts"] or result["lattice_collisions"]
    )
    return result


def _positive_steps(lattice: LatticeDesign) -> List[int]:
    steps = {abs(step) for step in lattice.step.values()} | set(lattice.magnitude.values())
    return sorted(steps)


def _components(n_molecules: int, edges: Sequence[Tuple[int, int]]) -> List[List[int]]:
    adjacency = defaultdict(list)
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    seen: set = set()
    components = []
    for start in range(n_molecules):
        if start in seen:
            continue
        seen.add(start)
        queue, component = deque([start]), []
        while queue:
            current = queue.popleft()
            component.append(current)
            for neighbour in adjacency[current]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    queue.append(neighbour)
        components.append(component)
    return components


def iter_snapshots(sim_dir: Path) -> Iterator[Tuple[int, List[dict]]]:
    """Yield (timestep, complexes) for every bonded-complex snapshot, in time order."""
    snapshot_dir = Path(sim_dir) / "DATA" / "COMPLEXES"
    for path in sorted(snapshot_dir.glob("*.json"), key=lambda p: int(p.stem)):
        yield int(path.stem), json.loads(path.read_text())


def analyse_run(sim_dir: Path, system_json: Path, clash_nm: float = 3.0) -> dict:
    """Summarise one NERDSS run: mis-assembly, stretched bonds and final-state quality."""
    sim_dir = Path(sim_dir)
    model = Model.from_nerdss_dir(sim_dir)
    lattice = designed_lattice(system_json)

    summary = {
        "lattice_available": lattice is not None,
        "snapshots": 0,
        "mis_assembled_snapshots": 0,
        "first_mis_assembled_step": None,
        "largest_mis_assembled": 0,
        "stretched_bond_snapshots": 0,
        "first_stretched_bond_step": None,
        "worst_min_com_separation_nm": float("inf"),
        "self_binding_bonds_final": 0,
        "largest_complex_final": 1,
        "complex_sizes_final": [],
        "lattice_bond_completeness_final": None,
    }

    for step, complexes in iter_snapshots(sim_dir):
        results = [analyse_complex(c, model, lattice, clash_nm) for c in complexes]
        summary["snapshots"] += 1
        if any(r["mis_assembled"] for r in results):
            summary["mis_assembled_snapshots"] += 1
            summary["first_mis_assembled_step"] = summary["first_mis_assembled_step"] or step
            summary["largest_mis_assembled"] = max(
                summary["largest_mis_assembled"],
                max(r["n_molecules"] for r in results if r["mis_assembled"]),
            )
        if any(r["stretched_bonds"] for r in results):
            summary["stretched_bond_snapshots"] += 1
            summary["first_stretched_bond_step"] = summary["first_stretched_bond_step"] or step
        summary["worst_min_com_separation_nm"] = min(
            summary["worst_min_com_separation_nm"],
            min((r["min_com_separation_nm"] for r in results), default=float("inf")),
        )
        # Overwritten each snapshot, so the last one survives as the final state.
        summary["self_binding_bonds_final"] = sum(r["self_binding_bonds"] for r in results)
        summary["largest_complex_final"] = max((r["n_molecules"] for r in results), default=1)
        summary["complex_sizes_final"] = sorted((r["n_molecules"] for r in results), reverse=True)
        formed = sum(r["lattice_bonds_formed"] or 0 for r in results)
        possible = sum(r["lattice_bonds_possible"] or 0 for r in results)
        summary["lattice_bond_completeness_final"] = round(formed / possible, 3) if possible else None

    return summary


def bound_dimer_transforms(sim_dir: Path) -> Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]]:
    """Rigid transform of each bond type, taken from NERDSS's placement of isolated dimers.

    Only dimers are used: a two-molecule complex was placed by a diffusive association,
    which positions it exactly, never by a loop closure, which does not move anything.
    Returned as (rotation, translation) carrying the partner into the frame of the
    molecule holding the first-listed site.
    """
    model = Model.from_nerdss_dir(sim_dir)
    transforms: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]] = {}
    for _step, complexes in iter_snapshots(sim_dir):
        for complex_json in complexes:
            if len(complex_json["coords"]) != 2:
                continue
            bond = complex_json["bonds"][0]
            i, j, site_i, site_j = bond["molindex1"], bond["molindex2"], bond["site1"], bond["site2"]
            gap = model.site_gap(complex_json, bond)
            if gap is None or gap > STRETCHED_BOND_FACTOR:
                continue
            rotations = [quaternion_to_matrix(q) for q in complex_json["rotations"]]
            coords = np.asarray(complex_json["coords"], float)
            transforms.setdefault(
                (site_i, site_j),
                (rotations[i].T @ rotations[j], rotations[i].T @ (coords[j] - coords[i])),
            )
    return transforms


def closure_report(sim_dir: Path, system_json: Path) -> Optional[dict]:
    """How far the model is from being able to close its own rings.

    For a lattice design, a run of ``k`` short steps has to land where one long step of
    ``k`` times the size does. The report gives the mismatch between the two routes and
    whether the remaining site gap is inside NERDSS's loop-closure window, which is what
    decides whether the ring ever closes in a simulation.
    """
    lattice = designed_lattice(system_json)
    if lattice is None:
        return None
    model = Model.from_nerdss_dir(Path(sim_dir))
    transforms = bound_dimer_transforms(Path(sim_dir))
    sites = {site: step for site, step in lattice.step.items() if step > 0}
    if len(sites) < 2:
        return None

    short_site = min(sites, key=lambda site: sites[site])
    short_step = sites[short_site]
    report = {"short_site": short_site, "short_step": short_step}
    short = next((t for key, t in transforms.items() if key[0] == short_site), None)
    if short is None:
        return None

    for long_site, long_step in sorted(sites.items(), key=lambda item: item[1]):
        if long_site == short_site or long_step % short_step:
            continue
        long = next((t for key, t in transforms.items() if key[0] == long_site), None)
        if long is None:
            continue
        repeats = long_step // short_step

        rotation, translation = np.eye(3), np.zeros(3)
        for _ in range(repeats):
            rotation, translation = rotation @ short[0], translation + rotation @ short[1]

        partner_site = next(key[1] for key in transforms if key[0] == long_site)
        mol_name = next(iter(model.sites))
        gap = float(np.linalg.norm(
            model.sites[mol_name][long_site] - (translation + rotation @ model.sites[mol_name][partner_site])
        ))
        radius = model.radius(mol_name, long_site, mol_name, partner_site) or float("nan")
        cosine = (np.trace(rotation @ long[0].T) - 1.0) / 2.0
        report.update({
            "long_site": long_site,
            "long_step": long_step,
            "translation_error_nm": round(float(np.linalg.norm(translation - long[1])), 4),
            "rotation_error_deg": round(float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))), 3),
            "loop_gap_nm": round(gap, 4),
            "loop_closure_limit_nm": round(LOOP_CLOSURE_FACTOR * radius, 4),
            "closes": bool(gap <= LOOP_CLOSURE_FACTOR * radius),
        })
        break
    return report

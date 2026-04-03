"""
Utilities for "lego assembly" structure validation.

This validation mode reduces the designed assembly to one representative copy
per exported molecule type, turns binding effectively irreversible by forcing
all off-rates to zero, injects titration reactions so subunits can appear
gradually, and compares a final assembled structure against the designed
coarse-grained target with rigid alignment + RMSD.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, MutableMapping, Optional, Sequence, Tuple, Union
import json
import re
import shutil
import warnings
from collections import Counter, defaultdict

import numpy as np

from ionerdss.model.components.system import System
from ionerdss.model.components.instances import MoleculeInstance
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.utils.rigid_transform import apply_rigid_transform, rigid_transform_3d


CoordinateInput = Union[Mapping[str, Sequence[float]], np.ndarray, Sequence[Sequence[float]]]
DesignedStructureRecord = Dict[str, Any]


@dataclass(frozen=True)
class StructureValidationConfig:
    """Configuration for the coarse-grained structure validation protocol."""

    box_nm: Tuple[float, float, float] = (100.0, 100.0, 100.0)
    initial_molecule_count: int = 1
    titration_on_rate: Union[float, Dict[str, float]] = 1.0e-5
    target_filename: str = "structure_validation_target.json"
    titration_parms_filename: str = "parms_titrate.inp"


@dataclass(frozen=True)
class StructureValidationArtifacts:
    """Files and metadata generated for the validation simulation."""

    molecule_counts: Dict[str, int]
    target_counts: Dict[str, int]
    designed_coordinates: Dict[str, Tuple[float, float, float]]
    target_file: Path
    nerdss_files: Dict[str, Path]


@dataclass(frozen=True)
class StructureAlignmentResult:
    """Rigid-alignment result for a designed vs observed assembly."""

    labels: Tuple[str, ...]
    rmsd: float
    rotation_matrix: np.ndarray
    translation_vector: np.ndarray
    designed_coordinates: np.ndarray
    observed_coordinates: np.ndarray
    aligned_observed_coordinates: np.ndarray


@dataclass(frozen=True)
class StructureValidationSimulationResult:
    """Result of running a validation NERDSS simulation."""

    simulation_dir: Path
    histogram_file: Path
    final_coords_file: Path
    system_psf_file: Path
    restart_file: Path
    full_assembly_found: bool
    warning_message: Optional[str]
    first_full_assembly_time: Optional[float]
    observed_coordinates: Optional[Dict[str, Tuple[float, float, float]]]


def _as_xyz_array(coords: CoordinateInput, labels: Optional[Iterable[str]] = None) -> Tuple[Tuple[str, ...], np.ndarray]:
    """Normalize coordinate input into an ordered label tuple and an Nx3 array."""
    if isinstance(coords, Mapping):
        ordered_labels = tuple(sorted(coords))
        points = np.asarray([coords[label] for label in ordered_labels], dtype=float)
        return ordered_labels, points

    points = np.asarray(coords, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("Coordinate arrays must have shape (N, 3).")

    if labels is None:
        ordered_labels = tuple(str(i) for i in range(points.shape[0]))
    else:
        ordered_labels = tuple(labels)
        if len(ordered_labels) != points.shape[0]:
            raise ValueError("The number of labels must match the number of coordinates.")

    return ordered_labels, points


def get_representative_instances(system: System) -> Dict[str, MoleculeInstance]:
    """Choose one representative instance per molecule type.

    The representative is the instance with the largest number of interfaces.
    Ties are broken by molecule instance name to keep selection deterministic.
    """
    representatives: Dict[str, MoleculeInstance] = {}

    for molecule_instance in system.molecule_instances:
        molecule_type = molecule_instance.molecule_type
        if molecule_type is None:
            continue

        mol_name = molecule_type.name
        current = representatives.get(mol_name)
        if current is None:
            representatives[mol_name] = molecule_instance
            continue

        current_key = (len(current.interfaces_neighbors_map), current.name)
        candidate_key = (len(molecule_instance.interfaces_neighbors_map), molecule_instance.name)
        if candidate_key > current_key:
            representatives[mol_name] = molecule_instance

    return dict(sorted(representatives.items()))


def get_structure_validation_counts(system: System) -> Dict[str, int]:
    """Return the exact stoichiometry required for the designed assembly validation."""
    counts = {}
    for inst in system.molecule_instances:
        if inst.molecule_type:
            name = inst.molecule_type.name
            counts[name] = counts.get(name, 0) + 1
    return dict(sorted(counts.items()))


def build_validation_molecule_counts(system: System, initial_molecule_count: int = 1) -> Dict[str, int]:
    """Return validation counts with a configurable initial copy number per molecule type."""
    target_counts = get_structure_validation_counts(system)
    return {
        mol_name: int(initial_molecule_count) * count
        for mol_name, count in target_counts.items()
    }


def get_designed_structure(system: System) -> Dict[str, DesignedStructureRecord]:
    """Return the designed assembly with global COM and interface coordinates for each instance."""
    designed: Dict[str, DesignedStructureRecord] = {}

    for molecule_instance in system.molecule_instances:
        molecule_type = molecule_instance.molecule_type
        if molecule_type is None:
            continue

        interfaces = []
        for interface_instance, partner_instance in molecule_instance.interfaces_neighbors_map.items():
            if partner_instance is None:
                continue

            partner_type = getattr(partner_instance, "molecule_type", None)
            interface_coord = getattr(interface_instance, "absolute_coord", None)
            if interface_coord is None:
                continue

            interfaces.append(
                {
                    "interface_instance": interface_instance.get_name(),
                    "interface_type": interface_instance.interface_type.get_name()
                    if getattr(interface_instance, "interface_type", None) is not None
                    else None,
                    "binding_partner_instance": partner_instance.name,
                    "binding_partner_type": partner_type.name if partner_type is not None else None,
                    "global_interface_coord": tuple(float(value) for value in interface_coord),
                }
            )

        interfaces.sort(
            key=lambda item: (
                item["interface_instance"],
                item["interface_type"] or "",
                item["binding_partner_instance"],
                item["binding_partner_type"] or "",
                item["global_interface_coord"],
            )
        )

        designed[molecule_instance.name] = {
            "instance": molecule_instance.name,
            "type": molecule_type.name,
            "global_com_coord": tuple(float(value) for value in molecule_instance.com),
            "interfaces": interfaces,
        }

    return dict(sorted(designed.items()))


def _designed_structure_to_coordinate_map(
    designed_structure: Mapping[str, DesignedStructureRecord],
) -> Dict[str, Tuple[float, float, float]]:
    """Extract the global COM coordinate map used by validation target serialization."""
    return {
        instance_name: tuple(float(value) for value in record["global_com_coord"])
        for instance_name, record in designed_structure.items()
    }


def write_structure_validation_target(
    system: System,
    output_path: Union[str, Path],
    molecule_counts: Optional[Mapping[str, int]] = None,
    designed_coordinates: Optional[Mapping[str, Sequence[float]]] = None,
) -> Path:
    """Write the designed validation target to JSON."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "molecule_counts": dict(molecule_counts or get_structure_validation_counts(system)),
        "designed_coordinates": {
            key: tuple(float(value) for value in coords)
            for key, coords in (
                designed_coordinates
                or _designed_structure_to_coordinate_map(get_designed_structure(system))
            ).items()
        },
    }

    with open(output, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    return output


def prepare_structure_validation(
    system: System,
    workspace_manager=None,
    config: Optional[StructureValidationConfig] = None,
    parms_overrides: Optional[MutableMapping[str, object]] = None,
    designed_coordinates: Optional[Mapping[str, Sequence[float]]] = None,
) -> StructureValidationArtifacts:
    """Export the special NERDSS input deck for structure validation."""
    config = config or StructureValidationConfig()
    molecule_counts = build_validation_molecule_counts(
        system,
        initial_molecule_count=config.initial_molecule_count,
    )
    target_counts = get_structure_validation_counts(system)
    
    final_designed_coordinates = {
        key: tuple(float(value) for value in coords)
        for key, coords in (
            designed_coordinates
            or _designed_structure_to_coordinate_map(get_designed_structure(system))
        ).items()
    }

    exporter = NERDSSExporter(system, workspace_manager)
    export_overrides = dict(parms_overrides or {})
    export_overrides["force_off_ratekb"] = 0.0
    export_overrides["titration_on_rate_3d_ka"] = config.titration_on_rate

    nerdss_files = exporter.export_all(
        molecule_counts=molecule_counts,
        box_nm=config.box_nm,
        parms_overrides=export_overrides,
    )

    if workspace_manager is not None:
        target_file = workspace_manager.workspace_path / "nerdss_files" / config.target_filename
    else:
        target_file = Path("nerdss_files") / config.target_filename

    target_file = write_structure_validation_target(
        system=system,
        output_path=target_file,
        molecule_counts=target_counts,
        designed_coordinates=final_designed_coordinates,
    )

    parms_path = nerdss_files.get("parms")
    if parms_path is not None and config.titration_parms_filename:
        titration_parms_path = parms_path.with_name(config.titration_parms_filename)
        if titration_parms_path != parms_path:
            shutil.copyfile(parms_path, titration_parms_path)
            nerdss_files["parms_titrate"] = titration_parms_path

    return StructureValidationArtifacts(
        molecule_counts=molecule_counts,
        target_counts=target_counts,
        designed_coordinates=final_designed_coordinates,
        target_file=target_file,
        nerdss_files=nerdss_files,
    )


def align_structure_to_design(
    designed_coordinates: CoordinateInput,
    observed_coordinates: CoordinateInput,
    *,
    labels: Optional[Iterable[str]] = None,
    backend: str = "kabsch",
    plot: bool = False,
) -> StructureAlignmentResult:
    """Rigidly align an observed structure onto the designed target and compute RMSD."""
    designed_labels, designed_xyz = _as_xyz_array(designed_coordinates, labels=labels)
    observed_labels, observed_xyz = _as_xyz_array(observed_coordinates, labels=labels)

    if designed_labels != observed_labels:
        raise ValueError(
            "Designed and observed structures must have the same ordered labels. "
            "Pass dictionaries keyed by molecule type to match automatically."
        )

    if designed_xyz.shape != observed_xyz.shape:
        raise ValueError("Designed and observed structures must have the same shape.")

    if backend == "kabsch":
        rotation_matrix, translation_vector = rigid_transform_3d(observed_xyz, designed_xyz)
    elif backend == "biopython":
        try:
            from Bio.SVDSuperimposer import SVDSuperimposer
        except ImportError as exc:
            raise ImportError(
                "Biopython is required for backend='biopython'. Install biopython or use backend='kabsch'."
            ) from exc

        superimposer = SVDSuperimposer()
        superimposer.set(designed_xyz, observed_xyz)
        superimposer.run()
        rotation_matrix, translation_vector = superimposer.get_rotran()
    else:
        raise ValueError("backend must be 'kabsch' or 'biopython'.")
    aligned_observed = apply_rigid_transform(rotation_matrix, translation_vector, observed_xyz)
    deltas = aligned_observed - designed_xyz
    rmsd = float(np.sqrt(np.mean(np.sum(deltas * deltas, axis=1))))

    if plot:
        try:
            import matplotlib.pyplot as plt
            
            fig = plt.figure(figsize=(12, 10))
            
            # Plot 1: Unaligned Observed vs Designed
            ax1 = fig.add_subplot(221, projection='3d')
            ax1.scatter(designed_xyz[:, 0], designed_xyz[:, 1], designed_xyz[:, 2], 
                        c='blue', label='Designed', marker='o', alpha=0.7)
            ax1.scatter(observed_xyz[:, 0], observed_xyz[:, 1], observed_xyz[:, 2], 
                        c='red', label='Observed (Unaligned)', marker='^', alpha=0.7)
            ax1.set_title('Before Alignment')
            ax1.legend()
            
            # Plot 2: Aligned Observed vs Designed
            ax2 = fig.add_subplot(222, projection='3d')
            ax2.scatter(designed_xyz[:, 0], designed_xyz[:, 1], designed_xyz[:, 2], 
                        c='blue', label='Designed', marker='o', alpha=0.7)
            ax2.scatter(aligned_observed[:, 0], aligned_observed[:, 1], aligned_observed[:, 2], 
                        c='green', label='Observed (Aligned)', marker='^', alpha=0.7)
            ax2.set_title(f'After Alignment (RMSD: {rmsd:.4f} nm)')
            ax2.legend()
            
            # Plot 3: Designed Only
            ax3 = fig.add_subplot(223, projection='3d')
            ax3.scatter(designed_xyz[:, 0], designed_xyz[:, 1], designed_xyz[:, 2], 
                        c='blue', label='Designed', marker='o', alpha=0.7)
            ax3.set_title('After Alignment (Designed Only)')
            ax3.legend()

            # Plot 4: Observed Aligned Only
            ax4 = fig.add_subplot(224, projection='3d')
            ax4.scatter(aligned_observed[:, 0], aligned_observed[:, 1], aligned_observed[:, 2], 
                        c='green', label='Observed (Aligned)', marker='^', alpha=0.7)
            ax4.set_title('After Alignment (Observed Only)')
            ax4.legend()
            
            plt.tight_layout()
            plt.show()
        except ImportError:
            import logging
            logging.getLogger(__name__).warning("matplotlib is required for plotting structure alignment.")

    return StructureAlignmentResult(
        labels=designed_labels,
        rmsd=rmsd,
        rotation_matrix=rotation_matrix,
        translation_vector=translation_vector,
        designed_coordinates=designed_xyz,
        observed_coordinates=observed_xyz,
        aligned_observed_coordinates=aligned_observed,
    )


def _parse_psf_com_records(system_psf_file: Union[str, Path]) -> list[Tuple[int, int, str]]:
    """Return `(atom_index, mol_id, mol_name)` entries for COM atoms in PSF order."""
    records: list[Tuple[int, int, str]] = []
    natom = None
    with open(system_psf_file, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.endswith("!NATOM"):
                natom = int(stripped.split()[0])
                continue
            if natom is None or natom <= 0:
                continue

            parts = line.split()
            if len(parts) >= 5:
                atom_index = int(parts[0]) - 1
                mol_name = parts[1]
                mol_id = int(parts[2])
                atom_name = parts[3]
                if atom_name == "COM":
                    records.append((atom_index, mol_id, mol_name))
                natom -= 1
                if natom == 0:
                    break
    return records


def _parse_xyz_coordinates(xyz_file: Union[str, Path]) -> np.ndarray:
    """Parse a single XYZ frame into an Nx3 coordinate array."""
    with open(xyz_file, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    atom_count = int(lines[0].strip())
    coord_lines = lines[2:2 + atom_count]
    coords = []
    for line in coord_lines:
        parts = line.split()
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.asarray(coords, dtype=float)


def _parse_restart_snapshot(
    restart_file: Union[str, Path],
) -> Tuple[Dict[int, set[int]], Dict[int, Tuple[float, float, float]]]:
    """Parse molecule connectivity and COM coordinates from a NERDSS restart file."""
    lines = Path(restart_file).read_text(encoding="utf-8", errors="replace").splitlines()

    start_idx = None
    for idx, line in enumerate(lines):
        if line.strip() == "#All Molecules and coordinates":
            start_idx = idx
            break
    if start_idx is None:
        raise ValueError(f"Could not find molecule section in restart file: {restart_file}")

    header_parts = lines[start_idx + 1].split()
    if not header_parts:
        raise ValueError(f"Malformed molecule section header in restart file: {restart_file}")
    molecule_count = int(header_parts[0])

    adjacency: Dict[int, set[int]] = defaultdict(set)
    restart_coords: Dict[int, Tuple[float, float, float]] = {}
    idx = start_idx + 2
    for _ in range(molecule_count):
        header = lines[idx].split()
        if len(header) < 5:
            raise ValueError(f"Malformed molecule block header in restart file near line {idx + 1}")
        mol_id = int(header[0])
        idx += 1

        metadata_tokens = lines[idx].split()
        if len(metadata_tokens) < 6:
            raise ValueError(f"Malformed molecule metadata line in restart file near line {idx + 1}")
        idx += 1

        coord_tokens = lines[idx].split()
        if len(coord_tokens) < 3:
            raise ValueError(f"Malformed coordinate line in restart file near line {idx + 1}")
        restart_coords[mol_id] = (
            float(coord_tokens[0]),
            float(coord_tokens[1]),
            float(coord_tokens[2]),
        )
        idx += 1

        free_list_line = lines[idx].split()
        if not free_list_line:
            raise ValueError(f"Malformed free interface list in restart file near line {idx + 1}")
        free_list_size = int(free_list_line[0])
        idx += 1

        bound_list_line = lines[idx].split()
        if not bound_list_line:
            raise ValueError(f"Malformed bound interface list in restart file near line {idx + 1}")
        bound_list_size = int(bound_list_line[0])
        idx += 1

        partner_line = lines[idx].split()
        if not partner_line:
            raise ValueError(f"Malformed bound partner list in restart file near line {idx + 1}")
        partner_count = int(partner_line[0])
        partner_ids = [int(value) for value in partner_line[1:1 + partner_count]]
        idx += 1

        iface_count_line = lines[idx].split()
        if not iface_count_line:
            raise ValueError(f"Malformed interface count in restart file near line {idx + 1}")
        iface_count = int(iface_count_line[0])
        idx += 1

        if bound_list_size != partner_count:
            raise ValueError(
                "Mismatch between bound interface count and bound partner count in "
                f"{restart_file} near line {idx}."
            )

        for partner_id in partner_ids:
            if partner_id != mol_id:
                adjacency[mol_id].add(partner_id)
                adjacency[partner_id].add(mol_id)

        for _iface_idx in range(iface_count):
            iface_line = lines[idx].split()
            if len(iface_line) < 6:
                raise ValueError(f"Malformed interface header in restart file near line {idx + 1}")
            is_bound = int(iface_line[5])
            idx += 1

            iface_coord_line = lines[idx].split()
            if len(iface_coord_line) < 3:
                raise ValueError(f"Malformed interface coordinate line in restart file near line {idx + 1}")
            idx += 1

            if is_bound:
                bound_partner_line = lines[idx].split()
                if len(bound_partner_line) < 3:
                    raise ValueError(f"Malformed bound interface payload in restart file near line {idx + 1}")
                idx += 1

        for list_name in ("prevlist", "prevmyface", "prevpface", "prevnorm", "ps_prev", "prevsep"):
            list_line = lines[idx].split()
            if not list_line:
                raise ValueError(f"Malformed {list_name} list in restart file near line {idx + 1}")
            list_size = int(list_line[0])
            if len(list_line) < 1 + list_size:
                raise ValueError(f"Truncated {list_name} list in restart file near line {idx + 1}")
            idx += 1

    return adjacency, restart_coords


def _parse_restart_molecule_partners(restart_file: Union[str, Path]) -> Dict[int, set[int]]:
    """Parse final-frame molecule connectivity from a NERDSS restart file."""
    adjacency, _ = _parse_restart_snapshot(restart_file)
    return adjacency


def _restart_snapshot_sort_key(restart_path: Path) -> Tuple[Tuple[int, ...], str]:
    """Sort restart snapshot files from earliest to latest by embedded numeric suffixes."""
    numeric_parts = tuple(int(part) for part in re.findall(r"\d+", restart_path.stem))
    return numeric_parts, restart_path.name


def _iter_restart_snapshot_candidates(primary_restart_file: Union[str, Path]) -> list[Path]:
    """Return candidate restart snapshots, starting with DATA/restart.dat then RESTART newest-to-oldest."""
    primary = Path(primary_restart_file)
    candidates: list[Path] = []
    if primary.exists():
        candidates.append(primary)

    simulation_dir = primary.parent.parent
    restart_dir = simulation_dir / "RESTART"
    if restart_dir.is_dir():
        restart_files = [
            path for path in restart_dir.iterdir()
            if path.is_file() and path.suffix == ".dat"
        ]
        for path in sorted(restart_files, key=_restart_snapshot_sort_key, reverse=True):
            if path not in candidates:
                candidates.append(path)

    return candidates


def _extract_observed_com_coordinates(
    system_psf_file: Union[str, Path],
    final_coords_file: Optional[Union[str, Path]],
    restart_file: Union[str, Path],
    target_counts: Mapping[str, int],
) -> Dict[str, Tuple[float, float, float]]:
    """Extract COM coordinates from one connected final-frame assembly."""
    com_records = _parse_psf_com_records(system_psf_file)
    xyz_coords = None
    if final_coords_file is not None and Path(final_coords_file).exists():
        xyz_coords = _parse_xyz_coordinates(final_coords_file)
    adjacency, restart_coords = _parse_restart_snapshot(restart_file)

    mol_id_to_name: Dict[int, str] = {}
    mol_id_to_coord: Dict[int, Tuple[float, float, float]] = {}
    for atom_index, mol_id, mol_name in com_records:
        mol_id_to_name[mol_id] = mol_name
        if mol_id in restart_coords:
            mol_id_to_coord[mol_id] = restart_coords[mol_id]
        elif xyz_coords is not None:
            mol_id_to_coord[mol_id] = tuple(xyz_coords[atom_index].tolist())
        else:
            raise ValueError(f"Missing coordinates for molecule id {mol_id} in restart snapshot {restart_file}")
        adjacency.setdefault(mol_id, set())

    matching_components: list[list[int]] = []
    visited: set[int] = set()
    for mol_id in sorted(mol_id_to_name):
        if mol_id in visited:
            continue

        component: list[int] = []
        stack = [mol_id]
        visited.add(mol_id)
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in sorted(adjacency.get(current, ())):
                if neighbor in mol_id_to_name and neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)

        component_counts = Counter(mol_id_to_name[node_id] for node_id in component)
        if dict(component_counts) == dict(target_counts):
            matching_components.append(sorted(component))

    if not matching_components:
        raise ValueError(
            "No connected component in the final restart snapshot matches the designed "
            f"target composition {dict(target_counts)}."
        )

    selected_component = min(matching_components, key=lambda component: (component[0], component))
    if len(matching_components) > 1:
        warnings.warn(
            "Multiple full assemblies matching the validation target were present in the final "
            f"snapshot; selecting the component with the smallest molecule id set {selected_component}.",
            RuntimeWarning,
        )

    observed = {}
    type_counts: Dict[str, int] = {}
    for mol_id in selected_component:
        mol_name = mol_id_to_name[mol_id]
        copy_idx = type_counts.get(mol_name, 0)
        type_counts[mol_name] = copy_idx + 1
        
        if target_counts.get(mol_name, 1) > 1:
            key = f"{mol_name}_{copy_idx}"
        else:
            key = mol_name
            
        observed[key] = mol_id_to_coord[mol_id]

    missing = sorted(set(target_counts) - set(type_counts))
    if missing:
        raise ValueError(f"Missing COM coordinates for molecule types: {missing}")

    return observed


def _find_observed_com_coordinates_in_restart_snapshots(
    system_psf_file: Union[str, Path],
    final_coords_file: Optional[Union[str, Path]],
    restart_file: Union[str, Path],
    target_counts: Mapping[str, int],
) -> Tuple[Dict[str, Tuple[float, float, float]], Path]:
    """Search DATA and RESTART snapshots for a matching connected full assembly."""
    errors: list[str] = []
    for candidate_restart in _iter_restart_snapshot_candidates(restart_file):
        try:
            observed = _extract_observed_com_coordinates(
                system_psf_file=system_psf_file,
                final_coords_file=final_coords_file if candidate_restart == Path(restart_file) else None,
                restart_file=candidate_restart,
                target_counts=target_counts,
            )
            return observed, candidate_restart
        except Exception as exc:
            errors.append(f"{candidate_restart.name}: {exc}")

    raise ValueError(
        "No matching connected assembly was found in DATA/restart.dat or any RESTART snapshot. "
        + " | ".join(errors)
    )


def run_structure_validation_simulation(
    artifacts: StructureValidationArtifacts,
    nerdss_dir: Union[str, Path],
    *,
    sim_index: int = 1,
    sim_dir_name: str = "validation_output",
) -> StructureValidationSimulationResult:
    """Run a real NERDSS validation simulation and extract one full assembly if present."""
    from ionerdss.analysis.io.parser import parse_complex_histogram
    from ionerdss.nerdss_simulation import Simulation

    work_dir = Path(artifacts.nerdss_files["parms"]).parent
    simulation = Simulation(str(work_dir))
    simulation.parmfile = artifacts.nerdss_files.get("parms_titrate", artifacts.nerdss_files["parms"]).name
    simulation.run_new_simulations(
        sim_indices=[sim_index],
        sim_dir=str(work_dir / sim_dir_name),
        nerdss_dir=str(Path(nerdss_dir).expanduser()),
        parallel=False,
        progress=False,
        verbose=False,
    )

    simulation_dir = work_dir / sim_dir_name / str(sim_index)
    histogram_file = simulation_dir / "DATA" / "histogram_complexes_time.dat"
    final_coords_file = simulation_dir / "DATA" / "final_coords.xyz"
    system_psf_file = simulation_dir / "DATA" / "system.psf"
    restart_file = simulation_dir / "DATA" / "restart.dat"

    hist_times, hist_comps, hist_matrix = parse_complex_histogram(histogram_file)
    target_comp = dict(artifacts.target_counts)

    first_full_assembly_time = None
    full_assembly_found = False
    if len(hist_times) > 0:
        for col_idx, comp in enumerate(hist_comps):
            if comp == target_comp:
                counts = np.asarray(hist_matrix[:, [col_idx]].todense()).ravel()
                hits = np.where(counts > 0)[0]
                if len(hits) > 0:
                    full_assembly_found = True
                    first_full_assembly_time = float(hist_times[hits[0]])
                break

    warning_message = None
    observed_coordinates = None
    selected_restart_file = restart_file
    if full_assembly_found:
        try:
            observed_coordinates, selected_restart_file = _find_observed_com_coordinates_in_restart_snapshots(
                system_psf_file=system_psf_file,
                final_coords_file=final_coords_file,
                restart_file=restart_file,
                target_counts=artifacts.target_counts,
            )
        except Exception as exc:
            raise ValueError(
                "The target composition appeared in the histogram, but no matching connected "
                "assembly was found in DATA/restart.dat or any RESTART snapshot. "
                f"{exc}"
            ) from exc
    else:
        warning_message = (
            "Validation warning: no full assembly matching the designed one-copy target "
            f"{target_comp} was found in {histogram_file}."
        )
        warnings.warn(warning_message, RuntimeWarning)

    return StructureValidationSimulationResult(
        simulation_dir=simulation_dir,
        histogram_file=histogram_file,
        final_coords_file=final_coords_file,
        system_psf_file=system_psf_file,
        restart_file=selected_restart_file,
        full_assembly_found=full_assembly_found,
        warning_message=warning_message,
        first_full_assembly_time=first_full_assembly_time,
        observed_coordinates=observed_coordinates,
    )

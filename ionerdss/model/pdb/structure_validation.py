"""
Utilities for structure validation.

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
from itertools import permutations, product

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
    preflight_warning_message: Optional[str] = None


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
    largest_observed_assembly_size: int
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


def _strip_designed_label_to_type(label: str) -> str:
    """Reduce ionerdss-style labels like `chain_type` to `type`."""
    if "_" not in label:
        return label
    return label.split("_", 1)[1]


def _strip_observed_label_to_type(label: str) -> str:
    """Reduce NERDSS-style labels like `type_0` to `type`."""
    if "_" not in label:
        return label

    prefix, suffix = label.rsplit("_", 1)
    if suffix.isdigit():
        return prefix
    return label


def _compute_alignment(
    designed_xyz: np.ndarray,
    observed_xyz: np.ndarray,
    *,
    backend: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return rigid transform, aligned coordinates, and RMSD."""
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
    return rotation_matrix, translation_vector, aligned_observed, rmsd


def _match_coordinate_maps(
    designed_coordinates: CoordinateInput,
    observed_coordinates: CoordinateInput,
    *,
    labels: Optional[Iterable[str]] = None,
    backend: str,
) -> Tuple[Tuple[str, ...], np.ndarray, np.ndarray]:
    """Return labels and coordinate arrays with homomer-aware key matching for mappings."""
    designed_labels, designed_xyz = _as_xyz_array(designed_coordinates, labels=labels)
    observed_labels, observed_xyz = _as_xyz_array(observed_coordinates, labels=labels)

    if designed_xyz.shape != observed_xyz.shape:
        raise ValueError("Designed and observed structures must have the same shape.")

    if designed_labels == observed_labels:
        return designed_labels, designed_xyz, observed_xyz

    if not isinstance(designed_coordinates, Mapping) or not isinstance(observed_coordinates, Mapping):
        raise ValueError(
            "Designed and observed structures must have the same ordered labels. "
            "Pass dictionaries keyed by molecule labels to enable automatic matching."
        )

    designed_type_by_label = {
        label: _strip_designed_label_to_type(label)
        for label in designed_labels
    }
    observed_type_by_label = {
        label: _strip_observed_label_to_type(label)
        for label in observed_labels
    }

    designed_type_counts = Counter(designed_type_by_label.values())
    observed_type_counts = Counter(observed_type_by_label.values())
    if designed_type_counts != observed_type_counts:
        raise ValueError(
            "Designed and observed structures do not describe the same molecule-type composition after "
            "normalizing homomer labels."
        )

    designed_labels_by_type: Dict[str, list[str]] = defaultdict(list)
    observed_labels_by_type: Dict[str, list[str]] = defaultdict(list)
    for label in designed_labels:
        designed_labels_by_type[designed_type_by_label[label]].append(label)
    for label in observed_labels:
        observed_labels_by_type[observed_type_by_label[label]].append(label)

    for type_name in designed_type_counts:
        designed_labels_by_type[type_name].sort()
        observed_labels_by_type[type_name].sort()

    permutation_sets = [
        list(permutations(observed_labels_by_type[type_name]))
        for type_name, count in sorted(designed_type_counts.items())
        if count > 1
    ]
    repeated_types = [
        type_name
        for type_name, count in sorted(designed_type_counts.items())
        if count > 1
    ]

    best_labels: Optional[Tuple[str, ...]] = None
    best_observed_xyz: Optional[np.ndarray] = None
    best_rmsd: Optional[float] = None

    permutation_products = product(*permutation_sets) if permutation_sets else [()]
    for perm_choice in permutation_products:
        observed_order_by_type = {
            type_name: list(observed_labels_by_type[type_name])
            for type_name in designed_type_counts
        }
        for type_name, permuted_labels in zip(repeated_types, perm_choice):
            observed_order_by_type[type_name] = list(permuted_labels)

        matched_labels = tuple(designed_type_by_label[label] for label in designed_labels)
        matched_observed_labels = []
        type_offsets: Dict[str, int] = defaultdict(int)
        for designed_label in designed_labels:
            type_name = designed_type_by_label[designed_label]
            idx = type_offsets[type_name]
            matched_observed_labels.append(observed_order_by_type[type_name][idx])
            type_offsets[type_name] += 1

        candidate_observed_xyz = np.asarray(
            [observed_coordinates[label] for label in matched_observed_labels],
            dtype=float,
        )
        _, _, _, candidate_rmsd = _compute_alignment(
            designed_xyz,
            candidate_observed_xyz,
            backend=backend,
        )

        if best_rmsd is None or candidate_rmsd < best_rmsd:
            best_labels = matched_labels
            best_observed_xyz = candidate_observed_xyz
            best_rmsd = candidate_rmsd

    assert best_labels is not None
    assert best_observed_xyz is not None
    return best_labels, designed_xyz, best_observed_xyz


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


def _get_designed_connected_components(system: System) -> list[tuple[str, ...]]:
    """Return connected components of the designed assembly using molecule instance names."""
    adjacency: Dict[str, set[str]] = defaultdict(set)
    for molecule_instance in system.molecule_instances:
        adjacency.setdefault(molecule_instance.name, set())
        for _interface_instance, partner_instance in molecule_instance.interfaces_neighbors_map.items():
            if partner_instance is None:
                continue
            adjacency[molecule_instance.name].add(partner_instance.name)
            adjacency[partner_instance.name].add(molecule_instance.name)

    components: list[tuple[str, ...]] = []
    visited: set[str] = set()
    for instance_name in sorted(adjacency):
        if instance_name in visited:
            continue

        stack = [instance_name]
        visited.add(instance_name)
        component: list[str] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in sorted(adjacency[current]):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)

        components.append(tuple(sorted(component)))

    return components


def _format_disconnected_design_warning(components: Sequence[Sequence[str]]) -> Optional[str]:
    """Describe disconnected designed subunit groups for early validation feedback."""
    if len(components) <= 1:
        return None
    formatted_components = [
        ", ".join(component)
        for component in sorted((tuple(component) for component in components), key=lambda comp: comp)
    ]
    if len(formatted_components) == 2:
        return (
            "Validation preflight warning: the designed assembly graph is disconnected, so it cannot form a "
            f"single N-mer. Subunits {formatted_components[0]} are disconnected from subunits "
            f"{formatted_components[1]}."
        )

    return (
        "Validation preflight warning: the designed assembly graph is disconnected, so it cannot form a "
        f"single N-mer. Connected subunit groups: {'; '.join(formatted_components)}."
    )


def get_disconnected_design_message(system: System, *, prefix: str) -> Optional[str]:
    """Return a formatted disconnected-assembly message for the given system."""
    components = _get_designed_connected_components(system)
    if len(components) <= 1:
        return None

    formatted_components = [
        ", ".join(component)
        for component in sorted((tuple(component) for component in components), key=lambda comp: comp)
    ]
    if len(formatted_components) == 2:
        return (
            f"{prefix}: the designed assembly graph is disconnected, so it cannot form a single N-mer. "
            f"Subunits {formatted_components[0]} are disconnected from subunits {formatted_components[1]}."
        )

    return (
        f"{prefix}: the designed assembly graph is disconnected, so it cannot form a single N-mer. "
        f"Connected subunit groups: {'; '.join(formatted_components)}."
    )


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
    preflight_warning_message = _format_disconnected_design_warning(
        _get_designed_connected_components(system)
    )
    
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

    if preflight_warning_message:
        warnings.warn(preflight_warning_message, RuntimeWarning)

    return StructureValidationArtifacts(
        molecule_counts=molecule_counts,
        target_counts=target_counts,
        designed_coordinates=final_designed_coordinates,
        target_file=target_file,
        nerdss_files=nerdss_files,
        preflight_warning_message=preflight_warning_message,
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
    designed_labels, designed_xyz, observed_xyz = _match_coordinate_maps(
        designed_coordinates,
        observed_coordinates,
        labels=labels,
        backend=backend,
    )
    rotation_matrix, translation_vector, aligned_observed, rmsd = _compute_alignment(
        designed_xyz,
        observed_xyz,
        backend=backend,
    )

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


def _parse_restart_template_names_by_iface_count(lines: list[str]) -> Dict[int, str]:
    """Infer molecule names from the MolTemplates section keyed by unique interface counts."""
    molecule_section_idx = None
    for idx, line in enumerate(lines):
        stripped = line.strip()
        if stripped == "#All Molecules and coordinates":
            molecule_section_idx = idx
            break

    if molecule_section_idx is None:
        return {}

    names_by_iface_count: Dict[int, set[str]] = defaultdict(set)
    for line in lines[:molecule_section_idx]:
        parts = line.split()
        if (
            len(parts) != 3
            or not parts[0].isdigit()
            or not parts[2].isdigit()
            or not parts[1].isalpha()
        ):
            continue
        names_by_iface_count[int(parts[2])].add(parts[1])

    return {
        iface_count: next(iter(names))
        for iface_count, names in names_by_iface_count.items()
        if len(names) == 1
    }


def _parse_restart_snapshot(
    restart_file: Union[str, Path],
) -> Tuple[Dict[int, set[int]], Dict[int, Tuple[float, float, float]], Dict[int, str]]:
    """Parse molecule connectivity, COM coordinates, and inferred molecule names from a NERDSS restart file."""
    lines = Path(restart_file).read_text(encoding="utf-8", errors="replace").splitlines()
    template_names_by_iface_count = _parse_restart_template_names_by_iface_count(lines)

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
    restart_mol_names: Dict[int, str] = {}
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
        restart_mol_names[mol_id] = template_names_by_iface_count.get(iface_count, str(iface_count))
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

    return adjacency, restart_coords, restart_mol_names


def _parse_restart_molecule_partners(restart_file: Union[str, Path]) -> Dict[int, set[int]]:
    """Parse final-frame molecule connectivity from a NERDSS restart file."""
    adjacency, _, _ = _parse_restart_snapshot(restart_file)
    return adjacency


def _find_matching_components(
    adjacency: Mapping[int, set[int]],
    mol_id_to_name: Mapping[int, str],
    target_counts: Mapping[str, int],
) -> list[list[int]]:
    """Return connected components whose composition matches the target counts."""
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

    return matching_components


def _build_observed_coordinate_map(
    mol_names: Sequence[str],
    coords: Sequence[Sequence[float]],
    target_counts: Mapping[str, int],
) -> Dict[str, Tuple[float, float, float]]:
    """Build the observed coordinate map with stable per-copy labels for repeated molecule types."""
    observed: Dict[str, Tuple[float, float, float]] = {}
    type_counts: Dict[str, int] = {}
    for mol_name, coord in zip(mol_names, coords):
        copy_idx = type_counts.get(mol_name, 0)
        type_counts[mol_name] = copy_idx + 1
        key = f"{mol_name}_{copy_idx}" if target_counts.get(mol_name, 1) > 1 else mol_name
        observed[key] = tuple(float(value) for value in coord)

    missing = sorted(set(target_counts) - set(type_counts))
    if missing:
        raise ValueError(f"Missing COM coordinates for molecule types: {missing}")

    return observed


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


def _complex_snapshot_sort_key(snapshot_path: Path) -> Tuple[Tuple[int, ...], str]:
    """Sort COMPLEXES JSON snapshots by embedded numeric suffixes."""
    numeric_parts = tuple(int(part) for part in re.findall(r"\d+", snapshot_path.stem))
    return numeric_parts, snapshot_path.name


def _iter_complex_json_candidates(complexes_dir: Union[str, Path]) -> list[Path]:
    """Return COMPLEXES JSON snapshots newest-to-oldest."""
    directory = Path(complexes_dir)
    if not directory.is_dir():
        return []

    json_files = [path for path in directory.iterdir() if path.is_file() and path.suffix == ".json"]
    return sorted(json_files, key=_complex_snapshot_sort_key, reverse=True)


def _extract_observed_com_coordinates_from_complex_json(
    complex_json_file: Union[str, Path],
    target_counts: Mapping[str, int],
) -> Dict[str, Tuple[float, float, float]]:
    """Extract one matching complex directly from a COMPLEXES JSON snapshot."""
    payload = json.loads(Path(complex_json_file).read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise ValueError(f"Malformed COMPLEXES JSON snapshot: {complex_json_file}")

    matching_complexes: list[tuple[list[str], list[Sequence[float]]]] = []
    for complex_record in payload:
        if not isinstance(complex_record, dict):
            continue
        names = complex_record.get("names")
        coords = complex_record.get("coords")
        if not isinstance(names, list) or not isinstance(coords, list) or len(names) != len(coords):
            continue
        if dict(Counter(str(name) for name in names)) == dict(target_counts):
            matching_complexes.append(([str(name) for name in names], coords))

    if not matching_complexes:
        raise ValueError(
            "No complex in the COMPLEXES JSON snapshot matches the designed "
            f"target composition {dict(target_counts)}."
        )

    selected_names, selected_coords = matching_complexes[0]
    if len(matching_complexes) > 1:
        warnings.warn(
            "Multiple full assemblies matching the validation target were present in the COMPLEXES "
            f"snapshot {Path(complex_json_file).name}; selecting the first matching complex.",
            RuntimeWarning,
        )

    return _build_observed_coordinate_map(selected_names, selected_coords, target_counts)


def _find_observed_com_coordinates_in_complex_json_snapshots(
    complexes_dir: Union[str, Path],
    target_counts: Mapping[str, int],
) -> Tuple[Dict[str, Tuple[float, float, float]], Path]:
    """Search COMPLEXES JSON snapshots newest-to-oldest for a matching full assembly."""
    errors: list[str] = []
    for candidate_json in _iter_complex_json_candidates(complexes_dir):
        try:
            observed = _extract_observed_com_coordinates_from_complex_json(candidate_json, target_counts)
            return observed, candidate_json
        except Exception as exc:
            errors.append(f"{candidate_json.name}: {exc}")

    raise ValueError(
        "No matching assembly was found in DATA/COMPLEXES JSON snapshots. " + " | ".join(errors)
    )


def _extract_observed_com_coordinates(
    system_psf_file: Union[str, Path],
    final_coords_file: Optional[Union[str, Path]],
    restart_file: Union[str, Path],
    target_counts: Mapping[str, int],
) -> Dict[str, Tuple[float, float, float]]:
    """Extract COM coordinates from one connected final-frame assembly using restart-native typing only."""
    del system_psf_file, final_coords_file

    adjacency, restart_coords, restart_mol_names = _parse_restart_snapshot(restart_file)
    for mol_id in restart_coords:
        adjacency.setdefault(mol_id, set())

    matching_components = _find_matching_components(adjacency, restart_mol_names, target_counts)

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
    selected_names = [restart_mol_names[mol_id] for mol_id in selected_component]
    selected_coords = [restart_coords[mol_id] for mol_id in selected_component]
    observed = _build_observed_coordinate_map(selected_names, selected_coords, target_counts)
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


def _get_largest_complex_json_size(complex_json_file: Union[str, Path]) -> int:
    """Return the largest complex size present in a COMPLEXES JSON snapshot."""
    payload = json.loads(Path(complex_json_file).read_text(encoding="utf-8"))
    if not isinstance(payload, list):
        raise ValueError(f"Malformed COMPLEXES JSON snapshot: {complex_json_file}")

    largest_size = 0
    for complex_record in payload:
        if not isinstance(complex_record, dict):
            continue
        names = complex_record.get("names")
        coords = complex_record.get("coords")
        if not isinstance(names, list) or not isinstance(coords, list) or len(names) != len(coords):
            continue
        largest_size = max(largest_size, len(names))

    return largest_size


def _get_largest_complex_json_size_in_snapshots(complexes_dir: Union[str, Path]) -> int:
    """Return the largest complex size seen across COMPLEXES JSON snapshots."""
    largest_size = 0
    for candidate_json in _iter_complex_json_candidates(complexes_dir):
        try:
            largest_size = max(largest_size, _get_largest_complex_json_size(candidate_json))
        except Exception:
            continue
    return largest_size


def _get_largest_restart_component_size(restart_file: Union[str, Path]) -> int:
    """Return the largest connected component size in a restart snapshot."""
    adjacency, restart_coords, _restart_mol_names = _parse_restart_snapshot(restart_file)
    for mol_id in restart_coords:
        adjacency.setdefault(mol_id, set())

    largest_size = 0
    visited: set[int] = set()
    for mol_id in sorted(restart_coords):
        if mol_id in visited:
            continue

        stack = [mol_id]
        visited.add(mol_id)
        component_size = 0
        while stack:
            current = stack.pop()
            component_size += 1
            for neighbor in sorted(adjacency.get(current, ())):
                if neighbor in restart_coords and neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)

        largest_size = max(largest_size, component_size)

    return largest_size


def _get_largest_restart_component_size_in_snapshots(primary_restart_file: Union[str, Path]) -> int:
    """Return the largest connected component size seen across DATA and RESTART snapshots."""
    largest_size = 0
    for candidate_restart in _iter_restart_snapshot_candidates(primary_restart_file):
        try:
            largest_size = max(largest_size, _get_largest_restart_component_size(candidate_restart))
        except Exception:
            continue
    return largest_size


def run_structure_validation_simulation(
    artifacts: StructureValidationArtifacts,
    nerdss_dir: Union[str, Path],
    *,
    sim_index: int = 1,
    sim_dir_name: str = "validation_output",
    env: Optional[Mapping[str, str]] = None,
) -> StructureValidationSimulationResult:
    """Run a real NERDSS validation simulation and extract one full assembly if present.

    ``env`` holds environment variables the NERDSS executable needs, for example
    ``{"LD_LIBRARY_PATH": "/path/to/gsl/lib"}``. The entries are merged on top of the
    current ``os.environ``, so only the overrides need to be passed.
    """
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
        env=dict(env) if env is not None else None,
    )

    simulation_dir = work_dir / sim_dir_name / str(sim_index)
    histogram_file = simulation_dir / "DATA" / "histogram_complexes_time.dat"
    complexes_dir = simulation_dir / "DATA" / "COMPLEXES"
    final_coords_file = simulation_dir / "DATA" / "final_coords.xyz"
    system_psf_file = simulation_dir / "DATA" / "system.psf"
    restart_file = simulation_dir / "DATA" / "restart.dat"

    hist_times = np.asarray([])
    hist_comps = []
    hist_matrix = None
    if histogram_file.exists():
        hist_times, hist_comps, hist_matrix = parse_complex_histogram(histogram_file)
    target_comp = dict(artifacts.target_counts)

    first_full_assembly_time = None
    full_assembly_found = False
    if len(hist_times) > 0 and hist_matrix is not None:
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
    largest_observed_assembly_size = 0
    complex_json_candidates = _iter_complex_json_candidates(complexes_dir)
    if complex_json_candidates:
        largest_observed_assembly_size = _get_largest_complex_json_size_in_snapshots(complexes_dir)
        try:
            observed_coordinates, selected_restart_file = _find_observed_com_coordinates_in_complex_json_snapshots(
                complexes_dir=complexes_dir,
                target_counts=artifacts.target_counts,
            )
            full_assembly_found = True
        except Exception as exc:
            warning_message = (
                "Validation warning: no full assembly matching the designed one-copy target "
                f"{target_comp} was found in {complexes_dir}."
            )
            warnings.warn(warning_message, RuntimeWarning)
    else:
        fallback_warning = (
            "Validation warning: no COMPLEXES JSON snapshots were found in "
            f"{complexes_dir}; falling back to restart snapshots."
        )
        warnings.warn(fallback_warning, RuntimeWarning)
        warning_message = fallback_warning
        largest_observed_assembly_size = _get_largest_restart_component_size_in_snapshots(restart_file)
        try:
            observed_coordinates, selected_restart_file = _find_observed_com_coordinates_in_restart_snapshots(
                system_psf_file=system_psf_file,
                final_coords_file=final_coords_file,
                restart_file=restart_file,
                target_counts=artifacts.target_counts,
            )
            full_assembly_found = True
        except Exception:
            full_assembly_found = False

    if not full_assembly_found:
        no_assembly_warning = (
            "Validation warning: no full assembly matching the designed one-copy target "
            f"{target_comp} was found in {complexes_dir if complex_json_candidates else restart_file}."
        )
        warning_message = (
            no_assembly_warning
            if warning_message is None
            else f"{warning_message} {no_assembly_warning}"
        )
        if complex_json_candidates:
            warnings.warn(no_assembly_warning, RuntimeWarning)
    else:
        if warning_message is not None and not warning_message.endswith("snapshots."):
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
        largest_observed_assembly_size=largest_observed_assembly_size,
        observed_coordinates=observed_coordinates,
    )

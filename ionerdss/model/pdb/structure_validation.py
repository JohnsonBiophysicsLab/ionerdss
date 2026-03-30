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
from typing import Dict, Iterable, Mapping, MutableMapping, Optional, Sequence, Tuple, Union
import json
import shutil
import warnings

import numpy as np

from ionerdss.model.components.system import System
from ionerdss.model.components.instances import MoleculeInstance
from ionerdss.model.pdb.nerdss_exporter import NERDSSExporter
from ionerdss.utils.rigid_transform import apply_rigid_transform, rigid_transform_3d


CoordinateInput = Union[Mapping[str, Sequence[float]], np.ndarray, Sequence[Sequence[float]]]


@dataclass(frozen=True)
class StructureValidationConfig:
    """Configuration for the coarse-grained structure validation protocol."""

    box_nm: Tuple[float, float, float] = (100.0, 100.0, 100.0)
    initial_molecule_count: int = 1
    titration_on_rate: float = 1.0e-5
    target_filename: str = "structure_validation_target.json"
    titration_parms_filename: str = "parms_titrate.inp"


@dataclass(frozen=True)
class StructureValidationArtifacts:
    """Files and metadata generated for the validation simulation."""

    molecule_counts: Dict[str, int]
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
    """Return the validation stoichiometry: one copy per exported molecule type."""
    return {mol_name: 1 for mol_name in get_representative_instances(system)}


def build_validation_molecule_counts(system: System, initial_molecule_count: int = 1) -> Dict[str, int]:
    """Return validation counts with a configurable initial copy number per molecule type."""
    return {
        mol_name: int(initial_molecule_count)
        for mol_name in get_representative_instances(system)
    }


def get_designed_structure_coordinates(system: System) -> Dict[str, Tuple[float, float, float]]:
    """Return representative COM coordinates for the designed coarse-grained assembly."""
    representatives = get_representative_instances(system)
    return {
        mol_name: tuple(np.asarray(instance.com, dtype=float).tolist())
        for mol_name, instance in representatives.items()
    }


def write_structure_validation_target(
    system: System,
    output_path: Union[str, Path],
    molecule_counts: Optional[Mapping[str, int]] = None,
) -> Path:
    """Write the designed validation target to JSON."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "molecule_counts": dict(molecule_counts or get_structure_validation_counts(system)),
        "designed_coordinates": get_designed_structure_coordinates(system),
    }

    with open(output, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    return output


def prepare_structure_validation(
    system: System,
    workspace_manager=None,
    config: Optional[StructureValidationConfig] = None,
    parms_overrides: Optional[MutableMapping[str, object]] = None,
) -> StructureValidationArtifacts:
    """Export the special NERDSS input deck for structure validation."""
    config = config or StructureValidationConfig()
    molecule_counts = build_validation_molecule_counts(
        system,
        initial_molecule_count=config.initial_molecule_count,
    )

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
        molecule_counts=molecule_counts,
    )

    parms_path = nerdss_files.get("parms")
    if parms_path is not None and config.titration_parms_filename:
        titration_parms_path = parms_path.with_name(config.titration_parms_filename)
        if titration_parms_path != parms_path:
            shutil.copyfile(parms_path, titration_parms_path)
            nerdss_files["parms_titrate"] = titration_parms_path

    return StructureValidationArtifacts(
        molecule_counts=molecule_counts,
        designed_coordinates=get_designed_structure_coordinates(system),
        target_file=target_file,
        nerdss_files=nerdss_files,
    )


def align_structure_to_design(
    designed_coordinates: CoordinateInput,
    observed_coordinates: CoordinateInput,
    *,
    labels: Optional[Iterable[str]] = None,
    backend: str = "kabsch",
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

    return StructureAlignmentResult(
        labels=designed_labels,
        rmsd=rmsd,
        rotation_matrix=rotation_matrix,
        translation_vector=translation_vector,
        designed_coordinates=designed_xyz,
        observed_coordinates=observed_xyz,
        aligned_observed_coordinates=aligned_observed,
    )


def _parse_psf_com_records(system_psf_file: Union[str, Path]) -> list[Tuple[int, str]]:
    """Return `(atom_index, mol_name)` entries for COM atoms in PSF order."""
    records: list[Tuple[int, str]] = []
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
                atom_name = parts[3]
                if atom_name == "COM":
                    records.append((atom_index, mol_name))
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


def _extract_observed_com_coordinates(
    system_psf_file: Union[str, Path],
    final_coords_file: Union[str, Path],
    target_counts: Mapping[str, int],
) -> Dict[str, Tuple[float, float, float]]:
    """Extract COM coordinates keyed by molecule type from final NERDSS outputs."""
    if any(count != 1 for count in target_counts.values()):
        raise ValueError("Observed coordinate extraction currently requires one copy per molecule type.")

    com_records = _parse_psf_com_records(system_psf_file)
    xyz_coords = _parse_xyz_coordinates(final_coords_file)

    observed: Dict[str, Tuple[float, float, float]] = {}
    for atom_index, mol_name in com_records:
        if mol_name in target_counts and mol_name not in observed:
            observed[mol_name] = tuple(xyz_coords[atom_index].tolist())

    missing = sorted(set(target_counts) - set(observed))
    if missing:
        raise ValueError(f"Missing COM coordinates for molecule types: {missing}")

    return observed


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

    hist_times, hist_comps, hist_matrix = parse_complex_histogram(histogram_file)
    target_comp = dict(artifacts.molecule_counts)

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
    if full_assembly_found:
        observed_coordinates = _extract_observed_com_coordinates(
            system_psf_file=system_psf_file,
            final_coords_file=final_coords_file,
            target_counts=artifacts.molecule_counts,
        )
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
        full_assembly_found=full_assembly_found,
        warning_message=warning_message,
        first_full_assembly_time=first_full_assembly_time,
        observed_coordinates=observed_coordinates,
    )

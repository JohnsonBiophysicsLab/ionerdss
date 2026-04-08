"""
User-facing structure validation helpers for the PDB pipeline.

This module exposes the coarse-grained structure validation workflow as a
stable `ionerdss.model.pdb.validation` entry point.

Current logic:

1. Define the target composition from the designed one-copy validation system.
In structure_validation.py, get_structure_validation_counts() builds the expected full assembly as one representative copy of each molecule type, for example {"A": 1, "H": 1, "L": 1}.

2. Run the actual NERDSS validation simulation with that target in mind.
run_structure_validation_simulation(...) uses parms_titrate.inp, runs NERDSS, then looks for a matching full assembly in `DATA/COMPLEXES/*.json`. These JSON snapshots are the primary source for both existence checks and observed COM extraction.

3. If no COMPLEXES JSON snapshots exist, fall back to restart snapshots.
The code emits a warning and then scans `DATA/restart.dat` and any `RESTART/*.dat` snapshots for a connected component whose composition matches the target.

"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Sequence, Union

from ionerdss.model.components.system import System
from .structure_validation import (
    CoordinateInput,
    StructureAlignmentResult,
    StructureValidationArtifacts,
    StructureValidationConfig,
    StructureValidationSimulationResult,
    align_structure_to_design,
    build_validation_molecule_counts,
    get_designed_structure,
    get_representative_instances,
    get_structure_validation_counts,
    prepare_structure_validation,
    run_structure_validation_simulation,
    write_structure_validation_target,
)


def prepare(
    system: System,
    *,
    workspace_manager=None,
    box_nm: Sequence[float] = (100.0, 100.0, 100.0),
    titration_on_rate: float = 1.0e-5,
    target_filename: str = "structure_validation_target.json",
    parms_overrides: Optional[Dict[str, Any]] = None,
    designed_coordinates: Optional[Mapping[str, Sequence[float]]] = None,
) -> StructureValidationArtifacts:
    """Prepare the one-copy-per-type irreversible validation simulation."""
    config = StructureValidationConfig(
        box_nm=tuple(float(v) for v in box_nm),
        titration_on_rate=titration_on_rate,
        target_filename=target_filename,
    )
    return prepare_structure_validation(
        system=system,
        workspace_manager=workspace_manager,
        config=config,
        parms_overrides=parms_overrides,
        designed_coordinates=designed_coordinates,
    )


def setup_simulation(
    system: System,
    *,
    workspace_manager=None,
    box_nm: Sequence[float] = (100.0, 100.0, 100.0),
    initial_molecule_count: int = 1,
    titration_on_rate: Union[float, Dict[str, float]] = 1.0e-5,
    target_filename: str = "structure_validation_target.json",
    titration_parms_filename: str = "parms_titrate.inp",
    parms_overrides: Optional[Dict[str, Any]] = None,
    designed_coordinates: Optional[Mapping[str, Sequence[float]]] = None,
) -> StructureValidationArtifacts:
    """Set up the validation simulation with one of each, titration, and irreversible binding."""
    config = StructureValidationConfig(
        box_nm=tuple(float(v) for v in box_nm),
        initial_molecule_count=int(initial_molecule_count),
        titration_on_rate=titration_on_rate,
        target_filename=target_filename,
        titration_parms_filename=titration_parms_filename,
    )
    return prepare_structure_validation(
        system=system,
        workspace_manager=workspace_manager,
        config=config,
        parms_overrides=parms_overrides,
        designed_coordinates=designed_coordinates,
    )


def compare(
    designed_coordinates: Union[CoordinateInput, Mapping[str, Sequence[float]]],
    observed_coordinates: Union[CoordinateInput, Mapping[str, Sequence[float]]],
    *,
    backend: str = "kabsch",
    plot: bool = False,
) -> StructureAlignmentResult:
    """Align an observed structure onto the designed target and compute RMSD."""
    return align_structure_to_design(
        designed_coordinates=designed_coordinates,
        observed_coordinates=observed_coordinates,
        backend=backend,
        plot=plot,
    )


def align_structure(
    designed_coordinates: Union[CoordinateInput, Mapping[str, Sequence[float]]],
    observed_coordinates: Union[CoordinateInput, Mapping[str, Sequence[float]]],
    *,
    backend: str = "kabsch",
    plot: bool = False,
) -> StructureAlignmentResult:
    """Align an observed validation structure onto the designed target."""
    return align_structure_to_design(
        designed_coordinates=designed_coordinates,
        observed_coordinates=observed_coordinates,
        backend=backend,
        plot=plot,
    )


def run_simulation(
    artifacts: StructureValidationArtifacts,
    nerdss_dir: Union[str, Sequence[str]],
    *,
    sim_index: int = 1,
    sim_dir_name: str = "validation_output",
) -> StructureValidationSimulationResult:
    """Run the validation NERDSS job and extract one full assembly if it forms."""
    return run_structure_validation_simulation(
        artifacts=artifacts,
        nerdss_dir=nerdss_dir,
        sim_index=sim_index,
        sim_dir_name=sim_dir_name,
    )


__all__ = [
    "StructureAlignmentResult",
    "StructureValidationArtifacts",
    "StructureValidationConfig",
    "StructureValidationSimulationResult",
    "align_structure",
    "align_structure_to_design",
    "build_validation_molecule_counts",
    "compare",
    "get_designed_structure",
    "get_representative_instances",
    "get_structure_validation_counts",
    "prepare",
    "prepare_structure_validation",
    "run_simulation",
    "setup_simulation",
    "write_structure_validation_target",
]

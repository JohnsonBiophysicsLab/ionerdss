"""
User-facing structure validation helpers for the PDB pipeline.

This module exposes the coarse-grained structure validation workflow as a
stable `ionerdss.model.pdb.validation` entry point.

Current logic:

1. Define the target composition from the designed assembly.
In structure_validation.py, get_structure_validation_counts() counts every molecule instance of the designed system, so the expected full assembly carries the deposited stoichiometry, for example {"A": 1, "H": 1, "L": 1} for 8ERQ or {"A": 8} for 6BNO. The deck starts with these counts (times initial_molecule_count), titration adds further subunits, and the deck is written to its own directory, `structure_validation/` in the workspace by default.

2. Run the actual NERDSS validation simulation with that target in mind.
run_structure_validation_simulation(...) uses parms_titrate.inp, runs NERDSS, then looks for a matching full assembly in `DATA/COMPLEXES/*.json`. These JSON snapshots are the primary source for both existence checks and observed COM extraction; the deck sets `bondedComplexWrite` to nItr / 100 so that NERDSS writes them, unless `parms_overrides` sets it.

3. If no COMPLEXES JSON snapshots exist, fall back to restart snapshots.
The code emits a warning and then scans `DATA/restart.dat` and any `RESTART/*.dat` snapshots for a connected component whose composition matches the target.

"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Union

from ionerdss.model.components.system import System
from ionerdss.nerdss_simulation.executable import _merge_deprecated_nerdss_dir
from .structure_validation import (
    CoordinateInput,
    StructureAlignmentResult,
    StructureValidationArtifacts,
    StructureValidationConfig,
    StructureValidationSimulationResult,
    align_structure_to_design,
    build_validation_molecule_counts,
    classify_site_layout,
    collect_structure_validation_results,
    get_box_fit_message,
    get_degenerate_site_layouts,
    get_designed_assembly_extent,
    get_designed_structure,
    get_interface_com_proximity_message,
    get_near_com_interface_sites,
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
    interface_com_proximity_threshold_nm: Optional[float] = None,
    deck_dir: Union[str, Path] = "structure_validation",
) -> StructureValidationArtifacts:
    """Prepare the irreversible, titrated validation simulation for one copy of the designed assembly.

    ``interface_com_proximity_threshold_nm`` sets how close to its molecule's centre
    of mass a reacting interface site may sit before the preflight check reports it;
    ``None`` uses the hyperparameters in ``parms_overrides['hyperparams']`` or the
    module default. ``deck_dir`` is the directory the deck is written to, relative to
    the workspace (or the current directory without a ``workspace_manager``).
    """
    config = StructureValidationConfig(
        box_nm=tuple(float(v) for v in box_nm),
        titration_on_rate=titration_on_rate,
        target_filename=target_filename,
        interface_com_proximity_threshold_nm=interface_com_proximity_threshold_nm,
        deck_dir=deck_dir,
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
    interface_com_proximity_threshold_nm: Optional[float] = None,
    deck_dir: Union[str, Path] = "structure_validation",
) -> StructureValidationArtifacts:
    """Set up the titrated, irreversible validation simulation of the designed assembly.

    The deck starts with every molecule type as often as it occurs in the designed
    assembly, multiplied by ``initial_molecule_count``.

    ``interface_com_proximity_threshold_nm`` sets how close to its molecule's centre
    of mass a reacting interface site may sit before the preflight check reports it;
    ``None`` uses the hyperparameters in ``parms_overrides['hyperparams']`` or the
    module default. ``deck_dir`` is the directory the deck is written to, relative to
    the workspace (or the current directory without a ``workspace_manager``).
    """
    config = StructureValidationConfig(
        box_nm=tuple(float(v) for v in box_nm),
        initial_molecule_count=int(initial_molecule_count),
        titration_on_rate=titration_on_rate,
        target_filename=target_filename,
        titration_parms_filename=titration_parms_filename,
        interface_com_proximity_threshold_nm=interface_com_proximity_threshold_nm,
        deck_dir=deck_dir,
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
    nerdss_path: Optional[Union[str, Path]] = None,
    *,
    sim_index: int = 1,
    sim_dir_name: str = "validation_output",
    env: Optional[Mapping[str, str]] = None,
    nerdss_dir: Optional[Union[str, Path]] = None,
) -> StructureValidationSimulationResult:
    """Run the validation NERDSS job and extract one full assembly if it forms.

    ``nerdss_path`` is the NERDSS executable, or a directory containing ``nerdss`` or
    ``nerdss_mpi`` directly or in its ``bin/`` (for example a NERDSS checkout); ``None``
    looks the executable up on ``PATH``. The executable is run in place. ``nerdss_dir``
    is a deprecated alias for it.

    ``env`` holds environment variables the NERDSS executable needs, for example
    ``{"LD_LIBRARY_PATH": "/path/to/gsl/lib"}``. The entries are merged on top of the
    current ``os.environ``, so only the overrides need to be passed.
    """
    nerdss_path = _merge_deprecated_nerdss_dir(nerdss_path, nerdss_dir, "run_simulation")
    return run_structure_validation_simulation(
        artifacts=artifacts,
        nerdss_path=nerdss_path,
        sim_index=sim_index,
        sim_dir_name=sim_dir_name,
        env=env,
    )


def collect_results(
    artifacts: StructureValidationArtifacts,
    simulation_dir: Union[str, Path],
) -> StructureValidationSimulationResult:
    """Read the results of a validation NERDSS run that was launched outside ioNERDSS.

    Use this when the NERDSS executable is driven directly -- a manual
    ``subprocess.run``, a cluster submission, or a run from an earlier session --
    instead of through :func:`run_simulation`. It returns the same
    :class:`StructureValidationSimulationResult` that :func:`run_simulation` returns, so
    the downstream reporting and :func:`align_structure` steps are unchanged.

    ``simulation_dir`` is the directory the NERDSS run wrote its ``DATA`` directory
    into (the ``cwd`` of the external run); passing the ``DATA`` directory itself also
    works.

    Example:
        subprocess.run(f"{nerdss_cmd} -f parms_titrate.inp", shell=True, cwd=run_dir, env=env)
        sim_result = pdb.validation.collect_results(artifacts, simulation_dir=run_dir)
    """
    return collect_structure_validation_results(
        artifacts=artifacts,
        simulation_dir=simulation_dir,
    )


__all__ = [
    "StructureAlignmentResult",
    "StructureValidationArtifacts",
    "StructureValidationConfig",
    "StructureValidationSimulationResult",
    "align_structure",
    "align_structure_to_design",
    "build_validation_molecule_counts",
    "classify_site_layout",
    "collect_results",
    "collect_structure_validation_results",
    "compare",
    "get_box_fit_message",
    "get_degenerate_site_layouts",
    "get_designed_assembly_extent",
    "get_designed_structure",
    "get_interface_com_proximity_message",
    "get_near_com_interface_sites",
    "get_representative_instances",
    "get_structure_validation_counts",
    "prepare",
    "prepare_structure_validation",
    "run_simulation",
    "setup_simulation",
    "write_structure_validation_target",
]

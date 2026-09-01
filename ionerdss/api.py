"""
ionerdss.api - Simplified Public API

Provides convenience wrappers for common workflows.
"""

from typing import Optional, Dict, Any, Mapping, Sequence, Union
from pathlib import Path

from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.components.system import System
from ionerdss.model.pdb.structure_validation import (
    StructureAlignmentResult,
    StructureValidationArtifacts,
    StructureValidationConfig,
    align_structure_to_design,
    prepare_structure_validation,
)


def build_system_from_pdb(
    source: str,
    workspace_path: Optional[str] = None,
    fetch_format: Optional[str] = None,
    molecule_counts: Optional[Dict[str, int]] = None,
    structure_validation: bool = False,
    structure_validation_options: Optional[Dict[str, Any]] = None,
    **hyperparams_kwargs
) -> System:
    """Build ionerdss System from PDB structure (simplified API).
    
    This is a convenience function that combines PDBModelBuilder initialization
    and system building into a single call. All hyperparameter options can be
    passed as keyword arguments.
    
    Args:
        source: PDB ID (e.g., "4v6x") or path to PDB/mmCIF file.
        workspace_path: Workspace directory path. Defaults to "{source}_dir".
        fetch_format: Format for downloading structures ('pdb' or 'mmcif'). If None, uses hyperparameter default (usually 'bioassembly1').
        molecule_counts: Molecule counts for NERDSS export. Default 10 per type.
        structure_validation: Export the one-copy validation setup during build.
        structure_validation_options: Options passed to the validation export.
        **hyperparams_kwargs: Any PDBModelHyperparameters field as keyword arguments.
            Common options:
                - interface_detect_distance_cutoff: float (default 0.9)
                - generate_nerdss_files: bool (default True)
                - nerdss_water_box: list[float] (default [100, 100, 100])
                - ode_enabled: bool (default False)
                - ode_time_span: tuple[float, float]
                - ode_solver_method: str
                - ode_plot: bool
                - ode_save_csv: bool
    
    Returns:
        Complete System object ready for simulation.
    
    Examples:
        >>> # Simple usage with PDB ID
        >>> from ionerdss import build_system_from_pdb
        >>> system = build_system_from_pdb("4v6x")
        
        >>> # With custom parameters
        >>> system = build_system_from_pdb(
        ...     source="4v6x",
        ...     workspace_path="my_workspace",
        ...     interface_detect_distance_cutoff=1.0,
        ...     nerdss_water_box=[500, 500, 500],
        ...     ode_enabled=True,
        ...     ode_time_span=(0.0, 10.0)
        ... )
        
        >>> # From local file
        >>> system = build_system_from_pdb(
        ...     source="/path/to/structure.cif",
        ...     ode_enabled=True
        ... )
    """
    # Default workspace path
    if workspace_path is None:
        source_name = Path(source).stem if '/' in str(source) or '\\' in str(source) else source
        workspace_path = f"{source_name}_dir"
    
    # Create hyperparameters from kwargs
    hyperparams = PDBModelHyperparameters(**hyperparams_kwargs) if hyperparams_kwargs else None
    
    # Create builder
    builder = PDBModelBuilder(source=source, fetch_format=fetch_format, hyperparams=hyperparams)
    
    # Build system
    system = builder.build_system(
        workspace_path=workspace_path,
        hyperparams=hyperparams,
        molecule_counts=molecule_counts,
        structure_validation=structure_validation,
        structure_validation_options=structure_validation_options,
    )
    
    return system


def prepare_structure_validation_for_system(
    system: System,
    *,
    box_nm: Sequence[float] = (100.0, 100.0, 100.0),
    titration_on_rate: float = 1.0e-5,
    target_filename: str = "structure_validation_target.json",
    parms_overrides: Optional[Dict[str, Any]] = None,
) -> StructureValidationArtifacts:
    """Prepare the irreversible, titrated one-copy-per-type validation setup."""
    config = StructureValidationConfig(
        box_nm=tuple(float(v) for v in box_nm),
        titration_on_rate=titration_on_rate,
        target_filename=target_filename,
    )
    return prepare_structure_validation(
        system=system,
        workspace_manager=None,
        config=config,
        parms_overrides=parms_overrides,
    )


def compare_structure_to_design(
    designed_coordinates: Union[Mapping[str, Sequence[float]], Sequence[Sequence[float]]],
    observed_coordinates: Union[Mapping[str, Sequence[float]], Sequence[Sequence[float]]],
) -> StructureAlignmentResult:
    """Align an observed coarse-grained assembly onto the designed target and report RMSD."""
    return align_structure_to_design(
        designed_coordinates=designed_coordinates,
        observed_coordinates=observed_coordinates,
    )

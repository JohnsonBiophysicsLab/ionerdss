"""
ionerdss.api - Simplified Public API

Provides convenience wrappers for common workflows.
"""

from typing import Optional, Dict, Any
from pathlib import Path

from ionerdss.model.pdb.main import PDBModelBuilder
from ionerdss.model.pdb.hyperparameters import PDBModelHyperparameters
from ionerdss.model.components.system import System


def build_system_from_pdb(
    source: str,
    workspace_path: Optional[str] = None,
    fetch_format: Optional[str] = None,
    molecule_counts: Optional[Dict[str, int]] = None,
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
        **hyperparams_kwargs: Any PDBModelHyperparameters field as keyword arguments.
            Common options:
                - interface_detect_distance_cutoff: float (default 0.6)
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
        molecule_counts=molecule_counts
    )
    
    return system

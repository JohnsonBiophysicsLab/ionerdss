"""
ionerdss.model.pdb.api

High-level API for PDB model configuration and hyperparameter management.

This module provides convenient functions for setting up and managing
hyperparameters without requiring users to directly import or instantiate
the PDBModelHyperparameters class. Each function takes the PDBModelBuilder
whose hyperparameters it manages as its first argument; PDBModelBuilder has
methods of the same names that pass the builder itself.

## Quick Start

```python
from ionerdss.model import pdb

builder = pdb.PDBModelBuilder("1ABC")

# Set hyperparameters on the builder (creates or updates)
builder.set_hyperparameters(
    interface_detect_distance_cutoff=0.8,
    interface_detect_n_residue_cutoff=5,
    chain_grouping_matching_mode="sequence"
)
# equivalent: pdb.set_hyperparameters(builder, interface_detect_distance_cutoff=0.8, ...)

# Build model (the builder's hyperparameters are used)
system = builder.build_system(workspace_path="./workspace")
```

## Configuration Management

```python
# Save configuration
pdb.export_hyperparameters(builder, "config.json")

# Load configuration
pdb.import_hyperparameters(builder, "config.json")

# View current configuration
pdb.print_hyperparameters(builder)
```
"""

from typing import Optional, Dict, Any, TYPE_CHECKING
from pathlib import Path
from dataclasses import fields
import json

from .hyperparameters import PDBModelHyperparameters

if TYPE_CHECKING:
    from .main import PDBModelBuilder


# Sections of the hyperparameter reference in set_hyperparameters.__doc__ and of
# print_hyperparameters, following the field order of PDBModelHyperparameters.
# Both list any field missing from every section under "Other", so a new field
# cannot drop out of either.
_HYPERPARAMETER_SECTIONS = {
    "Core Detection": [
        "interface_detect_distance_cutoff",
        "interface_detect_n_residue_cutoff",
        "min_chain_length",
        "include_modified_residues",
    ],
    "Interface Type Assignment": [
        "interface_type_assignment_distance_threshold",
        "interface_type_assignment_angle_threshold",
        "split_spatial_interface_patches",
    ],
    "Interface Site Placement and Preflight": [
        "interface_com_proximity_threshold",
        "interface_site_placement",
    ],
    "Chain Grouping": [
        "chain_grouping_rmsd_threshold",
        "chain_grouping_seq_threshold",
        "chain_grouping_custom_aligner",
        "chain_grouping_matching_mode",
    ],
    "Steric Clash Detection": [
        "steric_clash_mode",
    ],
    "Template Building": [
        "signature_precision",
        "homodimer_distance_threshold",
        "homodimer_angle_threshold",
    ],
    "Homotypic Detection": [
        "homotypic_detection",
        "homotypic_detection_residue_similarity_threshold",
        "homotypic_detection_interface_radius",
    ],
    "Geometric Regularization": [
        "geometric_regularization",
        "symmetry_fold_tolerance",
        "com_shift_cap_ang",
        "is_on_sphere",
    ],
    "Template Regularization": [
        "template_regularization_strength",
    ],
    "Output Options": [
        "generate_visualizations",
        "generate_nerdss_files",
    ],
    "NERDSS Simulation": [
        "nerdss_water_box",
        "nerdss_total_molecule_count",
        "nerdss_time_step",
        "nerdss_n_itr",
        "nerdss_overlap_sep_limit",
        "disable_overlap_sep_limit_check",
    ],
    "ProAffinity Binding Energy Prediction": [
        "predict_affinity",
        "adfr_path",
        "proaffinity_backend",
        "proaffinity_python",
    ],
    "Structure Input": [
        "pdb_file_format",
    ],
    "ODE Pipeline": [
        "ode_enabled",
        "max_complex_size_ode",
        "ode_time_span",
        "ode_solver_method",
        "ode_atol",
        "ode_plot",
        "ode_save_csv",
        "ode_initial_concentrations",
        "ode_plot_species_indices",
        "ode_plot_sample_points",
        "ode_species_labels",
    ],
    "Kinetic Parameters": [
        "default_on_rate_3d_ka",
    ],
    "Transition Matrix": [
        "count_transition",
        "transition_matrix_size",
        "transition_write",
    ],
}


def _hyperparameter_sections() -> Dict[str, list]:
    """Return _HYPERPARAMETER_SECTIONS plus an "Other" section for unlisted fields.

    The internal ``units`` field is never listed.
    """
    listed = {name for names in _HYPERPARAMETER_SECTIONS.values() for name in names}
    unlisted = [
        field_info.name for field_info in fields(PDBModelHyperparameters)
        if field_info.name != 'units' and field_info.name not in listed
    ]
    sections = dict(_HYPERPARAMETER_SECTIONS)
    if unlisted:
        sections["Other"] = unlisted
    return sections


def _generate_hyperparameters_docstring() -> str:
    """Generate docstring for set_hyperparameters from PDBModelHyperparameters metadata.
    
    Returns:
        Complete docstring with parameter descriptions extracted from field metadata.
    """
    categories = _hyperparameter_sections()
    
    # Build field metadata dictionary
    field_metadata = {}
    for field_info in fields(PDBModelHyperparameters):
        if field_info.name == 'units':
            continue  # Skip units field
        field_metadata[field_info.name] = {
            'type': field_info.type,
            'default': field_info.default if field_info.default is not field_info.default_factory else field_info.default_factory(),
            'metadata': field_info.metadata
        }
    
    # Build docstring
    lines = [
        "Set or update hyperparameters for a PDBModelBuilder instance.",
        "",
        "Creates new hyperparameters if none exist on the builder, or updates existing ones.",
        "These hyperparameters are automatically used when calling builder.build_system().",
        "Keyword names must match PDBModelHyperparameters fields exactly.",
        "",
        "Hyperparameters Reference",
        "-------------------------",
        "",
    ]
    
    # Add parameter documentation by category
    for category, field_names in categories.items():
        lines.append(f"**{category}:**")
        for field_name in field_names:
            if field_name not in field_metadata:
                continue
            
            fmeta = field_metadata[field_name]
            type_str = str(fmeta['type']).replace('typing.', '').replace('<class \'', '').replace('\'>', '')
            default_val = fmeta['default']
            
            # Format default value
            if default_val is None:
                default_str = "None"
            elif isinstance(default_val, str):
                default_str = f'"{default_val}"'
            elif isinstance(default_val, tuple):
                default_str = str(default_val)
            else:
                default_str = str(default_val)
            
            # Get description and unit from metadata
            description = fmeta['metadata'].get('description', 'No description')
            unit = fmeta['metadata'].get('unit', '')
            
            if unit:
                param_line = f"- {field_name} ({type_str}, default={default_str}): {description} [{unit}]"
            else:
                param_line = f"- {field_name} ({type_str}, default={default_str}): {description}"
            
            lines.append(param_line)
        lines.append("")
    
    # Add Args, Returns, Examples sections
    lines.extend([
        "Args:",
        "    builder: PDBModelBuilder instance to configure (omitted when called as",
        "        builder.set_hyperparameters(...)).",
        "    **kwargs: Hyperparameter field names and values to set or update.",
        "",
        "Returns:",
        "    The updated PDBModelHyperparameters instance.",
        "",
        "Examples:",
        "    >>> from ionerdss.model import pdb",
        "    >>> ",
        "    >>> # Create builder",
        "    >>> builder = pdb.PDBModelBuilder('1ABC')",
        "    >>> ",
        "    >>> # Set hyperparameters with defaults",
        "    >>> builder.set_hyperparameters()",
        "    >>> ",
        "    >>> # Customize specific parameters",
        "    >>> builder.set_hyperparameters(",
        "    ...     interface_detect_distance_cutoff=0.8,",
        "    ...     interface_detect_n_residue_cutoff=5,",
        "    ...     chain_grouping_matching_mode='sequence'",
        "    ... )",
        "    >>> ",
        "    >>> # The module-level function takes the builder first",
        "    >>> pdb.set_hyperparameters(builder, nerdss_total_molecule_count=150)",
        "    >>> ",
        "    >>> # Enable advanced features",
        "    >>> builder.set_hyperparameters(",
        "    ...     steric_clash_mode='auto',",
        "    ...     geometric_regularization='auto',",
        "    ...     homotypic_detection='signature',",
        "    ...     ode_enabled=True,",
        "    ...     predict_affinity=True",
        "    ... )",
        "    >>> ",
        "    >>> # Build model (hyperparameters automatically used)",
        "    >>> system = builder.build_system(workspace_path='./workspace')",
        "",
        "Note:",
        "    These hyperparameters are automatically used by builder.build_system()",
        "    so you don't need to explicitly provide them.",
    ])
    
    return "\n".join(lines)


def set_hyperparameters(builder: 'PDBModelBuilder', **kwargs) -> PDBModelHyperparameters:
    if builder.hyperparams is None:
        # Create new hyperparameters
        builder.hyperparams = PDBModelHyperparameters(**kwargs)
    else:
        # Update existing hyperparameters
        current_config = builder.hyperparams.to_dict()
        current_config.update(kwargs)
        builder.hyperparams = PDBModelHyperparameters.from_dict(current_config)
    
    return builder.hyperparams


# Dynamically set docstring from field metadata
set_hyperparameters.__doc__ = _generate_hyperparameters_docstring()


def export_hyperparameters(builder: 'PDBModelBuilder', filepath: str) -> Dict[str, Any]:
    """Export builder's hyperparameters to JSON file.
    
    Args:
        builder: PDBModelBuilder instance.
        filepath: Path to save JSON file.
    
    Returns:
        Dictionary representation of hyperparameters.
    
    Raises:
        ValueError: If no hyperparameters have been set.
    
    Examples:
        >>> from ionerdss.model import pdb
        >>> 
        >>> # Create builder and set hyperparameters
        >>> builder = pdb.PDBModelBuilder("1ABC")
        >>> builder.set_hyperparameters(interface_detect_distance_cutoff=0.8)
        >>> 
        >>> # Export to file
        >>> builder.export_hyperparameters("config.json")
    """
    if builder.hyperparams is None:
        raise ValueError("No hyperparameters have been set. Call set_hyperparameters() first.")
    
    config = builder.hyperparams.to_dict()
    
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        json.dump(config, f, indent=2)
    
    return config


def import_hyperparameters(builder: 'PDBModelBuilder', filepath: str) -> PDBModelHyperparameters:
    """Import hyperparameters from JSON file and set on builder.
    
    Args:
        builder: PDBModelBuilder instance.
        filepath: Path to JSON file containing hyperparameters.
    
    Returns:
        The loaded PDBModelHyperparameters instance.
    
    Examples:
        >>> from ionerdss.model import pdb
        >>> 
        >>> # Create builder and load configuration from file
        >>> builder = pdb.PDBModelBuilder("1ABC")
        >>> builder.import_hyperparameters("config.json")
        >>> 
        >>> # Build model (loaded hyperparameters automatically used)
        >>> system = builder.build_system(workspace_path="./workspace")
    """
    path = Path(filepath)
    with open(path, 'r') as f:
        config = json.load(f)
    
    builder.hyperparams = PDBModelHyperparameters.from_dict(config)
    return builder.hyperparams


def print_hyperparameters(builder: 'PDBModelBuilder') -> str:
    """Print builder's hyperparameters in a human-readable format.
    
    Args:
        builder: PDBModelBuilder instance.
    
    Returns:
        String representation of hyperparameters.
    
    Raises:
        ValueError: If no hyperparameters have been set.
    
    Examples:
        >>> from ionerdss.model import pdb
        >>> 
        >>> # Create builder and set hyperparameters
        >>> builder = pdb.PDBModelBuilder("1ABC")
        >>> builder.set_hyperparameters()
        >>> 
        >>> # View configuration
        >>> print(builder.print_hyperparameters())
    """
    if builder.hyperparams is None:
        raise ValueError("No hyperparameters have been set. Call set_hyperparameters() first.")
    
    config = builder.hyperparams.to_dict()
    lines = []
    lines.append("PDB Model Hyperparameters")
    lines.append("=" * 50)
    lines.append("")
    
    # Group parameters by category
    categories = _hyperparameter_sections()
    
    for category, params in categories.items():
        lines.append(f"{category}:")
        for param in params:
            if param in config and param != 'chain_grouping_custom_aligner' and param != 'units':
                value = config[param]
                lines.append(f"  {param}: {value}")
        lines.append("")
    
    result = "\n".join(lines)
    print(result)
    return result

__all__ = [
    'set_hyperparameters',
    'export_hyperparameters',
    'import_hyperparameters',
    'print_hyperparameters',
]
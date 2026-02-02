"""
ODE Pipeline for ionerdss

This module provides functionality to calculate ODE solutions for molecular assembly
reactions before running NERDSS simulations. It integrates the graph-based reaction
network generator with the ODE solver.

Author: ionerdss team
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Union
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import csv

from typing import Any
from ionerdss.ode.reaction_string_parser import ReactionStringParser
from ionerdss.ode.reaction_ode_solver import solve_reaction_ode, dydt


@dataclass
class ODEPipelineConfig:
    """
    Configuration for ODE pipeline calculations.
    
    Attributes:
        t_span: Time span for integration [start, end] (default: [0.0, 10.0])
        initial_concentrations: Initial concentrations for species as dict {species_name: concentration}
                               If None, assumes first complex (monomer) at 1.0, others at 0.0
        solver_method: ODE solver method (default: "BDF" for stiff systems)
        atol: Absolute tolerance for solver (default: 1e-4)
        plot: Whether to generate plots (default: True)
        plot_species_indices: Indices of species to plot. If None, plots all (default: None)
        plot_sample_points: Number of points for plotting (default: 1000)
        save_csv: Whether to save results to CSV (default: True)
        species_labels: Custom labels for species in plots (default: None)
    """
    t_span: Tuple[float, float] = (0.0, 10.0)
    initial_concentrations: Optional[Dict[str, float]] = None
    solver_method: str = "BDF"
    atol: float = 1e-4
    plot: bool = True
    plot_species_indices: Optional[List[int]] = None
    plot_sample_points: int = 1000
    save_csv: bool = True
    species_labels: Optional[Dict[int, str]] = None


def calculate_ode_solution(
    complex_reaction_system: Any,
    config: Union[ODEPipelineConfig, Dict] = None,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Calculate ODE solution for a complex reaction system.
    
    Args:
        complex_reaction_system: The reaction system generated from PDB model
        config: Configuration for ODE calculation (ODEPipelineConfig or dict)
    
    Returns:
        Tuple of (time, concentrations, species_names)
        - time: 1D array of time points
        - concentrations: 2D array of shape (n_timepoints, n_species)
        - species_names: List of species names corresponding to concentration columns
    
    Example:
        >>> time, conc, species = calculate_ode_solution(reaction_system)
        >>> plt.plot(time, conc[:, 0], label=species[0])
    """
    # Handle config
    if config is None:
        config = ODEPipelineConfig()
    elif isinstance(config, dict):
        config = ODEPipelineConfig(**config)
    
    # Initialize parser
    rsp = ReactionStringParser()
    
    # Extract reaction information
    reaction_strings = [reaction.expression for reaction in complex_reaction_system.reactions]
    rate_constants = [reaction.rate for reaction in complex_reaction_system.reactions]
    
    # Parse reactions to get matrices
    species_names, rate_constant_names, reactant_matrix, product_matrix = \
        rsp.parse_reaction_strings(reaction_strings)
    
    # Setup initial concentrations
    n_species = len(species_names)
    y_init = np.zeros(n_species)
    
    if config.initial_concentrations:
        # Use user-provided initial concentrations
        for i, species in enumerate(species_names):
            val = config.initial_concentrations.get(species)
            if val is not None:
                y_init[i] = val
            elif '_' in species:
                # Try to map monomer names (e.g. A1_hash -> A)
                # Assumes species name format "Composition_Hash" where Composition is like "A1"
                composition = species.split('_')[0]
                # Check if it looks like a monomer "X1"
                if composition.endswith('1'):
                    # Extract type name (everything before '1')
                    mol_type = composition[:-1]
                    # Check if this type has an initial concentration
                    if mol_type in config.initial_concentrations and not any(c.isdigit() for c in mol_type):
                         # Ensure we didn't just strip a digit from "A12" -> "A1" (unlikely given sort, but safety)
                         # Actually just checking strict "X1" pattern matching against keys
                         val = config.initial_concentrations.get(mol_type)
                         if val is not None:
                             y_init[i] = val
    else:
        # Default: first species (monomer) at 1.0, rest at 0.0
        y_init[0] = 1.0
    
    # Solve ODE
    time, concentrations, species_names = solve_reaction_ode(
        dydt,
        config.t_span,
        y_init,
        reactant_matrix=reactant_matrix,
        product_matrix=product_matrix,
        k=rate_constants,
        plotting=False,  # We'll handle plotting separately
        method=config.solver_method,
        atol=config.atol,
        plotting_sample_points=config.plot_sample_points,
        species_names=species_names
    )
    
    return time, concentrations, species_names


    return saved_files

def _aggregate_by_composition(
    concentrations: np.ndarray,
    species_names: List[str]
) -> Tuple[np.ndarray, List[str]]:
    """
    Aggregate concentrations by species composition (prefix before '_').
    
    Args:
        concentrations: Array of shape (n_timepoints, n_species)
        species_names: List of species names
        
    Returns:
        Tuple of (aggregated_concentrations, unique_names)
    """
    # Find unique compositions and map indices
    composition_map = {} # name -> list of indices
    
    for i, name in enumerate(species_names):
        # Extract composition (e.g. "A2" from "A2_hash")
        # If no underscore, use full name
        if '_' in name:
            comp = name.split('_')[0]
        else:
            comp = name
            
        if comp not in composition_map:
            composition_map[comp] = []
        composition_map[comp].append(i)
    
    # Sort unique names for consistency
    unique_names = sorted(composition_map.keys())
    n_unique = len(unique_names)
    n_timepoints = concentrations.shape[0]
    
    agg_concentrations = np.zeros((n_timepoints, n_unique))
    
    for i, name in enumerate(unique_names):
        indices = composition_map[name]
        # Sum concentrations of all species with this composition
        agg_concentrations[:, i] = np.sum(concentrations[:, indices], axis=1)
        
    return agg_concentrations, unique_names

def save_ode_results(
    time: np.ndarray,
    concentrations: np.ndarray,
    species_names: List[str],
    output_dir: Path,
    config: ODEPipelineConfig = None,
    filename_prefix: str = "ode_results"
) -> Dict[str, Path]:
    """
    Save ODE results to files (CSV and optional plots).
    Also saves simplified results aggregated by composition.
    
    Args:
        time: Time points array
        concentrations: Concentration array (n_timepoints, n_species)
        species_names: List of species names
        output_dir: Directory to save results
        config: ODE pipeline configuration
        filename_prefix: Prefix for output files
    
    Returns:
        Dictionary with paths to saved files
    """
    if config is None:
        config = ODEPipelineConfig()
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_files = {}
    
    # ---------------------------------------------------------
    # Helper to save CSV and Plot for a given dataset
    # ---------------------------------------------------------
    def _save_dataset(time_arr, conc_arr, names, prefix):
        prefix_files = {}
        
        # Save CSV
        if config.save_csv:
            csv_path = output_dir / f"{prefix}.csv"
            with open(csv_path, 'w', newline='') as f:
                writer = csv.writer(f)
                species_list = list(names) if hasattr(names, '__iter__') else names
                writer.writerow(['time'] + species_list)
                for i, t in enumerate(time_arr):
                    writer.writerow([t] + conc_arr[i, :].tolist())
            prefix_files['csv'] = csv_path
            print(f"ODE results saved to: {csv_path}")
            
        # Generate Plot
        if config.plot:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            # Simple plotting for aggregated results (plot all)
            # For original results, use config indices
            if prefix == filename_prefix and config.plot_species_indices is not None:
                indices = config.plot_species_indices
            else:
                indices = range(len(names))
            
            for idx in indices:
                if idx < len(names):
                    label = names[idx]
                    # Use config labels only for exact matches on original
                    if prefix == filename_prefix and config.species_labels:
                        label = config.species_labels.get(idx, label)
                    
                    ax.plot(time_arr, conc_arr[:, idx], label=label, linewidth=2)
            
            ax.set_xlabel('Time (s)', fontsize=12)
            ax.set_ylabel(r'Concentration $(\mu\mathrm{M})$', fontsize=12)
            title = 'ODE Solution: Complex Assembly Kinetics'
            if prefix != filename_prefix:
                title += ' (Simplified)'
            ax.set_title(title, fontsize=14)
            ax.legend(loc='best')
            ax.grid(True, alpha=0.3)
            
            plot_path = output_dir / f"{prefix}.png"
            fig.savefig(plot_path, dpi=300, bbox_inches='tight')
            prefix_files['plot'] = plot_path
            print(f"ODE plot saved to: {plot_path}")
            plt.close(fig)
            
        return prefix_files

    # ---------------------------------------------------------
    # 1. Save Original Results
    # ---------------------------------------------------------
    original_files = _save_dataset(time, concentrations, species_names, filename_prefix)
    saved_files.update(original_files)
    
    # ---------------------------------------------------------
    # 2. Save Simplified Results (Aggregated)
    # ---------------------------------------------------------
    # Check if aggregation is possible (look for underscores)
    has_underscore = any('_' in name for name in species_names)
    
    if has_underscore:
        agg_conc, agg_names = _aggregate_by_composition(concentrations, species_names)
        
        # Determine prefix for simplified files
        # If original is "ode_solution", simplified is "ode_solution_simple"
        simple_prefix = f"{filename_prefix}_simple"
        
        simple_files = _save_dataset(time, agg_conc, agg_names, simple_prefix)
        
        # Add to saved_files with distinct keys if needed, 
        # or just assume user knows where to look. 
        # API return usually expects just 'csv', 'plot'. 
        # We can add 'csv_simple', 'plot_simple'
        if 'csv' in simple_files:
            saved_files['csv_simple'] = simple_files['csv']
        if 'plot' in simple_files:
            saved_files['plot_simple'] = simple_files['plot']
            
    return saved_files


def run_ode_pipeline(
    complex_reaction_system: Any,
    output_dir: Path,
    config: Union[ODEPipelineConfig, Dict] = None,
    filename_prefix: str = "ode_results"
) -> Tuple[np.ndarray, np.ndarray, List[str], Dict[str, Path]]:
    """
    Run complete ODE pipeline: calculate and save results.
    
    This is the main convenience function that combines calculation and saving.
    
    Args:
        complex_reaction_system: The reaction system from PDB model
        output_dir: Directory to save results
        config: ODE pipeline configuration
        filename_prefix: Prefix for output files
    
    Returns:
        Tuple of (time, concentrations, species_names, saved_files)
    """
    # Calculate ODE solution
    time, concentrations, species_names = calculate_ode_solution(
        complex_reaction_system, config
    )
    
    # Handle config for saving
    if config is None:
        config = ODEPipelineConfig()
    elif isinstance(config, dict):
        config = ODEPipelineConfig(**config)
    
    # Save results
    saved_files = save_ode_results(
        time, concentrations, species_names,
        output_dir, config, filename_prefix
    )
    
    return time, concentrations, species_names, saved_files

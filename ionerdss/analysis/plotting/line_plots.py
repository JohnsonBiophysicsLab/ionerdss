"""
Line plot functions for the ionerdss package.
These functions are focused on creating line plots of various metrics over time or other variables.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Optional, Tuple, Dict, Any

from ..data.core import Data

# Import utils from processors
from ..data.processors.utils import (
    compute_average_assembly_size,
    eval_condition
)


def plot_line_speciescopy_vs_time(
    data: Data,
    save_dir: str,
    simulations_index: list,
    legend: list,
    user_file_name: str = None,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6),
):
    """
    Plot species copy number vs. time for selected simulations.

    Parameters:
        save_dir (str): The base directory where simulations are stored.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Species or groups of species to plot.
            - [['A(A1!1).A(A1!1)']] → plot 'A(A1!1).A(A1!1)'
            - [['A(A1!1).A(A1!1)'], ['A(A2!1).A(A2!1)']] → plot two species separately
            - [['A(A1!1).A(A1!1)', 'A(A2!1).A(A2!1)']] → plot their sum
        user_file_name(str): User defined file name for csv output
        show_type (str): Display mode, "both", "individuals", or "average".
        simulations_dir (list, optional): List of directories for each simulation.
        figure_size (tuple): Size of the plot figure. 
    """
    # Ensure the save path for processed data exists
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    # Get simulation directories
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # Read data using the data IO utility
    all_sim_data = data.get_copy_numbers_data(selected_dirs)
    
    # Filter out None values (failed reads)
    all_sim_data = [df for df in all_sim_data['dataframes'] if df is not None]

    if not all_sim_data:
        print("No valid simulation data found.")
        return
    
    # Align data to the shortest time series
    min_length = min(len(df) for df in all_sim_data)
    all_sim_data = [df.iloc[:min_length] for df in all_sim_data]

    # Compute average and standard deviation
    time_values = all_sim_data[0]["Time (s)"].values
    species_data = {}

    for species_list in legend:
        species_key = "+".join(species_list)
        values = np.array([df[species_list].sum(axis=1).values for df in all_sim_data])

        species_data[species_key] = {
            "mean": values.mean(axis=0),
            "std": values.std(axis=0),
            "raw": values
        }

    # Save processed data
    for species, data in species_data.items():
        # Check whether file name is too long
        if len(f"{species.replace('+', '_')}") > 100 and user_file_name is None:
            print("Error: Generated File Name Too Long. Define a file name using the plot_figure optional argument user_file_name.")
            return
        
        # Use user-defined name or auto-generate
        if user_file_name is None:
            save_path = os.path.join(plot_data_dir, f"{species.replace('+', '_')}.csv")
        else:
            save_path = os.path.join(plot_data_dir, f"{user_file_name}.csv")
        
        df_to_save = pd.DataFrame({
            "Time (s)": time_values,
            "Mean": data["mean"],
            "Std": data["std"]
        })
        df_to_save.to_csv(save_path, index=False)
        print(f"Processed data for {species} saved to {save_path}")
    
    # Plot results
    plt.figure(figsize=figure_size)
    sns.set_style("ticks")
    
    for species, data in species_data.items():
        if show_type in {"individuals", "both"}:
            for i, sim_values in enumerate(data["raw"]):
                plt.plot(time_values, sim_values, alpha=0.3, linestyle="dashed", 
                         label=f"{species} (simulation {i})" if show_type == "both" else None)

        if show_type in {"average", "both"}:
            plt.plot(time_values, data["mean"], label=f"{species} (average)", linewidth=2)
            plt.fill_between(time_values, data["mean"] - data["std"], data["mean"] + data["std"], alpha=0.2)

    plt.xlabel("Time (s)")
    plt.ylabel("Copy Number")
    plt.legend()
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(plot_data_dir, "species_vs_time_plot.svg")
    plt.savefig(plot_path, format="svg")
    print(f"Plot saved to {plot_path}")
    plt.show()


def plot_line_maximum_assembly_size_vs_time(
    data:Data,
    save_dir: str,
    simulations_index: list,
    legend: list,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6),
):
    """
    Plot the maximum assembly size vs. time based on species composition in complexes.

    Parameters:
        save_dir (str): The base directory where simulations are stored.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Species to consider in assembly size calculation.
        show_type (str): Display mode, "both", "individuals", or "average".
        simulations_dir (list): List of simulation directories.
        figure_size (tuple): Size of the figure. 
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    all_sim_data = []
    
    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # Read histogram complex data for each simulation
    for sim_dir in selected_dirs:
        single_data = data.get_histogram_data(sim_dir)
        if not data["Time (s)"]:
            continue
            
        time_series = single_data["Time (s)"]
        max_assembly_sizes = []
        
        for complexes in single_data["complexes"]:
            max_size = max([sum(complex_dict.values()) for count, complex_dict in complexes], default=0)
            max_assembly_sizes.append(max_size)
            
        if time_series and max_assembly_sizes:
            df = pd.DataFrame({"Time (s)": time_series, "Max Assembly Size": max_assembly_sizes})
            all_sim_data.append(df)

    if not all_sim_data:
        print("No valid simulation data found.")
        return

    # Align data to the shortest time series
    min_length = min(len(df) for df in all_sim_data)
    all_sim_data = [df.iloc[:min_length] for df in all_sim_data]

    time_values = all_sim_data[0]["Time (s)"].values
    max_sizes = np.array([df["Max Assembly Size"].values for df in all_sim_data])

    avg_max_size = max_sizes.mean(axis=0)
    std_max_size = max_sizes.std(axis=0)

    # Save processed data
    save_path = os.path.join(plot_data_dir, "max_assembly_size_vs_time.csv")
    df_to_save = pd.DataFrame({
        "Time (s)": time_values,
        "Mean Max Assembly Size": avg_max_size,
        "Std Max Assembly Size": std_max_size
    })
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    # Plot the results
    plt.figure(figsize=figure_size)

    if show_type in {"individuals", "both"}:
        for i, sim_values in enumerate(max_sizes):
            plt.plot(time_values, sim_values, alpha=0.3, linestyle="dashed", 
                     label=f"Individual run {i}" if show_type == "both" else None)

    if show_type in {"average", "both"}:
        plt.plot(time_values, avg_max_size, label="Average", linewidth=2)
        plt.fill_between(time_values, avg_max_size - std_max_size, 
                        avg_max_size + std_max_size, alpha=0.2)

    plt.xlabel("Time (s)")
    plt.ylabel("Max Assembly Size")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "max_assembly_size_vs_time.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


def plot_line_average_assembly_size_vs_time(
    data:Data,
    save_dir: str,
    simulations_index: list,
    legend: list,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6),
):
    """
    Plot the average assembly size vs. time based on species composition in complexes.

    Parameters:
        save_dir (str): The base directory where simulations are stored.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Conditions for computing average assembly size.
        show_type (str): Display mode, "both", "individuals", or "average".
        simulations_dir (list): List of simulation directories.
        figure_size (tuple): Size of the figure. 
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    all_sim_data = []
    
    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # Read data for each simulation
    for sim_dir in selected_dirs:
        single_data = data.get_histogram_data(sim_dir)
        if not data["Time (s)"]:
            continue
            
        time_series = single_data["Time (s)"]
        condition_results = {condition: [] for condition in legend}
        
        for complexes in single_data["complexes"]:
            avg_sizes = compute_average_assembly_size(complexes, legend)
            for cond in legend:
                condition_results[cond].append(avg_sizes.get(cond, 0))
        
        if time_series:
            df = pd.DataFrame({"Time (s)": time_series, **condition_results})
            all_sim_data.append(df)

    if not all_sim_data:
        print("No valid simulation data found.")
        return

    # Align data to the shortest time series
    min_length = min(len(df) for df in all_sim_data)
    all_sim_data = [df.iloc[:min_length] for df in all_sim_data]

    time_values = all_sim_data[0]["Time (s)"].values
    avg_data = {cond: np.array([df[cond].values for df in all_sim_data]) for cond in legend}

    # Compute mean and standard deviation
    mean_values = {cond: data.mean(axis=0) for cond, data in avg_data.items()}
    std_values = {cond: data.std(axis=0) for cond, data in avg_data.items()}

    # Save processed data
    save_path = os.path.join(plot_data_dir, "average_assembly_size_vs_time.csv")
    df_to_save = pd.DataFrame({"Time (s)": time_values, 
                               **{f"Mean {cond}": mean_values[cond] for cond in legend},
                               **{f"Std {cond}": std_values[cond] for cond in legend}})
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    # Plot the data
    plt.figure(figsize=figure_size)
    sns.set_style("ticks")

    for cond in legend:
        if show_type in {"individuals", "both"}:
            for i, sim_values in enumerate(avg_data[cond]):
                plt.plot(time_values, sim_values, alpha=0.3, linestyle="dashed",
                         label=f"Individual run {i} ({cond})" if show_type == "both" else None)

        if show_type in {"average", "both"}:
            plt.plot(time_values, mean_values[cond], label=f"Average ({cond})", linewidth=2)
            plt.fill_between(time_values, mean_values[cond] - std_values[cond], 
                            mean_values[cond] + std_values[cond], alpha=0.2)

    plt.xlabel("Time (s)")
    plt.ylabel("Average Assembly Size")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "average_assembly_size_vs_time.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


def plot_line_fraction_of_monomers_assembled_vs_time(
    data:Data,
    save_dir: str,
    simulations_index: list,
    legend: list,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6),
):
    """
    Plot the fraction of monomers assembled in complex vs. time based on species composition in complexes.
    
    Parameters:
        save_dir (str): The base directory where simulations are stored.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Conditions for computing assembly fractions (e.g., ["A>=2"]).
        show_type (str): Display mode, "both", "individuals", or "average".
        simulations_dir (list): List of simulation directories.
        figure_size (tuple): Size of the figure. 
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    all_sim_data = []
    
    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # Read data for each simulation
    for sim_dir in selected_dirs:
        single_data = data.get_histogram_data(sim_dir)
        if not single_data["Time (s)"]:
            continue
            
        time_series = single_data["Time (s)"]
        fraction_results = {condition: [] for condition in legend}
        
        for complexes in single_data["complexes"]:
            for cond in legend:
                selected_counts = 0
                total_counts = 0

                for count, complex_dict in complexes:
                    matches, target_species = eval_condition(complex_dict, cond)

                    if matches:
                        selected_counts += count * complex_dict.get(target_species, 0)  # Sum only the target species count

                    if target_species in complex_dict:
                        total_counts += count * complex_dict[target_species]  # Sum only in complexes where species exists

                fraction = selected_counts / total_counts if total_counts > 0 else 0
                fraction_results[cond].append(fraction)
        
        if time_series:
            df = pd.DataFrame({"Time (s)": time_series, **fraction_results})
            all_sim_data.append(df)
    
    if not all_sim_data:
        print("No valid simulation data found.")
        return
    
    # Align data to the shortest time series
    min_length = min(len(df) for df in all_sim_data)
    all_sim_data = [df.iloc[:min_length] for df in all_sim_data]
    
    time_values = all_sim_data[0]["Time (s)"].values
    fraction_data = {cond: np.array([df[cond].values for df in all_sim_data]) for cond in legend}
    
    # Compute mean and standard deviation
    mean_values = {cond: data.mean(axis=0) for cond, data in fraction_data.items()}
    std_values = {cond: data.std(axis=0) for cond, data in fraction_data.items()}
    
    # Save processed data
    save_path = os.path.join(plot_data_dir, "fraction_of_monomers_assembled_vs_time.csv")
    df_to_save = pd.DataFrame({"Time (s)": time_values, 
                               **{f"Mean {cond}": mean_values[cond] for cond in legend},
                               **{f"Std {cond}": std_values[cond] for cond in legend}})
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")
    
    # Plot the data
    plt.figure(figsize=figure_size)
    sns.set_style("ticks")
    
    for cond in legend:
        if show_type in {"individuals", "both"}:
            for i, sim_values in enumerate(fraction_data[cond]):
                plt.plot(time_values, sim_values, alpha=0.3, linestyle="dashed",
                         label=f"Individual run {i} ({cond})" if show_type == "both" else None)
        
        if show_type in {"average", "both"}:
            plt.plot(time_values, mean_values[cond], label=f"Average ({cond})", linewidth=2)
            plt.fill_between(time_values, mean_values[cond] - std_values[cond], 
                            mean_values[cond] + std_values[cond], alpha=0.2)
    
    plt.xlabel("Time (s)")
    plt.ylabel("Fraction of Monomers Assembled")
    plt.legend()
    plt.tight_layout()
    
    plot_path = os.path.join(plot_data_dir, "fraction_of_monomers_assembled_vs_time.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


def plot_complex_count_vs_time(
    data:Data,
    save_dir: str,
    simulations_index: list,
    target_complexes: list,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6),
):
    """
    Plot the count of specific complexes vs. time.
    
    Parameters:
        save_dir (str): The base directory where simulations are stored.
        simulations_index (list): Indices of the simulations to include.
        target_complexes (list): List of complex specifications to track (e.g., ["A: 1.", "A: 4.", "A: 2."]).
        show_type (str): Display mode, "both", "individuals", or "average".
        simulations_dir (list): List of simulation directories.
        figure_size (tuple): Size of the figure.
    """
    
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    all_sim_data = []
    
    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # Read data for each simulation
    for sim_dir in selected_dirs:
        single_data = data.get_histogram_data(sim_dir)
        if not single_data["Time (s)"]:
            continue
            
        time_series = single_data["Time (s)"]
        complex_counts = {complex_type: [] for complex_type in target_complexes}
        
        for complexes in single_data["complexes"]:
            # Initialize counts for this time point
            current_counts = {complex_type: 0 for complex_type in target_complexes}
            
            # Parse the complexes at this time point
            for count, complex_dict in complexes:
                # Convert complex_dict to string format like "A: 4."
                complex_str = format_complex_dict(complex_dict)
                
                if complex_str in target_complexes:
                    current_counts[complex_str] += count
            
            # Append counts for this time point
            for complex_type in target_complexes:
                complex_counts[complex_type].append(current_counts[complex_type])
        
        if time_series:
            df = pd.DataFrame({"Time (s)": time_series, **complex_counts})
            all_sim_data.append(df)
    
    if not all_sim_data:
        print("No valid simulation data found.")
        return
    
    # Align data to the shortest time series
    min_length = min(len(df) for df in all_sim_data)
    all_sim_data = [df.iloc[:min_length] for df in all_sim_data]
    
    time_values = all_sim_data[0]["Time (s)"].values
    count_data = {complex_type: np.array([df[complex_type].values for df in all_sim_data]) 
                  for complex_type in target_complexes}
    
    # Compute mean and standard deviation
    mean_values = {complex_type: data.mean(axis=0) for complex_type, data in count_data.items()}
    std_values = {complex_type: data.std(axis=0) for complex_type, data in count_data.items()}
    
    # Save processed data
    save_path = os.path.join(plot_data_dir, "complex_count_vs_time.csv")
    df_to_save = pd.DataFrame({"Time (s)": time_values,
                               **{f"Mean {complex_type}": mean_values[complex_type] for complex_type in target_complexes},
                               **{f"Std {complex_type}": std_values[complex_type] for complex_type in target_complexes}})
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")
    
    # Plot the data
    plt.figure(figsize=figure_size, dpi=600)
    sns.set_style("ticks")
    
    # Define colors for different complex types
    colors = plt.cm.tab10(np.linspace(0, 1, len(target_complexes)))
    
    for i, complex_type in enumerate(target_complexes):
        color = colors[i]
        
        if show_type in {"individuals", "both"}:
            for j, sim_values in enumerate(count_data[complex_type]):
                plt.plot(time_values, sim_values, alpha=0.3, linestyle="dashed", color=color,
                         label=f"Individual run {j} ({complex_type})" if show_type == "individuals" else None)
        
        if show_type in {"average", "both"}:
            plt.plot(time_values, mean_values[complex_type], label=f"Average ({complex_type})", 
                    linewidth=2, color=color)
            plt.fill_between(time_values, mean_values[complex_type] - std_values[complex_type], 
                            mean_values[complex_type] + std_values[complex_type], alpha=0.2, color=color)
    
    plt.xlabel("Time (s)")
    plt.ylabel("Complex Count")
    plt.legend()
    plt.tight_layout()
    
    plot_path = os.path.join(plot_data_dir, "complex_count_vs_time.svg")
    plt.savefig(plot_path, format="svg")
    # plt.show()
    print(f"Plot saved to {plot_path}")


def format_complex_dict(complex_dict):
    """
    Convert a complex dictionary to string format like "A: 4."
    
    Parameters:
        complex_dict (dict): Dictionary representing a complex (e.g., {'A': 4})
        
    Returns:
        str: Formatted string representation (e.g., "A: 4.")
    """
    if len(complex_dict) == 1:
        species, count = next(iter(complex_dict.items()))
        return f"{species}: {count}."
    else:
        # For multi-species complexes, sort by species name for consistency
        sorted_items = sorted(complex_dict.items())
        parts = [f"{species}: {count}" for species, count in sorted_items]
        return ", ".join(parts) + "."


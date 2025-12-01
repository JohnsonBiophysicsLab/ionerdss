"""
Histogram plot functions for the ionerdss package.
These functions create various types of histograms for analyzing species distributions.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from typing import List, Optional, Tuple, Dict, Any

from ..data.core import Data


from ..data.processors.utils import (
    filter_by_time_frame,
)

def plot_hist_complex_species_size(
    data:Data,
    save_dir: str,
    simulations_index: list,
    legend: list,
    bins: int = 10,
    time_frame: tuple = None,
    frequency: bool = False,
    normalize: bool = False,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6),
):
    """
    Plot a histogram of the average number or frequency of different complex species sizes over a time frame.

    Parameters:
        save_dir (str): The base directory where simulation results are stored.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Species to be counted in determining complex sizes.
        bins (int): Number of bins for the histogram.
        time_frame (tuple, optional): Time range (start, end) to consider for averaging.
        frequency (bool): Whether to plot frequency instead of absolute count.
        normalize (bool): Whether to normalize the histogram (ensuring total area = 1).
        show_type (str): Display mode - "both", "individuals", or "average".
        simulations_dir (list): List of directories for each simulation.
        figure_size (tuple): Size of the figure. 
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # Read data for each simulation
    all_sizes_per_sim = data.get_complex_sizes(selected_dirs)
    all_sizes_combined = []
    for sim_sizes in all_sizes_per_sim:
        all_sizes_combined.extend(sim_sizes)

    if not all_sizes_per_sim:
        print("No valid simulation data found.")
        return
    
    # Determine global bin edges
    global_hist, bin_edges = np.histogram(all_sizes_combined, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = bin_edges[1] - bin_edges[0]
    
    # Compute histograms for each simulation using the same bin edges
    hist_values_all = []

    for sizes in all_sizes_per_sim:
        hist_values, _ = np.histogram(sizes, bins=bin_edges)
        hist_values_all.append(hist_values)

    hist_values_all = np.array(hist_values_all)

    # Compute mean and standard deviation
    mean_values = np.mean(hist_values_all, axis=0)
    std_values = np.std(hist_values_all, axis=0)

    total = np.sum(mean_values)
    
    if frequency and total > 0:
        mean_values = mean_values / total
        std_values = std_values / total
    
    if normalize and total > 0:
        mean_values = mean_values / bin_width
        std_values = std_values / bin_width

    # Save data
    df_to_save = pd.DataFrame({
        "Bin Center": bin_centers,
        "Mean Count": mean_values,
        "Std Dev": std_values
    })
    save_path = os.path.join(plot_data_dir, "hist_average_number_vs_size.csv")
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    # Plot with error bars
    plt.figure(figsize=figure_size)
    plt.bar(bin_centers, mean_values, width=bin_width * 0.9, alpha=0.7, label="Mean")
    plt.errorbar(bin_centers, mean_values, yerr=std_values, fmt='o', color='black', capsize=5, label="Std Dev")

    species_all = "+".join(legend)
    plt.xlabel(f"Number of {species_all} in Complexes")
    plt.ylabel("Normalized Frequency" if normalize else ("Frequency" if frequency else "Complex Count"))
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "hist_average_number_vs_size.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


def plot_hist_monomer_counts_vs_complex_size(
    data:Data,
    save_dir: str,
    simulations_index: list,
    legend: list,
    bins: int = 10,
    time_frame: tuple = None,
    frequency: bool = False,
    normalize: bool = False,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6),
):
    """
    Plot a histogram of the total number of monomers as a function of complex size over a time frame.

    The X-axis represents complex species size (only considering species in the legend), 
    and the Y-axis represents the total number of monomers found in those complexes.

    Parameters:
        save_dir (str): Directory to save output plots.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Species to be counted in determining complex sizes.
        bins (int): Number of bins for the histogram.
        time_frame (tuple, optional): Time range (start, end) to consider for averaging.
        frequency (bool): Whether to plot frequency instead of absolute count.
        normalize (bool): Whether to normalize the histogram.
        show_type (str): Display mode - "both", "individuals", or "average".
        simulations_dir (list): List of directories for each simulation.
        figure_size (tuple): Size of the figure. 
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]

    # Step 1: Read data for each simulation
    all_sizes_per_sim = data.get_complex_sizes(selected_dirs)
    all_sizes_combined = []
    for sim_sizes in all_sizes_per_sim:
        all_sizes_combined.extend(sim_sizes)

    if not all_sizes_per_sim:
        print("No valid simulation data found.")
        return
    
    # Step 2: Determine global bin edges
    _, bin_edges = np.histogram(all_sizes_combined, bins=bins)  # Compute fixed bin edges
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = bin_edges[1] - bin_edges[0]

    # Step 3: Compute histograms for monomer counts using the same bin edges
    monomer_values_all = []

    for sizes in all_sizes_per_sim:
        monomer_values, _ = np.histogram(sizes, bins=bin_edges, weights=sizes)  # Weight = complex size
        monomer_values_all.append(monomer_values)

    monomer_values_all = np.array(monomer_values_all)

    # Step 4: Compute mean and standard deviation
    mean_values = np.mean(monomer_values_all, axis=0)
    std_values = np.std(monomer_values_all, axis=0)

    total = np.sum(mean_values)
    
    if frequency and total > 0:
        mean_values = mean_values / total
        std_values = std_values / total
    
    if normalize and total > 0:
        mean_values = mean_values / bin_width
        std_values = std_values / bin_width

    # Save data
    df_to_save = pd.DataFrame({
        "Bin Center": bin_centers,
        "Mean Monomer Count": mean_values,
        "Std Dev": std_values
    })
    save_path = os.path.join(plot_data_dir, "hist_monomer_count_vs_size.csv")
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    # Step 5: Plot with error bars
    plt.figure(figsize=figure_size)
    plt.bar(bin_centers, mean_values, width=bin_width * 0.9, alpha=0.7, label="Mean")
    plt.errorbar(bin_centers, mean_values, yerr=std_values, fmt='o', color='black', capsize=5, label="Std Dev")

    species_all = "+".join(legend)
    plt.xlabel(f"Number of {species_all} in Complexes")
    plt.ylabel("Normalized Frequency" if normalize else ("Frequency" if frequency else "Total Monomers in Complexes"))
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "hist_monomer_count_vs_size.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    
    print(f"Plot saved to {plot_path}")


def plot_stackedhist_complex_species_size(
    data:Data,
    save_dir: str,
    simulations_index: list,
    legend: list,
    bins: int = 10,
    time_frame: tuple = None,
    frequency: bool = False,
    normalize: bool = False,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6)
):
    """
    Plot a stacked histogram of complex species size over time.
    The X-axis represents complex species size (only considering species in the legend),
    and the Y-axis represents the number of complexes found with that size.
    Each color in the stack represents a different condition (species) from the legend.

    Parameters:
        save_dir (str): The base directory where simulation results are stored.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Species to be counted in determining complex sizes.
        bins (int): Number of bins for the histogram.
        time_frame (tuple, optional): Time range (start, end) to consider for averaging.
        frequency (bool): Whether to plot frequency instead of absolute count.
        normalize (bool): Whether to normalize the histogram.
        show_type (str): Display mode - "both", "individuals", or "average".
        simulations_dir (list): List of directories for each simulation.
        figure_size (tuple): Size of the figure.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)
    
    if not legend or ":" not in legend[0]:
        raise ValueError("Legend must be in format 'species:condition1,condition2,...'")
    
    x_species, y_conditions_str = legend[0].split(":")
    y_conditions = [cond.strip() for cond in y_conditions_str.split(",")]
    y_var = re.findall(r"[A-Za-z_]+", y_conditions[0])[0]

    all_histograms = []

    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # Read data from each simulation
    for sim_dir in selected_dirs:
        singledata = data.get_histogram_data(sim_dir)
        if not singledata["Time (s)"]:
            continue
            
        # Filter by time frame if specified
        if time_frame:
            singledata = filter_by_time_frame(singledata, time_frame)
            
        histogram = {cond: [] for cond in y_conditions}
        
        for complexes in singledata["complexes"]:
            for count, species_dict in complexes:
                x = species_dict.get(x_species, 0)
                y = species_dict.get(y_var, 0)
                
                for cond in y_conditions:
                    cond = cond.strip()
                    if any(op in cond for op in ['<', '>', '=', '<=', '>=']):
                        op_index = min(i for i, char in enumerate(cond) 
                                      if char in ['<', '>', '='])
                        cond_eval = cond[op_index:]
                        # replace = with == for eval
                        cond_eval = cond_eval.replace('=', '==') if '==' not in cond_eval else cond_eval
                        
                        # Evaluate the condition
                        if eval(f"{y}{cond_eval}"):
                            histogram[cond].extend([x] * count)
                            
        if any(histogram.values()):
            all_histograms.append(histogram)

    if not all_histograms:
        print("No valid simulation data found.")
        return

    # Calculate stacked counts for each condition
    stacked_counts = {cond: np.zeros(bins) for cond in y_conditions}
    all_data = []
    for hist in all_histograms:
        for cond in y_conditions:
            all_data.extend(hist.get(cond, []))
            
    if not all_data:
        print("No data to plot.")
        return
        
    bin_edges = np.histogram_bin_edges(all_data, bins=bins)

    for histogram in all_histograms:
        for cond, values in histogram.items():
            if values:
                hist, _ = np.histogram(values, bins=bin_edges)
                stacked_counts[cond] += hist

    # Average across simulations
    for cond in y_conditions:
        stacked_counts[cond] /= len(all_histograms)

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bottom = np.zeros_like(bin_centers)

    # Plot the stacked histogram
    plt.figure(figsize=figure_size)
    
    for cond in y_conditions:
        plt.bar(bin_centers, stacked_counts[cond], width=bin_edges[1]-bin_edges[0], 
                bottom=bottom, label=cond)
        bottom += stacked_counts[cond]

    plt.xlabel(f"Number of {x_species} in Complexes")
    plt.ylabel("Complex Count")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "stacked_hist_complex_species_size.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Stacked histogram saved to {plot_path}")

    # Save the data for further analysis
    stacked_df = pd.DataFrame(stacked_counts, index=bin_centers)
    stacked_df.to_csv(os.path.join(plot_data_dir, "stacked_hist_complex_species_size.csv"))
    print(f"Stacked histogram data saved to {os.path.join(plot_data_dir, 'stacked_hist_complex_species_size.csv')}")

def plot_hist_complex_species_size_3d(
    data:Data,
    save_dir: str,
    simulations_index: list,
    legend: list,
    bins: int = 10,
    time_bins: int = 10,
    frequency: bool = False,
    normalize: bool = False,
    simulations_dir: list = None,
    figure_size: tuple = (10, 8),
):
    """
    Plot a 3D histogram of complex species size over time.

    Parameters:
        save_dir (str): The base directory where simulation results are stored.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Species to be counted in determining complex sizes.
        bins (int): Number of bins for the histogram.
        time_bins (int): Number of time bins for the histogram.
        frequency (bool): Whether to plot frequency instead of absolute count.
        normalize (bool): Whether to normalize the histogram.
        simulations_dir (list): List of directories for each simulation.
        figure_size (tuple): Size of the figure. 
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    all_data = []

    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # First pass to collect all sizes and times
    for sim_dir in selected_dirs:
        single_data = data.get_histogram_data(sim_dir)
        if not single_data["Time (s)"]:
            continue
        
        sim_data = []
        for i, time in enumerate(single_data["Time (s)"]):
            for count, species_dict in single_data["complexes"][i]:
                size = sum(species_dict.get(s, 0) for s in legend if s in species_dict)
                sim_data.extend([(time, size)] * count)
                
        all_data.extend(sim_data)

    if not all_data:
        print("No valid data found.")
        return

    # Organize into time bins
    times, sizes = zip(*all_data)
    time_edges = np.linspace(min(times), max(times), time_bins + 1)
    size_edges = np.histogram_bin_edges(sizes, bins=bins)
    size_centers = (size_edges[:-1] + size_edges[1:]) / 2
    time_centers = (time_edges[:-1] + time_edges[1:]) / 2

    print(f"Time edges: {time_edges}")
    print(f"Time centers: {time_centers}")
    print(f"Size edges: {size_edges}")
    print(f"Size centers: {size_centers}")

    # Prepare 2D histogram: rows=time bins, cols=size bins
    hist2d = np.zeros((time_bins, bins))
    for t, s in all_data:
        t_idx = min(np.searchsorted(time_edges, t, side='right') - 1, time_bins - 1)
        s_idx = min(np.searchsorted(size_edges, s, side='right') - 1, bins - 1)
        if 0 <= t_idx < time_bins and 0 <= s_idx < bins:
            hist2d[t_idx, s_idx] += 1

    hist2d /= len(simulations_index)

    if frequency:
        with np.errstate(divide='ignore', invalid='ignore'):
            hist2d = hist2d / np.sum(hist2d, axis=1, keepdims=True)
            hist2d = np.nan_to_num(hist2d, nan=0.0)
        
    if normalize:
        bin_width = size_edges[1] - size_edges[0]
        hist2d = hist2d / bin_width

    # 3D bar plot
    fig = plt.figure(figsize=figure_size)
    ax = fig.add_subplot(111, projection='3d')

    xpos, ypos = np.meshgrid(size_centers, time_centers)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)
    dx = (size_edges[1] - size_edges[0]) * 0.9
    dy = (time_edges[1] - time_edges[0]) * 0.9
    dz = hist2d.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, shade=True)
    ax.set_xlabel("Complex Size")
    ax.set_ylabel("Time (s)")
    ax.set_zlabel("Normalized Frequency" if normalize else ("Frequency" if frequency else "Complex Count"))
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "3D_hist_complex_species.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"3D plot saved to {plot_path}")
    
    # Save the data for further analysis
    hist_data = pd.DataFrame(hist2d, index=time_centers, columns=size_centers)
    hist_data.to_csv(os.path.join(plot_data_dir, "hist_complex_species_size_3d.csv"))
    print(f"Histogram data saved to {os.path.join(plot_data_dir, 'hist_complex_species_size_3d.csv')}")


def plot_hist_monomer_counts_vs_complex_size_3d(
    data:Data,
    save_dir: str,
    simulations_index: list,
    legend: list,
    bins: int = 10,
    time_bins: int = 10,
    frequency: bool = False,
    normalize: bool = False,
    simulations_dir: list = None,
    figure_size: tuple = (10, 8)
):
    """
    Plot a 3D histogram of the total number of monomers as a function of complex size over time.
    The X-axis represents complex species size (only considering species in the legend),
    the Y-axis represents time intervals (in seconds),
    and the Z-axis represents the total number of monomers found in those complexes.

    Parameters:
        save_dir (str): The base directory where simulation results are stored.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Species to be counted in determining complex sizes.
        bins (int): Number of bins for the histogram.
        time_bins (int): Number of time bins for the histogram.
        frequency (bool): Whether to plot frequency instead of absolute count.
        normalize (bool): Whether to normalize the histogram.
        simulations_dir (list): List of directories for each simulation.
        figure_size (tuple): Size of the figure.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)
    
    all_data = []

    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # Read data from each simulation
    for sim_dir in selected_dirs:
        single_data = data.get_histogram_data(sim_dir)
        if not data["Time (s)"]:
            continue
        
        sim_data = []
        for i, time in enumerate(single_data["Time (s)"]):
            for count, species_dict in single_data["complexes"][i]:
                size = sum(species_dict.get(s, 0) for s in legend if s in species_dict)
                # Weight by size (number of monomers)
                sim_data.append((time, size, count * size))
                
        all_data.extend(sim_data)

    if not all_data:
        print("No valid data found.")
        return

    times, sizes, weights = zip(*all_data)
    time_edges = np.linspace(min(times), max(times), time_bins + 1)
    size_edges = np.histogram_bin_edges(sizes, bins=bins)
    size_centers = (size_edges[:-1] + size_edges[1:]) / 2
    time_centers = (time_edges[:-1] + time_edges[1:]) / 2

    print(f"Time edges: {time_edges}")
    print(f"Time centers: {time_centers}")
    print(f"Size edges: {size_edges}")
    print(f"Size centers: {size_centers}")

    hist2d = np.zeros((time_bins, bins))
    for t, s, w in all_data:
        t_idx = min(np.searchsorted(time_edges, t, side='right') - 1, time_bins - 1)
        s_idx = min(np.searchsorted(size_edges, s, side='right') - 1, bins - 1)
        if 0 <= t_idx < time_bins and 0 <= s_idx < bins:
            hist2d[t_idx, s_idx] += w

    hist2d /= len(simulations_index)

    if frequency:
        with np.errstate(divide='ignore', invalid='ignore'):
            hist2d = hist2d / np.sum(hist2d, axis=1, keepdims=True)
            hist2d = np.nan_to_num(hist2d, nan=0.0)
            
    if normalize:
        bin_width = size_edges[1] - size_edges[0]
        hist2d = hist2d / bin_width

    # 3D bar plot
    fig = plt.figure(figsize=figure_size)
    ax = fig.add_subplot(111, projection='3d')

    xpos, ypos = np.meshgrid(size_centers, time_centers)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)
    dx = (size_edges[1] - size_edges[0]) * 0.9
    dy = (time_edges[1] - time_edges[0]) * 0.9
    dz = hist2d.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, shade=True)
    ax.set_xlabel("Complex Size")
    ax.set_ylabel("Time (s)")
    ax.set_zlabel("Monomer Count" if not frequency else "Frequency")
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "3D_hist_monomer_species.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"3D plot saved to {plot_path}")

    # Save the data for further analysis
    hist_data = pd.DataFrame(hist2d, index=time_centers, columns=size_centers)
    hist_data.to_csv(os.path.join(plot_data_dir, "hist_monomer_count_vs_size_3d.csv"))
    print(f"Histogram data saved to {os.path.join(plot_data_dir, 'hist_monomer_count_vs_size_3d.csv')}")

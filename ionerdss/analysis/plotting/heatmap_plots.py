"""
Heatmap plot functions for the ionerdss package.
These functions create heatmap visualizations for complex species analysis.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from typing import List, Optional, Tuple, Dict, Any

from ..data.core import Data



def format_sig(x, sig=3):
    """Format a number with specified significant digits."""
    return f"{x:.{sig}g}"


def plot_heatmap_complex_species_size(
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
    Plot a 2D heatmap of the average number of different complex species sizes over time.
    The X-axis represents complex species size (only considering species in the legend),
    the Y-axis represents time intervals (in seconds),
    and the color in each box indicates the average number of corresponding complexes at each time interval.

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
                sim_data.extend([(time, size)] * count)
                
        all_data.extend(sim_data)

    if not all_data:
        print("No valid data found.")
        return

    times, sizes = zip(*all_data)
    time_edges = np.linspace(min(times), max(times), time_bins + 1)
    size_edges = np.histogram_bin_edges(sizes, bins=bins)

    hist2d, _, _ = np.histogram2d(times, sizes, bins=[time_edges, size_edges])
    hist2d /= len(simulations_index)

    if frequency:
        with np.errstate(divide='ignore', invalid='ignore'):
            hist2d = hist2d / np.sum(hist2d, axis=1, keepdims=True)
            hist2d = np.nan_to_num(hist2d, nan=0.0)
            
    if normalize:
        bin_width = size_edges[1] - size_edges[0]
        hist2d = hist2d / bin_width

    # Create DataFrame for the heatmap
    df = pd.DataFrame(hist2d, 
                      index=(time_edges[:-1] + time_edges[1:]) / 2,
                      columns=(size_edges[:-1] + size_edges[1:]) / 2)

    save_path = os.path.join(plot_data_dir, "heatmap_complex_species_size.csv")
    df.to_csv(save_path)
    print(f"Heatmap data saved to {save_path}")

    plt.figure(figsize=figure_size)
    heatmap = sns.heatmap(df, cmap="viridis", cbar_kws={
        "label": "Normalized Frequency" if normalize else ("Frequency" if frequency else "Complex Count")
    })
    plt.xlabel("Complex Size")
    plt.ylabel("Time (s)")
    
    # Format tick labels with appropriate precision
    heatmap.set_xticklabels([format_sig(x, 3) for x in df.columns])
    heatmap.set_yticklabels([format_sig(y, 3) for y in df.index])
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "heatmap_complex_species_size.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Heatmap plot saved to {plot_path}")


def plot_heatmap_monomer_counts_vs_complex_size(
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
    Plot a 2D heatmap of the average number of monomers as a function of complex size over time.
    The X-axis represents complex species size (only considering species in the legend),
    the Y-axis represents time intervals (in seconds),
    and the color in each box indicates the average number of monomers found in those complexes.

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

    # Weighted 2D histogram
    hist2d, _, _ = np.histogram2d(times, sizes, bins=[time_edges, size_edges], weights=weights)
    hist2d /= len(simulations_index)

    if frequency:
        with np.errstate(divide='ignore', invalid='ignore'):
            hist2d = hist2d / np.sum(hist2d, axis=1, keepdims=True)
            hist2d = np.nan_to_num(hist2d, nan=0.0)
            
    if normalize:
        bin_width = size_edges[1] - size_edges[0]
        hist2d = hist2d / bin_width

    # Create DataFrame for the heatmap
    df = pd.DataFrame(hist2d, 
                      index=(time_edges[:-1] + time_edges[1:]) / 2,
                      columns=(size_edges[:-1] + size_edges[1:]) / 2)

    save_path = os.path.join(plot_data_dir, "heatmap_monomer_counts_vs_complex_size.csv")
    df.to_csv(save_path)
    print(f"Heatmap data saved to {save_path}")

    plt.figure(figsize=figure_size)
    heatmap = sns.heatmap(df, cmap="viridis", cbar_kws={
        "label": "Normalized Frequency" if normalize else ("Frequency" if frequency else "Monomer Count")
    })
    plt.xlabel("Complex Size")
    plt.ylabel("Time (s)")
    
    # Format tick labels with appropriate precision
    heatmap.set_xticklabels([format_sig(x, 3) for x in df.columns])
    heatmap.set_yticklabels([format_sig(y, 3) for y in df.index])
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "heatmap_monomer_counts_vs_complex_size.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Heatmap plot saved to {plot_path}")


def plot_heatmap_species_a_vs_species_b(
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
    Plot a 2D heatmap of the average number of two selected species (species_a and species_b) in complexes over time.
    The X-axis represents the number of species_a, the Y-axis represents the number of species_b,
    and the color in each box indicates the average number of complexes containing those species.

    Parameters:
        save_dir (str): The base directory where simulation results are stored.
        simulations_index (list): Indices of the simulations to include.
        legend (list): Species to be counted in determining complex sizes, e.g., ["A", "B"].
        bins (int): Number of bins for the histogram.
        time_bins (int): Number of time bins for the histogram.
        frequency (bool): Whether to plot frequency instead of absolute count.
        normalize (bool): Whether to normalize the histogram.
        simulations_dir (list): List of directories for each simulation.
        figure_size (tuple): Size of the figure.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)
    
    if len(legend) < 2:
        raise ValueError("At least two species must be specified in the legend.")
    
    species_x, species_y = legend[0], legend[1]
    all_data = []

    # Get the simulation directories to process
    selected_dirs = [simulations_dir[idx] for idx in simulations_index]
    
    # Read data from each simulation
    for sim_dir in selected_dirs:
        single_data = data.get_histogram_data(sim_dir)
        if not data["Time (s)"]:
            continue
        
        sim_data = []
        for complexes in single_data["complexes"]:
            for count, species_dict in complexes:
                x = species_dict.get(species_x, 0)
                y = species_dict.get(species_y, 0)
                sim_data.extend([(x, y)] * count)
                
        all_data.extend(sim_data)

    if not all_data:
        print("No valid data found.")
        return

    x_vals, y_vals = zip(*all_data)
    heatmap, xedges, yedges = np.histogram2d(x_vals, y_vals, bins=bins)
    heatmap /= len(simulations_index)

    print(f"X edges: {xedges}")
    print(f"Y edges: {yedges}")

    if frequency:
        with np.errstate(divide='ignore', invalid='ignore'):
            heatmap = heatmap / np.sum(heatmap)
            heatmap = np.nan_to_num(heatmap, nan=0.0)
            
    if normalize:
        bin_area = (xedges[1] - xedges[0]) * (yedges[1] - yedges[0])
        heatmap = heatmap / bin_area

    # Create DataFrame for the heatmap
    df = pd.DataFrame(heatmap,
                      index=(xedges[:-1] + xedges[1:]) / 2,
                      columns=(yedges[:-1] + yedges[1:]) / 2)
                      
    save_path = os.path.join(plot_data_dir, "heatmap_species_a_vs_b.csv")
    df.to_csv(save_path)
    print(f"Heatmap data saved to {save_path}")

    plt.figure(figsize=figure_size)
    heatmap = sns.heatmap(df, cmap="viridis", cbar_kws={
        "label": "Normalized Frequency" if normalize else ("Frequency" if frequency else "Complex Count")
    })
    plt.xlabel(species_y)
    plt.ylabel(species_x)
    
    # Format tick labels with appropriate precision
    heatmap.set_xticklabels([format_sig(x, 3) for x in df.columns])
    heatmap.set_yticklabels([format_sig(y, 3) for y in df.index])
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "heatmap_species_a_vs_b.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Heatmap plot saved to {plot_path}")

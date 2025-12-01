"""
3D plot functions for the ionerdss package.
This module contains 3D histogram and related visualization functions.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from typing import List, Optional, Tuple, Dict, Any

from ..data.core import Data



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
        if not single_data["Time (s)"]:
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

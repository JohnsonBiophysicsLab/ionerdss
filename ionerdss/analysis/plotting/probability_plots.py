"""
Probability plot functions for the ionerdss package.
These functions create probability visualizations for cluster dynamics analysis.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Optional, Tuple, Dict, Any

from ..data.core import Data

def plot_line_free_energy(
    data: Data,
    save_dir: str,
    simulations_index: list,
    time_frame: tuple = None,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6),
):
    """
    Plot the change in free energy over a selected time frame for different sizes of complexes.

    The x-axis represents the size of the complex, and the y-axis represents the free energy in units of KbT.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    # Get transition data
    transition_data = data.get_transition_data(time_frame=time_frame)
    matrices = transition_data['matrices']
    
    all_free_energies = []
    
    for matrix in matrices:
        if matrix is None:
            continue
            
        counts_per_size = matrix.sum(axis=1)
        total = counts_per_size.sum()
        
        if total > 0:
            probabilities = counts_per_size / total
            with np.errstate(divide='ignore'):
                free_energy = -np.log(probabilities)
                free_energy[np.isinf(free_energy)] = np.nan
            all_free_energies.append(free_energy)

    if not all_free_energies:
        print("No valid simulation data found.")
        return

    # Align data to the shortest array
    min_length = min(len(arr) for arr in all_free_energies)
    all_free_energies = [arr[:min_length] for arr in all_free_energies]

    sizes = np.arange(1, min_length + 1)
    free_energy_array = np.array(all_free_energies)
    avg_free_energy = np.nanmean(free_energy_array, axis=0)
    std_free_energy = np.nanstd(free_energy_array, axis=0)

    # Save processed data
    df_to_save = pd.DataFrame({
        "Cluster Size": sizes,
        "Mean Free Energy (kBT)": avg_free_energy,
        "Std Free Energy": std_free_energy
    })
    save_path = os.path.join(plot_data_dir, "free_energy_vs_size.csv")
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    # Plot the data
    plt.figure(figsize=figure_size)

    if show_type in {"individuals", "both"}:
        for i, fe in enumerate(free_energy_array):
            plt.plot(sizes, fe, alpha=0.3, linestyle="dashed", 
                    label=f"Individual run {i}" if show_type == "both" else None)

    if show_type in {"average", "both"}:
        plt.plot(sizes, avg_free_energy, label="Average", linewidth=2)
        plt.fill_between(sizes, avg_free_energy - std_free_energy, 
                        avg_free_energy + std_free_energy, alpha=0.2)

    plt.xlabel("Cluster Size (n)")
    plt.ylabel(r"$\mathrm{Free\ Energy}\ (k_\mathrm{B}T)$")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "free_energy_vs_size.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")

def plot_line_symmetric_association_probability(
    data: Data,
    save_dir: str,
    simulations_index: list,
    legend: list = None,
    time_frame: tuple = None,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6)
):
    """
    This line plot represents the probability of association between complexes of different sizes.
    Each event is counted symmetrically from both participating sizes.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)
    
    conds = []
    if legend is None:
        legend = ["associate size > 2", "associate size = 2", "associate size < 2"]
        conds = [">2", "==2", "<2"]
    else:
        for l in legend:
            if not l.startswith("associate size"):
                raise ValueError(f"Legend '{l}' must start with 'associate size'.")
            cond = l.replace("associate size", "").strip()
            if not any(op in cond for op in [">=", "<=", "==", "!=", ">", "<", "="]):
                raise ValueError(f"Legend condition '{cond}' must contain a valid comparison operator.")
            cond = cond.replace("=", "==") if "=" in cond and "==" not in cond else cond
            conds.append(cond)

    # Get transition data
    transition_data = data.get_transition_data(time_frame=time_frame)
    matrices = transition_data['matrices']
    
    all_assoc_probs = []

    for matrix in matrices:
        if matrix is None:
            continue
            
        max_size = matrix.shape[0]
        assoc_probs = [[] for _ in conds]

        for n in range(max_size - 1):
            assoc_counts = []
            for m in range(n + 1, max_size):
                pair_size = m - n
                count = matrix[m, n]
                if pair_size == n + 1:
                    count /= 2
                assoc_counts.append((pair_size, count))

            total_assoc = sum(c for _, c in assoc_counts)
            for j, cond in enumerate(conds):
                selected = [c for s, c in assoc_counts if eval(f"{s}{cond}")]
                assoc_probs[j].append(sum(selected) / total_assoc if total_assoc > 0 else np.nan)

        all_assoc_probs.append(np.array(assoc_probs))

    if not all_assoc_probs:
        print("No valid simulation data found.")
        return

    min_length = min(prob.shape[1] for prob in all_assoc_probs)
    all_assoc_probs = [prob[:, :min_length] for prob in all_assoc_probs]

    cluster_sizes = np.arange(1, min_length + 1)
    prob_array = np.array(all_assoc_probs)

    avg_probs = np.nanmean(prob_array, axis=0)
    std_probs = np.nanstd(prob_array, axis=0)

    df_to_save = {"Cluster Size": cluster_sizes}
    for i, label in enumerate(legend):
        df_to_save[f"{label} (avg)"] = avg_probs[i]
        df_to_save[f"{label} (std)"] = std_probs[i]
    df_to_save = pd.DataFrame(df_to_save)

    save_path = os.path.join(plot_data_dir, "symmetric_association_probability.csv")
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    plt.figure(figsize=figure_size)

    if show_type in {"individuals", "both"}:
        for sim_idx, sim_probs in enumerate(prob_array):
            for i, label in enumerate(legend):
                plt.plot(cluster_sizes, sim_probs[i], linestyle="dashed", alpha=0.3, 
                         label=f"{label} (run {sim_idx})" if show_type == "both" else None)

    if show_type in {"average", "both"}:
        for i, label in enumerate(legend):
            plt.plot(cluster_sizes, avg_probs[i], label=f"{label} (avg)", linewidth=2)
            plt.fill_between(cluster_sizes, avg_probs[i] - std_probs[i], 
                             avg_probs[i] + std_probs[i], alpha=0.2)

    plt.xlabel("Cluster Size (n)")
    plt.ylabel("Association Probability")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "symmetric_association_probability.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


def plot_line_asymmetric_association_probability(
    data: Data,
    save_dir: str,
    simulations_index: list,
    legend: list = None,
    time_frame: tuple = None,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6)
):
    """
    This line plot represents the probability of association between complexes of different sizes.
    Each event is counted asymmetrically from the larger participating size.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)
    
    conds = []
    if legend is None:
        legend = ["associate size > 2", "associate size = 2", "associate size < 2"]
        conds = [">2", "==2", "<2"]
    else:
        for l in legend:
            if not l.startswith("associate size"):
                raise ValueError(f"Legend '{l}' must start with 'associate size'.")
            cond = l.replace("associate size", "").strip()
            if not any(op in cond for op in [">=", "<=", "==", "!=", ">", "<", "="]):
                raise ValueError(f"Legend condition '{cond}' must contain a valid comparison operator.")
            cond = cond.replace("=", "==") if "=" in cond and "==" not in cond else cond
            conds.append(cond)

    # Get transition data
    transition_data = data.get_transition_data(time_frame=time_frame)
    matrices = transition_data['matrices']
    
    all_assoc_probs = []

    for matrix in matrices:
        if matrix is None:
            continue
            
        max_size = matrix.shape[0]
        assoc_probs = [[] for _ in conds]

        for n in range(max_size - 1):
            assoc_counts = []
            for m in range(n + 1, max_size):
                pair_size = m - n
                count = matrix[m, n]
                if pair_size == n + 1:
                    count /= 2
                if pair_size > n + 1:
                    break
                assoc_counts.append((pair_size, count))

            total_assoc = sum(c for _, c in assoc_counts)
            for j, cond in enumerate(conds):
                selected = [c for s, c in assoc_counts if eval(f"{s}{cond}")]
                assoc_probs[j].append(sum(selected) / total_assoc if total_assoc > 0 else np.nan)

        all_assoc_probs.append(np.array(assoc_probs))

    if not all_assoc_probs:
        print("No valid simulation data found.")
        return

    min_length = min(prob.shape[1] for prob in all_assoc_probs)
    all_assoc_probs = [prob[:, :min_length] for prob in all_assoc_probs]

    cluster_sizes = np.arange(1, min_length + 1)
    prob_array = np.array(all_assoc_probs)

    avg_probs = np.nanmean(prob_array, axis=0)
    std_probs = np.nanstd(prob_array, axis=0)

    df_to_save = {"Cluster Size": cluster_sizes}
    for i, label in enumerate(legend):
        df_to_save[f"{label} (avg)"] = avg_probs[i]
        df_to_save[f"{label} (std)"] = std_probs[i]
    df_to_save = pd.DataFrame(df_to_save)

    save_path = os.path.join(plot_data_dir, "asymmetric_association_probability.csv")
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    plt.figure(figsize=figure_size)

    if show_type in {"individuals", "both"}:
        for sim_idx, sim_probs in enumerate(prob_array):
            for i, label in enumerate(legend):
                plt.plot(cluster_sizes, sim_probs[i], linestyle="dashed", alpha=0.3, 
                         label=f"{label} (run {sim_idx})" if show_type == "both" else None)

    if show_type in {"average", "both"}:
        for i, label in enumerate(legend):
            plt.plot(cluster_sizes, avg_probs[i], label=f"{label} (avg)", linewidth=2)
            plt.fill_between(cluster_sizes, avg_probs[i] - std_probs[i], 
                             avg_probs[i] + std_probs[i], alpha=0.2)

    plt.xlabel("Cluster Size (n)")
    plt.ylabel("Association Probability")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "asymmetric_association_probability.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


def plot_line_symmetric_dissociation_probability(
    data: Data,
    save_dir: str,
    simulations_index: list,
    legend: list = None,
    time_frame: tuple = None,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6)
):
    """
    This line plot represents the probability of dissociation between complexes of different sizes.
    Each event is counted symmetrically for both generated sizes.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)
    
    conds = []
    if legend is None:
        legend = ["dissociate size > 2", "dissociate size = 2", "dissociate size < 2"]
        conds = [">2", "==2", "<2"]
    else:
        for l in legend:
            if not l.startswith("dissociate size"):
                raise ValueError(f"Legend '{l}' must start with 'dissociate size'.")
            cond = l.replace("dissociate size", "").strip()
            if not any(op in cond for op in [">=", "<=", "==", "!=", ">", "<", "="]):
                raise ValueError(f"Legend condition '{cond}' must contain a valid comparison operator.")
            cond = cond.replace("=", "==") if "=" in cond and "==" not in cond else cond
            conds.append(cond)

    # Get transition data
    transition_data = data.get_transition_data(time_frame=time_frame)
    matrices = transition_data['matrices']

    all_dissoc_probs = []

    for matrix in matrices:
        if matrix is None:
            continue
            
        max_size = matrix.shape[0]
        dissoc_probs = [[] for _ in conds]

        for n in range(1, max_size):
            dissoc_counts = []
            # loop from n-1 to 0 to statistically count dissociations
            for m in range(n - 1, -1, -1):
                pair_size = n - m
                count = matrix[m, n]
                if pair_size == m + 1:
                    count /= 2
                dissoc_counts.append((pair_size, count))

            total_dissoc = sum(c for _, c in dissoc_counts)
            for j, cond in enumerate(conds):
                selected = [c for s, c in dissoc_counts if eval(f"{s}{cond}")]
                dissoc_probs[j].append(sum(selected) / total_dissoc if total_dissoc > 0 else np.nan)

        all_dissoc_probs.append(np.array(dissoc_probs))

    if not all_dissoc_probs:
        print("No valid simulation data found.")
        return

    min_length = min(prob.shape[1] for prob in all_dissoc_probs)
    all_dissoc_probs = [prob[:, :min_length] for prob in all_dissoc_probs]

    cluster_sizes = np.arange(2, min_length + 2)  # start from size 2 since we are looking at dissociations
    prob_array = np.array(all_dissoc_probs)

    avg_probs = np.nanmean(prob_array, axis=0)
    std_probs = np.nanstd(prob_array, axis=0)

    df_to_save = {"Cluster Size": cluster_sizes}
    for i, label in enumerate(legend):
        df_to_save[f"{label} (avg)"] = avg_probs[i]
        df_to_save[f"{label} (std)"] = std_probs[i]
    df_to_save = pd.DataFrame(df_to_save)

    save_path = os.path.join(plot_data_dir, "symmetric_dissociation_probability.csv")
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    plt.figure(figsize=figure_size)

    if show_type in {"individuals", "both"}:
        for sim_idx, sim_probs in enumerate(prob_array):
            for i, label in enumerate(legend):
                plt.plot(cluster_sizes, sim_probs[i], linestyle="dashed", alpha=0.3, 
                         label=f"{label} (run {sim_idx})" if show_type == "both" else None)

    if show_type in {"average", "both"}:
        for i, label in enumerate(legend):
            plt.plot(cluster_sizes, avg_probs[i], label=f"{label} (avg)", linewidth=2)
            plt.fill_between(cluster_sizes, avg_probs[i] - std_probs[i], 
                             avg_probs[i] + std_probs[i], alpha=0.2)

    plt.xlabel("Cluster Size (n)")
    plt.ylabel("Dissociation Probability")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "symmetric_dissociation_probability.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


def plot_line_asymmetric_dissociation_probability(
    data: Data,
    save_dir: str,
    simulations_index: list,
    legend: list = None,
    time_frame: tuple = None,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6)
):
    """
    This line plot represents the probability of dissociation between complexes of different sizes.
    Each dissociation event is counted once.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)
    
    conds = []
    if legend is None:
        legend = ["dissociate size > 2", "dissociate size = 2", "dissociate size < 2"]
        conds = [">2", "==2", "<2"]
    else:
        for l in legend:
            if not l.startswith("dissociate size"):
                raise ValueError(f"Legend '{l}' must start with 'dissociate size'.")
            cond = l.replace("dissociate size", "").strip()
            if not any(op in cond for op in [">=", "<=", "==", "!=", ">", "<", "="]):
                raise ValueError(f"Legend condition '{cond}' must contain a valid comparison operator.")
            cond = cond.replace("=", "==") if "=" in cond and "==" not in cond else cond
            conds.append(cond)

    # Get transition data
    transition_data = data.get_transition_data(time_frame=time_frame)
    matrices = transition_data['matrices']

    all_dissoc_probs = []

    for matrix in matrices:
        if matrix is None:
            continue
            
        max_size = matrix.shape[0]
        dissoc_probs = [[] for _ in conds]

        for n in range(1, max_size):
            dissoc_counts = []
            # loop from n-1 to 0 to statistically count dissociations
            for m in range(n - 1, -1, -1):
                pair_size = n - m
                count = matrix[m, n]
                if pair_size == m + 1:
                    count /= 2
                if pair_size > m + 1:
                    break
                dissoc_counts.append((pair_size, count))

            total_dissoc = sum(c for _, c in dissoc_counts)
            for j, cond in enumerate(conds):
                selected = [c for s, c in dissoc_counts if eval(f"{s}{cond}")]
                dissoc_probs[j].append(sum(selected) / total_dissoc if total_dissoc > 0 else np.nan)

        all_dissoc_probs.append(np.array(dissoc_probs))

    if not all_dissoc_probs:
        print("No valid simulation data found.")
        return

    min_length = min(prob.shape[1] for prob in all_dissoc_probs)
    all_dissoc_probs = [prob[:, :min_length] for prob in all_dissoc_probs]

    cluster_sizes = np.arange(2, min_length + 2)  # start from size 2 since we are looking at dissociations
    prob_array = np.array(all_dissoc_probs)

    avg_probs = np.nanmean(prob_array, axis=0)
    std_probs = np.nanstd(prob_array, axis=0)

    df_to_save = {"Cluster Size": cluster_sizes}
    for i, label in enumerate(legend):
        df_to_save[f"{label} (avg)"] = avg_probs[i]
        df_to_save[f"{label} (std)"] = std_probs[i]
    df_to_save = pd.DataFrame(df_to_save)

    save_path = os.path.join(plot_data_dir, "asymmetric_dissociation_probability.csv")
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    plt.figure(figsize=figure_size)

    if show_type in {"individuals", "both"}:
        for sim_idx, sim_probs in enumerate(prob_array):
            for i, label in enumerate(legend):
                plt.plot(cluster_sizes, sim_probs[i], linestyle="dashed", alpha=0.3, 
                         label=f"{label} (run {sim_idx})" if show_type == "both" else None)

    if show_type in {"average", "both"}:
        for i, label in enumerate(legend):
            plt.plot(cluster_sizes, avg_probs[i], label=f"{label} (avg)", linewidth=2)
            plt.fill_between(cluster_sizes, avg_probs[i] - std_probs[i], 
                             avg_probs[i] + std_probs[i], alpha=0.2)

    plt.xlabel("Cluster Size (n)")
    plt.ylabel("Dissociation Probability")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "asymmetric_dissociation_probability.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


def plot_line_growth_probability(
    data: Data,
    save_dir: str,
    simulations_index: list,
    legend: list = None,
    time_frame: tuple = None,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6)
):
    """
    This line plot represents the probability of growth between complexes of different sizes.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    # Get transition data
    transition_data = data.get_transition_data(time_frame=time_frame)
    matrices = transition_data['matrices']
    
    all_growth_probs = []

    for matrix in matrices:
        if matrix is None:
            continue
            
        max_size = matrix.shape[0]
        growth_probs = []

        for n in range(0, max_size):
            dissoc_counts = []
            assoc_counts = []
            
            # Count dissociation events for size n
            for m in range(n - 1, -1, -1):
                pair_size = n - m
                count = matrix[m, n]
                if pair_size == m + 1:
                    count /= 2
                if pair_size > m + 1:
                    break
                dissoc_counts.append(count)

            # Count association events for size n
            for m in range(n + 1, max_size):
                pair_size = m - n
                count = matrix[m, n]
                if pair_size == n + 1:
                    count /= 2
                assoc_counts.append(count)

            total_dissoc = sum(dissoc_counts) if dissoc_counts else 0
            total_assoc = sum(assoc_counts) if assoc_counts else 0
            growth_probs.append(total_assoc / (total_dissoc + total_assoc) if (total_dissoc + total_assoc) > 0 else np.nan)

        all_growth_probs.append(np.array(growth_probs))

    if not all_growth_probs:
        print("No valid simulation data found.")
        return

    min_length = min(len(prob) for prob in all_growth_probs)
    all_growth_probs = [prob[:min_length] for prob in all_growth_probs]

    cluster_sizes = np.arange(1, min_length + 1)
    prob_array = np.vstack(all_growth_probs)

    avg_probs = np.nanmean(prob_array, axis=0)
    std_probs = np.nanstd(prob_array, axis=0)

    df_to_save = pd.DataFrame({
        "Cluster Size": cluster_sizes,
        "Growth Probability (avg)": avg_probs,
        "Growth Probability (std)": std_probs
    })

    save_path = os.path.join(plot_data_dir, "growth_probability.csv")
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    plt.figure(figsize=figure_size)

    if show_type in {"individuals", "both"}:
        for i, sim_probs in enumerate(prob_array):
            plt.plot(cluster_sizes, sim_probs, linestyle="dashed", alpha=0.3, 
                     label=f"Run {i}" if show_type == "both" else None)

    if show_type in {"average", "both"}:
        plt.plot(cluster_sizes, avg_probs, label="Average", linewidth=2)
        plt.fill_between(cluster_sizes, avg_probs - std_probs, avg_probs + std_probs, alpha=0.2)

    plt.xlabel("Cluster Size (n)")
    plt.ylabel("Growth Probability")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "growth_probability.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


def plot_line_liftime(
    data: Data,
    save_dir: str,
    simulations_index: list,
    legend: list = None,
    time_frame: tuple = None,
    show_type: str = "both",
    simulations_dir: list = None,
    figure_size: tuple = (10, 6)
):
    """
    This line plot represents the average lifetime between complexes of different sizes.
    """
    plot_data_dir = os.path.join(save_dir, "figure_plot_data")
    os.makedirs(plot_data_dir, exist_ok=True)

    # Get transition data
    transition_data = data.get_transition_data(time_frame=time_frame)
    lifetimes_list = transition_data['lifetimes']

    all_lifetime_arrays = []
    max_cluster_size = 0

    for lifetime in lifetimes_list:
        if not lifetime:
            continue
            
        sizes = sorted(lifetime.keys())
        max_cluster_size = max(max_cluster_size, max(sizes, default=0))

        avg_lifetimes = []
        for size in range(1, max_cluster_size + 1):
            lifetimes = lifetime.get(size, [])
            avg_lifetimes.append(np.mean(lifetimes) if lifetimes else np.nan)

        all_lifetime_arrays.append(avg_lifetimes)

    if not all_lifetime_arrays:
        print("No valid simulation lifetime data found.")
        return
    
    # Pad all arrays to same length
    for i in range(len(all_lifetime_arrays)):
        diff = max_cluster_size - len(all_lifetime_arrays[i])
        if diff > 0:
            all_lifetime_arrays[i].extend([np.nan] * diff)

    cluster_sizes = np.arange(1, max_cluster_size + 1)
    lifetime_array = np.array(all_lifetime_arrays)

    avg_lifetime = np.nanmean(lifetime_array, axis=0)
    std_lifetime = np.nanstd(lifetime_array, axis=0)

    df_to_save = pd.DataFrame({
        "Cluster Size": cluster_sizes,
        "Mean Lifetime (s)": avg_lifetime,
        "Std Lifetime": std_lifetime
    })

    save_path = os.path.join(plot_data_dir, "lifetime_vs_size.csv")
    df_to_save.to_csv(save_path, index=False)
    print(f"Processed data saved to {save_path}")

    # Plotting
    plt.figure(figsize=figure_size)

    if show_type in {"individuals", "both"}:
        for i, lt in enumerate(lifetime_array):
            plt.plot(cluster_sizes, lt, linestyle="dashed", alpha=0.3, 
                     label=f"Sim {i}" if show_type == "both" else None)

    if show_type in {"average", "both"}:
        plt.plot(cluster_sizes, avg_lifetime, label="Average", linewidth=2)
        plt.fill_between(cluster_sizes, avg_lifetime - std_lifetime, 
                         avg_lifetime + std_lifetime, alpha=0.2)

    plt.xlabel("Cluster Size (n)")
    plt.ylabel("Average Lifetime (s)")
    plt.legend()
    plt.tight_layout()

    plot_path = os.path.join(plot_data_dir, "lifetime_vs_size.svg")
    plt.savefig(plot_path, format="svg")
    plt.show()
    print(f"Plot saved to {plot_path}")


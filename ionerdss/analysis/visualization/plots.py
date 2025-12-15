"""
Plotting functions for ionerdss.analysis.

These functions accept pandas DataFrames or numpy arrays and return matplotlib Axes objects.
They are decoupled from file I/O.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Optional, List, Union
from matplotlib.axes import Axes

from .config import PlotStyle

def plot_free_energy(
    df: pd.DataFrame, 
    ax: Optional[Axes] = None, 
    label: Optional[str] = None,
    **kwargs
) -> Axes:
    """
    Plots the Free Energy profile.
    
    Args:
        df: DataFrame with 'size' and 'free_energy' columns.
        ax: Matplotlib Axes object. If None, creates new figure.
        label: Legend label.
        **kwargs: Arguments passed to plt.plot.
    
    Returns:
        ax: The Axes object.
    """
    if ax is None:
        fig, ax = plt.subplots()
    
    if df.empty:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
        return ax
    
    defaults = PlotStyle.get_default_kwargs()
    defaults.update(kwargs)
    
    ax.plot(df['size'], df['free_energy'], label=label, marker='o', **defaults)
    
    ax.set_xlabel("Cluster Size")
    ax.set_ylabel("Free Energy ($k_B T$)")
    ax.set_title("Free Energy Profile")
    
    if label:
        ax.legend()
        
    return ax


def plot_size_distribution(
    df: pd.DataFrame,
    ax: Optional[Axes] = None,
    log_scale: bool = True,
    label: Optional[str] = None,
    **kwargs
) -> Axes:
    """
    Plots the cluster size probability distribution.
    
    Args:
        df: DataFrame with 'size' and 'probability' columns.
        log_scale: If True, use log scale for Y axis.
    """
    if ax is None:
        fig, ax = plt.subplots()
    
    if df.empty:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
        return ax
        
    defaults = PlotStyle.get_default_kwargs()
    defaults.update(kwargs)
    
    ax.bar(df['size'], df['probability'], label=label, alpha=0.7, **kwargs)
    
    if log_scale:
        ax.set_yscale('log')
        
    ax.set_xlabel("Cluster Size")
    ax.set_ylabel("Probability")
    ax.set_title("Cluster Size Distribution")
    
    if label:
        ax.legend()
        
    return ax


def plot_growth_probabilities(
    df: pd.DataFrame,
    ax: Optional[Axes] = None,
    **kwargs
) -> Axes:
    """
    Plots Growth (Association) vs Shrinkage (Dissociation) probabilities.
    """
    if ax is None:
        fig, ax = plt.subplots()
    
    if df.empty:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
        return ax
        
    defaults = PlotStyle.get_default_kwargs()
    defaults.update(kwargs)
    
    ax.plot(df['size'], df['growth_prob'], label='Growth (Assoc)', marker='^', color='green', **defaults)
    ax.plot(df['size'], df['shrink_prob'], label='Shrinkage (Dissoc)', marker='v', color='red', **defaults)
    
    ax.set_xlabel("Cluster Size")
    ax.set_ylabel("Probability")
    ax.set_title("Growth vs Shrinkage Probabilities")
    ax.legend()
    ax.set_ylim(0, 1.0)
    
    return ax


def plot_heatmap(
    matrix: np.ndarray,
    ax: Optional[Axes] = None,
    log_scale: bool = True,
    cmap: str = "viridis",
    title: str = "Transition Matrix"
) -> Axes:
    """
    Plots a heatmap of the transition matrix.
    """
    if ax is None:
        fig, ax = plt.subplots()
    
    if matrix.size == 0:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center', transform=ax.transAxes)
        return ax
        
    if log_scale:
        # Log scale handling: add small epsilon or mask zeros
        data = np.log1p(matrix)
        cbar_label = "Log(Count + 1)"
    else:
        data = matrix
        cbar_label = "Count"
        
    sns.heatmap(data, ax=ax, cmap=cmap, cbar_kws={'label': cbar_label})
    
    ax.set_xlabel("From Size (Index)")
    ax.set_ylabel("To Size (Index)")
    ax.set_title(title)
    # Invert Y axis to match matrix convention (0 at top) if desired, 
    # but usually size 1 at bottom left is preferred for physics plots?
    # Heatmap usually puts (0,0) at top-left.
    # If index 0 is size 1, and we want size 1 at bottom:
    ax.invert_yaxis()
    
    return ax


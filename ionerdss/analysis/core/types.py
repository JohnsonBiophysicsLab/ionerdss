"""
Type definitions for ionerdss.analysis.

This module defines the core data structures used throughout the library.
It uses TypedDict and dataclasses to ensure type safety without incurring
runtime overhead.
"""

from typing import Dict, List, Optional, TypedDict, Union, Any
from dataclasses import dataclass
import numpy as np
import pandas as pd
from pathlib import Path

# Use numpy.typing if available, otherwise fallback to Any
try:
    from numpy.typing import NDArray
    ArrayLike = NDArray[Any]
except ImportError:
    ArrayLike = np.ndarray


@dataclass
class SimulationMetadata:
    """Metadata associated with a simulation."""
    id: str
    path: Path
    box_size: Optional[List[float]] = None
    volume: Optional[float] = None
    params: Optional[Dict[str, Any]] = None


class TransitionData(TypedDict):
    """
    Container for transition matrix data.
    
    Attributes:
        matrix (ArrayLike): Square matrix (NxN) where M[i, j] is the count
            of transitions from size j+1 to size i+1 (or however the data is structured).
        time (float): The timestamp for this snapshot.
    """
    matrix: ArrayLike
    time: float


class LifetimeData(TypedDict):
    """
    Container for lifetime statistics.
    
    Attributes:
        lifetimes (Dict[int, List[float]]): Mapping from cluster size to list of observed lifetimes.
        time (float): The timestamp for this snapshot.
    """
    lifetimes: Dict[int, List[float]]
    time: float


@dataclass
class SimulationData:
    """
    aggregated data for a single simulation.
    """
    transitions: List[TransitionData]
    lifetimes: List[LifetimeData]
    # Raw dataframes/lists
    copy_numbers: Optional[pd.DataFrame] = None
    complex_histograms: Optional[List[Dict[str, Any]]] = None

    # Optional lazy-loaded DataFrames for caching processed results
    df_free_energy: Optional[pd.DataFrame] = None
    df_size_dist: Optional[pd.DataFrame] = None


class PlotConfig(TypedDict, total=False):
    """Configuration for plotting."""
    title: str
    xlabel: str
    ylabel: str
    figsize: List[float]
    style: str
    cmap: str
    normalize: bool
    log_scale: bool




"""
Main API for ionerdss.analysis.

This module exposes the Analyzer class, which is the primary entry point for users.
"""

from pathlib import Path
from typing import List, Union, Optional, Tuple, Dict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from .core.simulation import Simulation
from .io.loader import DataLoader
from .processing import transitions
from .visualization import plots, config

class Plotter:
    """
    Namespace for plotting functions bound to an Analyzer instance.
    """
    def __init__(self, analyzer: 'Analyzer'):
        self.analyzer = analyzer
        config.PlotStyle.set_style()

    def free_energy(self, simulation_index: int = 0, ax=None, **kwargs):
        """Plots free energy for a specific simulation."""
        sim = self.analyzer.get_simulation(simulation_index)
        df = self.analyzer.compute_free_energy(sim)
        return plots.plot_free_energy(df, ax=ax, label=sim.id, **kwargs)

    def size_distribution(self, simulation_index: int = 0, ax=None, **kwargs):
        """Plots size distribution for a specific simulation."""
        sim = self.analyzer.get_simulation(simulation_index)
        df = self.analyzer.compute_size_distribution(sim)
        return plots.plot_size_distribution(df, ax=ax, label=sim.id, **kwargs)

    def transitions(self, simulation_index: int = 0, ax=None, **kwargs):
        """Plots growth/shrinkage probabilities."""
        sim = self.analyzer.get_simulation(simulation_index)
        matrix = sim.get_transition_matrix()
        df = transitions.compute_transition_probabilities(matrix)
        return plots.plot_growth_probabilities(df, ax=ax, **kwargs)

    def heatmap(self, simulation_index: int = 0, ax=None, **kwargs):
        """Plots transition matrix heatmap."""
        sim = self.analyzer.get_simulation(simulation_index)
        matrix = sim.get_transition_matrix()
        return plots.plot_heatmap(matrix, ax=ax, **kwargs)


class Analyzer:
    """
    Main analysis controller.
    
    Usage:
        analyzer = Analyzer("./my_data")
        analyzer.plot.free_energy()
    """
    
    def __init__(self, root_dir: Union[str, Path]):
        self.root_dir = Path(root_dir)
        self.loader = DataLoader(self.root_dir)
        self.simulations: List[Simulation] = self.loader.discover_simulations()
        self.plot = Plotter(self)
        
    def get_simulation(self, index_or_id: Union[int, str]) -> Simulation:
        """Retrieves a simulation by index or ID."""
        if isinstance(index_or_id, int):
            if 0 <= index_or_id < len(self.simulations):
                return self.simulations[index_or_id]
            raise IndexError(f"Simulation index {index_or_id} out of range.")
        
        for sim in self.simulations:
            if sim.id == index_or_id:
                return sim
        raise KeyError(f"Simulation ID {index_or_id} not found.")

    def load_simulations(self, simulations: Optional[List[Union[int, str]]] = None, time_frame: Optional[Tuple[float, float]] = None) -> List[Simulation]:
        """
        Compatibility method to retrieve simulations.
        """
        if simulations is None:
             return self.simulations
        
        results = []
        for sim_id in simulations:
             try:
                 results.append(self.get_simulation(sim_id))
             except (IndexError, KeyError):
                 continue
        return results

    def compute_size_distribution(self, sim: Simulation) -> pd.DataFrame:
        """Computes size distribution for a simulation."""
        matrix = sim.get_transition_matrix()
        return transitions.compute_size_distribution(matrix)

    def compute_free_energy(self, sim: Simulation, temperature: float = 1.0) -> pd.DataFrame:
        """Computes free energy for a simulation."""
        # Check cache
        if sim.data.df_free_energy is not None:
            return sim.data.df_free_energy
            
        df_dist = self.compute_size_distribution(sim)
        df_fe = transitions.compute_free_energy(df_dist, temperature)
        
        # Cache result (careful with mutability if we allow partial updates later)
        sim.data.df_free_energy = df_fe
        return df_fe



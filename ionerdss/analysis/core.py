"""
Refactored Analysis class for the ionerdss package.
Main interface with improved separation of concerns.
"""

import os
from typing import List, Optional, Tuple, Dict, Any

from .data import Data
from .plotting import PlotConfigure
from .legacy import LegacyPlotInterface


class Analyzer:
    """
    Main interface for analyzing and visualizing NERDSS simulation results.
    
    Provides both new modular API and legacy compatibility.
    """
    
    def __init__(self, save_dir: str = None, verbose:bool=False):
        if save_dir is None:
            save_dir = os.getcwd()
        elif save_dir.startswith("~"):
            save_dir = os.path.expanduser(save_dir)
        
        self.save_dir = os.path.abspath(save_dir)
        self._discover_simulations(verbose)
        self._setup_directories()
        
        # Initialize subsystems
        self.data = Data()
        self.plot = PlotConfigure(self.save_dir)
        self._legacy = LegacyPlotInterface(self)
        
        # Configure data reading by default with all possibile directories
        self.load_simulations()
    
    def _discover_simulations(self, verbose:bool=True):
        """Discover simulation directories containing DATA folders."""
        if os.path.exists(os.path.join(self.save_dir, "DATA")):
            self.simulation_dirs = [self.save_dir]
            if verbose: print("Detected a single simulation directory.")
        else:
            self.simulation_dirs = []
            for root, dirs, _ in os.walk(self.save_dir):
                if "DATA" in dirs:
                    self.simulation_dirs.append(root)
                    # Don't recurse into directories that contain DATA
                    dirs[:] = [d for d in dirs if d != "DATA"]
            if verbose: print(f"Detected {len(self.simulation_dirs)} simulation directories.")
        
        if not self.simulation_dirs:
            raise ValueError(f"No simulation directories found in {self.save_dir}")
    
    def _setup_directories(self):
        """Create necessary output directories."""
        self.plot_data_dir = os.path.join(self.save_dir, "figure_plot_data")
        os.makedirs(self.plot_data_dir, exist_ok=True)
    
    def load_simulations(self, 
                 simulations: Optional[List[int]] = None,
                 species: Optional[List[str]] = None,
                 time_frame: Optional[Tuple[float, float]] = None) -> Data:
        """
        Get a configured Data object for processing specific simulations.
        
        Returns a Data object pre-configured with simulation subset and filters.
        The returned object can process different data types independently.
        """
        if simulations is None:
            simulations = list(range(len(self.simulation_dirs)))
        
        # Validate simulation indices
        max_idx = len(self.simulation_dirs) - 1
        invalid = [s for s in simulations if s < 0 or s > max_idx]
        if invalid:
            raise ValueError(f"Invalid simulation indices: {invalid}. Valid range: 0-{max_idx}")
        
        # Create and configure data object
        data = Data()
        data.configure(
            simulation_dirs=self.simulation_dirs,
            simulations=simulations,
            species=species,
            time_frame=time_frame,
            cache_dir=self.plot_data_dir
        )
        self.data = data
        return data
    
    def set_plot(self, **kwargs) -> PlotConfigure:
        """
        Get a configured PlotConfigure object with custom settings.
        
        Returns a PlotConfigure object pre-configured with plotting parameters.
        """
        plot_config = PlotConfigure(self.save_dir)
        plot_config.configure(**kwargs)
        self.plot = plot_config
        return plot_config
    
    # Legacy compatibility methods
    def plot_figure(self, figure_type: str, **kwargs):
        """Legacy interface for backward compatibility."""
        return self._legacy.plot_figure(figure_type, **kwargs)
    
    def visualize_trajectory(self, **kwargs):
        """Legacy interface for trajectory visualization."""
        return self._legacy.visualize_trajectory(**kwargs)
    
    # Quick access methods for common operations
    def quick_plot(self, plot_type: str, **kwargs):
        """
        Quick plotting method that handles data processing and plotting in one call.
        
        Combines load_simulations() and plotting for simple use cases.
        """
        # Extract data parameters
        data_params = {
            'simulations': kwargs.pop('simulations', None),
            'species': kwargs.pop('species', None),
            'time_frame': kwargs.pop('time_frame', None)
        }
        
        # Extract plot parameters
        plot_params = {k: v for k, v in kwargs.items() 
                      if k in ['figure_size', 'style', 'font_size', 'save_format']}
        
        # Get configured objects
        data = self.load_simulations(**data_params)
        plot_config = self.set_plot(**plot_params)
        
        # Dispatch to appropriate plotting method
        plot_method = getattr(plot_config, plot_type, None)
        if plot_method is None:
            raise ValueError(f"Unknown plot type: {plot_type}")
        
        return plot_method(data=data, **kwargs)
    
    def get_simulation_info(self) -> Dict[str, Any]:
        """Get information about discovered simulations."""
        return {
            'save_dir': self.save_dir,
            'num_simulations': len(self.simulation_dirs),
            'simulation_dirs': self.simulation_dirs,
            'plot_data_dir': self.plot_data_dir
        }
    
    def clear_cache(self):
        """Clear all cached data."""
        self.data.clear_cache()
        print("All cached data cleared.")
    
    # Context manager support for batch operations
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clear_cache()

"""
Legacy compatibility interface for ionerdss analysis.
Bridges between new modular API and existing plot_figure interface.
"""

import matplotlib.pyplot as plt
from typing import List, Any, Dict

class LegacyPlotInterface:
    """
    Legacy compatibility bridge for existing plot_figure interface.
    """
    
    def __init__(self, analyzer):
        self.analyzer = analyzer
        
    def plot_figure(self, figure_type: str, **kwargs):
        """
        Legacy plot_figure method with full backward compatibility.
        """
        simulations = kwargs.get('simulations', [0])
        x = kwargs.get('x', 'time')
        y = kwargs.get('y', 'species')
        
        # Map to new API
        if figure_type == "line":
            if x == "size" and y == "free_energy":
                return self.analyzer.plot.free_energy(simulation_index=simulations[0], **kwargs)
            elif x == "size" and "probability" in y:
                return self.analyzer.plot.transitions(simulation_index=simulations[0], **kwargs)
                
        elif figure_type == "hist":
            if x == "size" and "count" in y:
                return self.analyzer.plot.size_distribution(simulation_index=simulations[0], **kwargs)
                
        elif figure_type == "heatmap":
            return self.analyzer.plot.heatmap(simulation_index=simulations[0], **kwargs)

        print(f"Warning: Legacy plot type {figure_type} x={x} y={y} not fully supported in refactored version.")
        return None


"""
ionerdss.nerdss_analysis package - Enhanced analysis tools for NERDSS simulations.

This package provides comprehensive tools for analyzing and visualizing NERDSS 
simulation results with both modern modular API and legacy compatibility.

Main Classes:
    Analysis: Primary interface for simulation analysis and visualization
    Data: Advanced data processing with specialized processors  
    PlotConfigure: Centralized plotting configuration and execution

Usage:
```python
# New modular API (recommended)
import ionerdss as ion
analysis = ion.Analyzer("/path/to/simulations")
data = analysis.load_simulations(simulations=[0,1,2], species=["A","B"])
plot = analysis.set_plot(figure_size=(12,8))
plot.line_speciescopy_vs_time(data=data, legend=[["A"], ["B"]])

# Legacy API (backward compatible)
analysis = ion.Analyzer("/path/to/simulations") 
analysis.plot_figure(figure_type="line", x="time", y="species", 
                    legend=[["A"], ["B"]], simulations=[0,1,2])
```
"""

# =====================================================================
# Development Note
# 2025.06.11 
# Mankun Sang
# 
# Histogram reading and Copynumber reading now have a unified 
# reading pipeline: Processor.read().
# Transition matrix reading is not modified.
# TODO: Based on Histogram reading and Copynumber reading, further 
# develop data reading for other output files including transition 
# matrix.
# TODO: The goal is to remove dependency on DataIO which is a legacy 
# from the previous refactoring. It violates the structure of analysis 
# module. The structure should be:
# \_
#   |_ __init__.py (handel import)
#   |_ analysis.py (basic Analysis class)
#   |_ plot_figures.py (basic plotting class)
#   |_
#   |_ data
#       |_ __init__.py
#       |_ core.py (data handling class)
#       |_ processors
#            |_ data reading modules, each output file has a module
#            |_ ...
#   |_ plotting
#       |_ plotting modules, each type of figure has a module
#       |_ ...
#   |_ legacy
#       |_ API for backward compatibility
#       |_ ...
# =====================================================================

# Main interface
from .core import Analyzer

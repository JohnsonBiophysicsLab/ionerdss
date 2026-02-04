"""
ionerdss.analysis: A modular library for analyzing biophysical simulations.
"""

from __future__ import annotations
from .api import Analyzer
from .core.types import SimulationMetadata, SimulationData
from .legacy.interface import LegacyPlotInterface

__all__ = ["Analyzer", "SimulationMetadata", "SimulationData", "LegacyPlotInterface"]


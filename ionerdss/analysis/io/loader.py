"""
Data loading utilities.

This module handles the discovery and instantiation of Simulation objects
from a root directory structure.
"""

import os
from pathlib import Path
from typing import List, Optional, Generator
import logging

from ..core.simulation import Simulation

logger = logging.getLogger(__name__)

class DataLoader:
    """
    Handles discovery and loading of simulations.
    """
    
    def __init__(self, root_dir: Path):
        self.root_dir = Path(root_dir)

    def discover_simulations(self, max_depth: int = 3) -> List[Simulation]:
        """
        Recursively finds simulation directories.
        A directory is considered a simulation if it contains a 'DATA' subdirectory.
        
        Args:
            max_depth (int): Maximum recursion depth.
            
        Returns:
            List[Simulation]: List of discovered Simulation objects.
        """
        sims = []
        root_level = len(self.root_dir.parts)
        
        for root, dirs, files in os.walk(self.root_dir):
            current_path = Path(root)
            current_depth = len(current_path.parts) - root_level
            
            if current_depth > max_depth:
                # Stop recursing deeper than max_depth
                dirs[:] = [] # Clear dirs to stop walking this branch
                continue
                
            if "DATA" in dirs:
                # Found a simulation directory
                sim = Simulation(current_path)
                sims.append(sim)
                # Don't recurse into DATA
                dirs.remove("DATA")
        
        logger.info(f"Discovered {len(sims)} simulations in {self.root_dir}")
        return sims


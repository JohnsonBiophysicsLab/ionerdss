"""
Simulation object encapsulation.

This module defines the Simulation class, which acts as the primary interface
for accessing data from a single simulation run.
"""

from pathlib import Path
from typing import Optional, Dict, List, Tuple
import logging
import pandas as pd
import numpy as np
import numpy.typing as npt
import re

from .types import SimulationMetadata, SimulationData
from ..io.parser import parse_transition_file, parse_copy_numbers, parse_complex_histogram

logger = logging.getLogger(__name__)

class Simulation:
    """
    Represents a single simulation run.
    
    Attributes:
        path (Path): Root directory of the simulation.
        id (str): Unique identifier (usually the directory name).
        data (SimulationData): The loaded data container.
    """
    
    def __init__(self, path: Path, sim_id: Optional[str] = None):
        self.path = Path(path)
        self.id = sim_id or self.path.name
        self._data: Optional[SimulationData] = None
        self._metadata: SimulationMetadata = SimulationMetadata(
            id=self.id,
            path=self.path
        )

    @property
    def data(self) -> SimulationData:
        """Lazy-loads the simulation data."""
        if self._data is None:
            self.load()
        return self._data # type: ignore

    def load(self) -> None:
        """Loads data from the simulation directory."""
        logger.info(f"Loading simulation: {self.id}")
        
        # Define path to expected data file
        # Try standard locations
        tran_file = self.path / "DATA" / "transition_matrix_time.dat"
        copy_file = self.path / "DATA" / "copy_numbers_time.dat"
        hist_file = self.path / "DATA" / "histogram_complexes_time.dat"

        if not tran_file.exists():
            logger.warning(f"Transition matrix file not found at {tran_file}")
            transitions, lifetimes = [], []
        else:
            transitions, lifetimes = parse_transition_file(tran_file)
        
        copy_numbers = parse_copy_numbers(copy_file)
        hist_times, hist_comps, hist_matrix = parse_complex_histogram(hist_file)

        self._data = SimulationData(
            transitions=transitions,
            lifetimes=lifetimes,
            copy_numbers=copy_numbers,
            hist_times=hist_times,
            hist_comps=hist_comps,
            hist_matrix=hist_matrix
        )
        logger.info(f"Loaded {len(transitions)} transition matrix time points for {self.id}")
        logger.info(f"Loaded {len(hist_times)} histogram complex time points for {self.id}")

    def get_transition_matrix(self, time_range: Optional[Tuple[float, float]] = None) -> np.ndarray:
        """
        Aggregates transition matrices within a time range.
        
        Args:
            time_range: (start, end) tuple. If None, uses all data.
            
        Returns:
            np.ndarray: Summed transition matrix.
        """
        if not self.data.transitions:
            return np.array([])
            
        filtered = self.data.transitions
        if time_range:
            start, end = time_range
            filtered = [t for t in filtered if start <= t["time"] <= end]
            
        if not filtered:
            return np.array([])

        # Stack matrices and sum along axis 0
        # We need to handle potentially different shapes by padding
        matrices = [t["matrix"] for t in filtered]
        
        max_rows = max(m.shape[0] for m in matrices)
        max_cols = max(m.shape[1] for m in matrices)
        
        sum_matrix = np.zeros((max_rows, max_cols), dtype=int)
        
        for m in matrices:
            rows, cols = m.shape
            sum_matrix[:rows, :cols] += m
            
        return sum_matrix

    def get_lifetimes(self, cluster_size: int) -> List[float]:
        """
        Retrieves all lifetimes for a specific cluster size.
        """
        all_lifetimes = []
        for record in self.data.lifetimes:
            if cluster_size in record["lifetimes"]:
                all_lifetimes.extend(record["lifetimes"][cluster_size])
        return all_lifetimes

    def get_time_series(self, species_name: Union[str, list[str], dict, list[dict]]) -> tuple[npt.NDArray[np.float64], dict[Union[str,dict], npt.NDArray[np.float64]]]:
        """
        Time series for counts of given species (compositions) as found in the histogram complexes file.
        Accepts a single compisition or a list of compositions for species_name.
        Composition can be specified by string (e.g. "A: 84. c1: 75. L: 84.") or by dictionary (e.g. {"A": 84, "c1": 75, "L": 84}).

        Returns:
            time (1D numpy array length N), counts (1D array if single species, or MxN numpy array where M=len(species_name) )
        """
        if not isinstance(species_name, list):
            species_name = [species_name]

        species_ts = np.zeros((len(species_name),len(self.data.hist_times)))
        
        for i,species in enumerate(species_name):
            if isinstance(species,str):
                # Parse composition: "C: 1. A: 1.", regex finds "Key: Value" pairs
                comp_dict = {}
                for match in re.finditer(r"([A-Za-z0-9]+):\s*(\d+)", species):
                    comp_dict[match.group(1)] = int(match.group(2))
                try:
                    comp_ind = self.data.hist_comps.index(comp_dict)
                    species_ts[i,:] = self.data.hist_matrix[:,comp_ind].toarray().ravel()
                except ValueError:
                    logger.error(f"Complex {species} not found in histogram data for simulation {self.id}")
            elif isinstance(species, dict):
                try:
                    comp_ind = self.data.hist_comps.index(species)
                    species_ts[i,:] = self.data.hist_matrix[:,comp_ind].toarray().ravel()
                except ValueError:
                    logger.error(f"Complex {species} not found in histogram data for simulation {self.id}")
            else:
                logger.error(f"Individual complex compositions must be either str or dict. Received {type(species)}")
        
        return self.data.hist_times, species_ts.squeeze()

    def get_largest_size_time_series(self, include=None, exclude=None, only_count_includes=False):
        """
        Time series of the size of the largest size complex as found from the histogram file.
        
        The "include" kwarg can be used to only look at complexes containing particular monomers.
        
        The "exclude" kwarg can be used to ignore complexes containing particular other monomers,
        
        e.g. get_largest_size_time_series(["A","B"],["C","D"]) will return the time series of the size of
        the largest complex containing both A and B that does not containing any C or D molecules.

        If only_count_includes=True, then the calculated "size" only counts the monomers of "include" species,
        i.e. a complex A:4.B:1.C:3. with include=["A","B"] would have size 5 if only_count_includes=True and size 8 otherwise.

        Returns:
                time (1D numpy array length N), size (1D array)
        """
        if include is not None and not isinstance(include,list):
            include = [include]
        
        if exclude is None:
            exclude = [] # do NOT put empty list as default parameter value; leads to difficult bugs
        elif not isinstance(exclude,list):
            exclude = [exclude]

        target_comps = list(self.data.hist_comps) # make a copy that we will modify

        if include is not None:
            for inc in include:
                # restrict to complexes containing desired monomer
                target_comps = [comp for comp in target_comps if inc in comp]

        # restrict to complexes which do not contain any of the excluded monomers
        for e in exclude:
            target_comps = [comp for comp in target_comps if e not in comp]
        
        target_indices = [self.data.hist_comps.index(comp) for comp in target_comps]
        
        if only_count_includes:
            target_sizes = np.array([sum([comp[monomer] for monomer in include]) for comp in target_comps])
        else:
            target_sizes = np.array([sum([comp[monomer] for monomer in comp]) for comp in target_comps])

        largest_size_ts = np.zeros(len(self.data.hist_times))
        for i,row in enumerate(self.data.hist_matrix):
            largest_size_ts[i] = np.amax(target_sizes[row[target_indices].toarray()>0])

        return self.data.hist_times, largest_size_ts

    def __repr__(self) -> str:
        return f"<Simulation id={self.id} path={self.path}>"



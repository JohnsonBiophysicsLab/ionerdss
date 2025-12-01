"""
Transition matrix and lifetime data processor.
"""

import os
import re
import numpy as np
from typing import List, Dict, Any, Tuple, Optional
from collections import defaultdict
import logging

# Configure logging
logger = logging.getLogger(__name__)

class TransitionProcessor:
    """
    Specialized processor for transition matrix and lifetime data.
    
    Handles reading, parsing, matrix operations, probability calculations,
    and lifetime statistics for cluster dynamics analysis.
    """
    
    def __init__(self):
        self._cache = {}
        self._selected_dirs = []

    def configure(self, selected_dirs: List[str]):
        self._selected_dirs = selected_dirs

    def read(self, selected_dirs: Optional[List[str]] = None, config: Dict[str, Any] = None) -> Tuple[List[Tuple[np.ndarray, Dict[int, List[float]]]], str]:
        """
        Read transition matrix and lifetime data from simulation directories.
        
        Parameters:
            selected_dirs (Optional[List[str]]): List of simulation directories.
            config (Dict[str, Any]): Configuration dictionary, e.g. {'time_frame': (start, end)}.
            
        Returns:
            Tuple[List[Tuple[np.ndarray, Dict[int, List[float]]]], str]:
                A tuple containing the list of (matrix, lifetime) tuples and the mode ('Single' or 'Multiple').
        """
        config = config or {"time_frame": None}
        
        # parse selected directories
        if not selected_dirs:
            if not self._selected_dirs:
                raise FileNotFoundError("No directory selected for reading.")
            selected_dirs = self._selected_dirs

        if isinstance(selected_dirs, list):
            return self.read_multiple(selected_dirs, config), 'Multiple'
        elif isinstance(selected_dirs, str):
            return self.read_single(selected_dirs, config), 'Single'
        else:
            raise TypeError(f"selected_dirs must be list or str, got {type(selected_dirs)}")

    def read_single(self, sim_dir: str, config: Dict[str, Any] = None) -> Tuple[np.ndarray, Dict[int, List[float]]]:
        """
        Read transition matrix and lifetime data from a single simulation directory.
        
        Parameters:
            sim_dir (str): Path to the simulation directory.
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            Tuple[np.ndarray, Dict[int, List[float]]]: Transition matrix and lifetimes.
        """
        config = config or {}
        time_frame = config.get('time_frame')
        
        file_path = os.path.join(sim_dir, "DATA", "transition_matrix_time.dat")
        
        if not os.path.exists(file_path):
            logger.warning(f"Transition matrix file not found: {file_path}")
            return np.array([]), {}
        
        try:
            matrix, lifetime = self._parse_transition_lifetime_data(file_path, time_frame)
            logger.debug(f"Successfully read transition matrix from {file_path}")
            return matrix, lifetime
        except Exception as e:
            logger.error(f"Error processing transition matrix from {file_path}: {e}")
            return np.array([]), {}

    def read_multiple(self, selected_dirs: List[str], config: Dict[str, Any] = None) -> List[Tuple[np.ndarray, Dict[int, List[float]]]]:
        """
        Read transition matrices from multiple simulations.
        
        Parameters:
            selected_dirs (List[str]): List of simulation directories.
            config (Dict[str, Any]): Configuration dictionary.
            
        Returns:
            List[Tuple[np.ndarray, Dict[int, List[float]]]]: List of (matrix, lifetime) tuples.
        """
        cache_key = f"transitions_{hash(tuple(sorted(selected_dirs)))}"
        if cache_key in self._cache:
             return self._cache[cache_key]

        results = []
        for sim_dir in selected_dirs:
            result = self.read_single(sim_dir, config)
            results.append(result)
        
        self._cache[cache_key] = results
        return results

    def _parse_transition_lifetime_data(self, file_path: str, time_frame: Optional[Tuple[float, float]] = None) -> Tuple[np.ndarray, Dict[int, List[float]]]:
        """
        Parse transition matrix and lifetime data from a file.
        
        Parameters:
            file_path (str): Path to the transition matrix file.
            time_frame (Optional[Tuple[float, float]]): Time range (start, end) to consider.
            
        Returns:
            Tuple[np.ndarray, Dict[int, List[float]]]: 
                A tuple containing the transition matrix and a dictionary of lifetimes per cluster size.
        """
        try:
            with open(file_path, "r") as f:
                content = f.read()
        except Exception as e:
            logger.error(f"Error reading transition matrix file {file_path}: {e}")
            return np.array([]), {}

        time_blocks = re.split(r"time:\s*", content)[1:]
        if not time_blocks:
            logger.warning(f"No time blocks found in {file_path}")
            return np.array([]), {}

        time_data = []

        for block in time_blocks:
            try:
                lines = block.strip().splitlines()
                if not lines:
                    continue
                    
                try:
                    time_val = float(lines[0])
                except ValueError:
                    continue
                
                # Parse transition matrix
                tm_lines = []
                tm_start = None
                
                for i, line in enumerate(lines):
                    if "transion matrix for each mol type:" in line:
                        tm_start = i + 2
                        break
                
                if tm_start is None:
                    continue
                    
                for i in range(tm_start, len(lines)):
                    if lines[i].startswith("lifetime for each mol type:"):
                        break
                    if lines[i].strip() and not lines[i].startswith(('A', 'B', 'C')):
                        try:
                            row = [int(x) for x in lines[i].split()]
                            if row:
                                tm_lines.append(row)
                        except ValueError:
                            continue
                
                if not tm_lines:
                    continue
                    
                transition_matrix = np.array(tm_lines)

                # Parse lifetimes
                lifetime = defaultdict(list)
                lt_start = None
                
                for i, line in enumerate(lines):
                    if "lifetime for each mol type:" in line:
                        lt_start = i + 2
                        break
                
                if lt_start is not None:
                    cluster_size = None
                    for line in lines[lt_start:]:
                        if line.startswith("size of the cluster:"):
                            try:
                                cluster_size = int(line.split(":")[1])
                            except (ValueError, IndexError):
                                continue
                        elif cluster_size is not None and line.strip():
                            try:
                                lifetimes = [float(x) for x in line.strip().split()]
                                lifetime[cluster_size].extend(lifetimes)
                            except ValueError:
                                continue

                time_data.append((time_val, transition_matrix, lifetime))
                
            except Exception as e:
                logger.warning(f"Error parsing time block in {file_path}: {e}")
                continue

        if not time_data:
            logger.warning(f"No valid time data found in {file_path}")
            return np.array([]), {}

        # Sort by time
        time_data.sort(key=lambda x: x[0])

        if time_frame:
            start, end = time_frame
            # Find the nearest time points within the frame
            valid_data = [(t, tm, lt) for t, tm, lt in time_data if start <= t <= end]
            
            if len(valid_data) >= 2:
                t_start, tm_start, lt_start = valid_data[0]
                t_end, tm_end, lt_end = valid_data[-1]
                
                # Calculate the difference in transition counts
                if tm_end.shape != tm_start.shape:
                    # Pad smaller matrix to match larger one
                    max_rows = max(tm_end.shape[0], tm_start.shape[0])
                    max_cols = max(tm_end.shape[1], tm_start.shape[1])
                    
                    tm_end_padded = np.zeros((max_rows, max_cols))
                    tm_end_padded[:tm_end.shape[0], :tm_end.shape[1]] = tm_end
                    
                    tm_start_padded = np.zeros((max_rows, max_cols))
                    tm_start_padded[:tm_start.shape[0], :tm_start.shape[1]] = tm_start
                    
                    matrix_delta = tm_end_padded - tm_start_padded
                else:
                    matrix_delta = tm_end - tm_start
                
                lifetime_delta = defaultdict(list)
                for k in lt_end:
                    lt1_len = len(lt_start.get(k, []))
                    lt2 = lt_end.get(k, [])
                    lifetime_delta[k] = lt2[lt1_len:]
            else:
                # Not enough data points in range
                return np.array([]), {}
        else:
            # If no time frame, use the last available data point (assuming cumulative)
            matrix_delta = time_data[-1][1]
            lifetime_delta = time_data[-1][2]

        return matrix_delta, dict(lifetime_delta)
    
    def aggregate_matrices(self, transition_data: Dict[str, Any]) -> np.ndarray:
        """Aggregate transition matrices across simulations."""
        matrices = transition_data['matrices']
        if not matrices:
            return np.array([])
        
        # Ensure all matrices have the same shape
        max_size = max(matrix.shape[0] for matrix in matrices)
        
        # Pad smaller matrices with zeros
        padded_matrices = []
        for matrix in matrices:
            if matrix.shape[0] < max_size:
                padded = np.zeros((max_size, max_size))
                padded[:matrix.shape[0], :matrix.shape[1]] = matrix
                padded_matrices.append(padded)
            else:
                padded_matrices.append(matrix)
        
        # Sum all matrices
        return np.sum(padded_matrices, axis=0)
    
    def calculate_size_probabilities(self, transition_data: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """Calculate probability distribution for each cluster size."""
        aggregated_matrix = self.aggregate_matrices(transition_data)
        if aggregated_matrix.size == 0:
            return {'probabilities': np.array([]), 'sizes': np.array([])}
        
        # Calculate total counts per size (sum over rows)
        counts_per_size = aggregated_matrix.sum(axis=1)
        total_counts = counts_per_size.sum()
        
        if total_counts == 0:
            probabilities = np.zeros_like(counts_per_size)
        else:
            probabilities = counts_per_size / total_counts
        
        sizes = np.arange(1, len(probabilities) + 1)
        
        return {
            'probabilities': probabilities,
            'sizes': sizes,
            'counts': counts_per_size,
            'total_counts': total_counts
        }
    
    def calculate_free_energy(self, transition_data: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """Calculate free energy landscape from size probabilities."""
        prob_data = self.calculate_size_probabilities(transition_data)
        probabilities = prob_data['probabilities']
        
        # Calculate free energy: F = -ln(P)
        with np.errstate(divide='ignore', invalid='ignore'):
            free_energy = -np.log(probabilities)
            free_energy[np.isinf(free_energy)] = np.nan
        
        return {
            'free_energy': free_energy,
            'sizes': prob_data['sizes'],
            'probabilities': probabilities
        }
    
    def calculate_association_probabilities(self, 
                                          transition_data: Dict[str, Any], 
                                          symmetric: bool = True) -> Dict[str, List[np.ndarray]]:
        """Calculate association probabilities for each cluster size."""
        aggregated_matrix = self.aggregate_matrices(transition_data)
        if aggregated_matrix.size == 0:
            return {'probabilities': [], 'sizes': np.array([])}
        
        max_size = aggregated_matrix.shape[0]
        association_probs = []
        
        for n in range(max_size - 1):
            assoc_counts = []
            
            for m in range(n + 1, max_size):
                pair_size = m - n
                count = aggregated_matrix[m, n]
                
                # For symmetric counting, divide by 2 when pair_size == n + 1
                if symmetric and pair_size == n + 1:
                    count /= 2
                    
                assoc_counts.append(count)
            
            total_assoc = sum(assoc_counts)
            if total_assoc > 0:
                probs = np.array(assoc_counts) / total_assoc
            else:
                probs = np.zeros(len(assoc_counts))
            
            association_probs.append(probs)
        
        return {
            'probabilities': association_probs,
            'sizes': np.arange(1, len(association_probs) + 1)
        }
    
    def calculate_dissociation_probabilities(self, 
                                           transition_data: Dict[str, Any], 
                                           symmetric: bool = True) -> Dict[str, List[np.ndarray]]:
        """Calculate dissociation probabilities for each cluster size."""
        aggregated_matrix = self.aggregate_matrices(transition_data)
        if aggregated_matrix.size == 0:
            return {'probabilities': [], 'sizes': np.array([])}
        
        max_size = aggregated_matrix.shape[0]
        dissociation_probs = []
        
        for n in range(1, max_size):
            dissoc_counts = []
            
            for m in range(n - 1, -1, -1):
                pair_size = n - m
                count = aggregated_matrix[m, n]
                
                # For symmetric counting, divide by 2 when pair_size == m + 1
                if symmetric and pair_size == m + 1:
                    count /= 2
                    
                dissoc_counts.append(count)
            
            total_dissoc = sum(dissoc_counts)
            if total_dissoc > 0:
                probs = np.array(dissoc_counts) / total_dissoc
            else:
                probs = np.zeros(len(dissoc_counts))
            
            dissociation_probs.append(probs)
        
        return {
            'probabilities': dissociation_probs,
            'sizes': np.arange(2, len(dissociation_probs) + 2)
        }
    
    def calculate_growth_probabilities(self, transition_data: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """Calculate growth vs shrinkage probabilities for each cluster size."""
        aggregated_matrix = self.aggregate_matrices(transition_data)
        if aggregated_matrix.size == 0:
            return {'growth_probs': np.array([]), 'sizes': np.array([])}
        
        max_size = aggregated_matrix.shape[0]
        growth_probs = []
        
        for n in range(max_size):
            # Count dissociation events (shrinkage)
            dissoc_counts = []
            for m in range(n - 1, -1, -1):
                pair_size = n - m
                count = aggregated_matrix[m, n]
                if pair_size == m + 1:
                    count /= 2
                dissoc_counts.append(count)
            
            # Count association events (growth)
            assoc_counts = []
            for m in range(n + 1, max_size):
                pair_size = m - n
                count = aggregated_matrix[m, n]
                if pair_size == n + 1:
                    count /= 2
                assoc_counts.append(count)
            
            total_dissoc = sum(dissoc_counts)
            total_assoc = sum(assoc_counts)
            total_events = total_dissoc + total_assoc
            
            if total_events > 0:
                growth_prob = total_assoc / total_events
            else:
                growth_prob = np.nan
            
            growth_probs.append(growth_prob)
        
        return {
            'growth_probs': np.array(growth_probs),
            'sizes': np.arange(1, len(growth_probs) + 1)
        }
    
    def aggregate_lifetimes(self, transition_data: Dict[str, Any]) -> Dict[int, List[float]]:
        """Aggregate lifetime data across simulations."""
        all_lifetimes = defaultdict(list)
        
        for lifetime_dict in transition_data['lifetimes']:
            for size, lifetimes in lifetime_dict.items():
                all_lifetimes[size].extend(lifetimes)
        
        return dict(all_lifetimes)
    
    def calculate_lifetime_statistics(self, transition_data: Dict[str, Any]) -> Dict[int, Dict[str, float]]:
        """Calculate statistical measures for lifetimes at each cluster size."""
        aggregated_lifetimes = self.aggregate_lifetimes(transition_data)
        
        lifetime_stats = {}
        for size, lifetimes in aggregated_lifetimes.items():
            if lifetimes:
                lifetime_array = np.array(lifetimes)
                lifetime_stats[size] = {
                    'mean': float(np.mean(lifetime_array)),
                    'std': float(np.std(lifetime_array)),
                    'median': float(np.median(lifetime_array)),
                    'min': float(np.min(lifetime_array)),
                    'max': float(np.max(lifetime_array)),
                    'count': len(lifetimes)
                }
            else:
                lifetime_stats[size] = {
                    'mean': 0.0, 'std': 0.0, 'median': 0.0,
                    'min': 0.0, 'max': 0.0, 'count': 0
                }
        
        return lifetime_stats
    
    def find_dominant_pathways(self, 
                             transition_data: Dict[str, Any], 
                             min_count: int = 10) -> List[Tuple[int, int, int]]:
        """Find dominant transition pathways (size_from, size_to, count)."""
        aggregated_matrix = self.aggregate_matrices(transition_data)
        if aggregated_matrix.size == 0:
            return []
        
        pathways = []
        rows, cols = np.where(aggregated_matrix >= min_count)
        
        for row, col in zip(rows, cols):
            count = int(aggregated_matrix[row, col])
            pathways.append((col + 1, row + 1, count))  # +1 for 1-based indexing
        
        # Sort by count (descending)
        pathways.sort(key=lambda x: x[2], reverse=True)
        
        return pathways
    
    def calculate_pathway_flux(self, transition_data: Dict[str, Any]) -> Dict[str, float]:
        """Calculate net flux for association vs dissociation pathways."""
        aggregated_matrix = self.aggregate_matrices(transition_data)
        if aggregated_matrix.size == 0:
            return {'association_flux': 0.0, 'dissociation_flux': 0.0, 'net_flux': 0.0}
        
        # Association flux: transitions to larger sizes
        assoc_flux = 0.0
        for i in range(aggregated_matrix.shape[0]):
            for j in range(i + 1, aggregated_matrix.shape[1]):
                assoc_flux += aggregated_matrix[j, i]
        
        # Dissociation flux: transitions to smaller sizes
        dissoc_flux = 0.0
        for i in range(aggregated_matrix.shape[0]):
            for j in range(i):
                dissoc_flux += aggregated_matrix[j, i]
        
        return {
            'association_flux': assoc_flux,
            'dissociation_flux': dissoc_flux,
            'net_flux': assoc_flux - dissoc_flux
        }
    
    def clear_cache(self):
        """Clear processor cache."""
        self._cache.clear()


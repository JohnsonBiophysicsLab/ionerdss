"""
Copy numbers data processor for time series analysis.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Any, Tuple, Optional
from .utils import determine_target_time_points, align_NERDSS
import os

# Configure logging, 
import logging
# inherit from the global level (should be setup in main)
logger = logging.getLogger(__name__)


class CopyNumberProcessor:
    """
    Specialized processor for species copy number time series data.
    
    Handles time series alignment, statistical analysis,
    and species group calculations.
    """
    
    def __init__(self):
        self._cache = {}
        self._selected_dirs = []

    def configure(self, selected_dirs: List[str]):
        self._selected_dirs = selected_dirs

    def read(self, selected_dirs, config = {"time_frame":None}) -> List[Dict[str, Any]]:
        """Decide to read multiple or read single"""

        # parse selected directories
        if not selected_dirs:
            if not self._selected_dirs:
                raise FileNotFoundError("No directory selected for reading.")
            selected_dirs = self._selected_dirs

        if isinstance(selected_dirs, list):
            return self.read_multiple(selected_dirs, config), 'Multiple'
        elif isinstance(selected_dirs, str):
            return self.read_single(selected_dirs), 'Single'
    
    def read_single(self, sim_dir: str) -> Dict[str, Any]:
        """
        Read species copy numbers vs time data from a simulation directory.
        
        Parameters:
            sim_dir (str): Path to the simulation directory
            
        Returns:
            Optional[pd.DataFrame]: DataFrame containing the data, or None if file not found
        """
        data_file = os.path.join(sim_dir, "DATA", "copy_numbers_time.dat")
        
        if not os.path.exists(data_file):
            logger.warning(f"Copy numbers file not found: {data_file}")
            return None
        
        try:
            df = pd.read_csv(data_file)
            df.rename(columns=lambda x: x.strip(), inplace=True)
            df = df.rename(columns={'time_points': 'Time (s)'}) # unify the name of time points
            logger.debug(f"Successfully read copy numbers from {data_file}")
            return df
        except Exception as e:
            logger.error(f"Error reading copy numbers from {data_file}: {e}")
            return None
    
    def read_multiple(
            self, 
            selected_dirs: Optional[List[str]] = None, 
            config:Dict[str,Any] = {"time_frame":None}, # left here for unified format
        ) -> List[Dict[str, Any]]:
        
        # check cache
        cache_key = "all_data"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # parse selected directories
        if not selected_dirs:
            if not self._selected_dirs:
                raise FileNotFoundError("No directory selected for reading.")
            selected_dirs = self._selected_dirs

        # Load raw data
        dataframes = []
        for sim_dir in selected_dirs:
            df = self.read_single(sim_dir)
            if df is not None:
                # TODO: filter by time frame. 
                # Now function filter_by_time_frame does not support copy number data
                # if config['time_frame']:
                #     df = filter_by_time_frame(df, config['time_frame'])
                dataframes.append(df)

        self._cache[cache_key] = dataframes

        return dataframes
    
    def calculate_species_groups(self, 
                               copy_data: Dict[str, Any], 
                               species_groups: List[List[str]]) -> Dict[str, np.ndarray]:
        """Calculate combined copy numbers for species groups, aligned to common time points."""

        dataframes = copy_data['dataframes']
        
        # Get target time points from copy_data (assumed to be present as 'common_time_points')
        # If not present, re-calculate
        if 'common_time_points' in copy_data:
             target_time = copy_data['common_time_points']
        elif 'aligned_data' in copy_data: # Backward compatibility if core wasn't updated yet
             target_time = copy_data['aligned_data']
        else:
             target_time = determine_target_time_points(dataframes)

        group_data = {}
        for i, group in enumerate(species_groups):
            group_name = f"group_{i}" if len(species_groups) > 1 else "combined"
            
            # Sum species in each group for each simulation AND align
            aligned_group_values = []
            for df in dataframes:
                available_species = [s for s in group if s in df.columns]
                if available_species:
                    group_sum = df[available_species].sum(axis=1).values
                else:
                    group_sum = np.zeros(len(df))
                
                # Align to target time points
                if 'Time (s)' in df.columns:
                    t_points = df['Time (s)'].values
                else:
                    # Fallback or error
                    t_points = np.arange(len(group_sum)) # Should not happen if read correctly
                
                aligned_sum, _ = align_NERDSS(t_points, group_sum, None, target_time)
                aligned_group_values.append(aligned_sum)
            
            group_data[group_name] = np.array(aligned_group_values)
        
        return group_data
    
    def compute_statistics(self, 
                         copy_data: Dict[str, Any], 
                         species_groups: List[List[str]]) -> Dict[str, Dict[str, np.ndarray]]:
        """Compute mean, std, min, max for species groups across simulations."""
        group_data = self.calculate_species_groups(copy_data, species_groups)
        
        stats = {}
        for group_name, values in group_data.items():
            stats[group_name] = {
                'mean': np.mean(values, axis=0),
                'std': np.std(values, axis=0),
                'min': np.min(values, axis=0),
                'max': np.max(values, axis=0),
                'median': np.median(values, axis=0),
                'raw_values': values
            }
        
        return stats
    
    def calculate_trends(self, 
                       copy_data: Dict[str, Any], 
                       species_groups: List[List[str]], 
                       window_size: int = 10) -> Dict[str, Dict[str, np.ndarray]]:
        """Calculate trends and derivatives for species groups."""
        group_data = self.calculate_species_groups(copy_data, species_groups)
        
        if 'common_time_points' in copy_data:
             time_points = copy_data['common_time_points']
        elif 'aligned_data' in copy_data:
             time_points = copy_data['aligned_data']
        else:
             time_points = determine_target_time_points(copy_data['dataframes'])
        
        trends = {}
        for group_name, values in group_data.items():
            mean_values = np.mean(values, axis=0)
            
            # Calculate moving average
            if len(mean_values) >= window_size:
                moving_avg = np.convolve(mean_values, np.ones(window_size)/window_size, mode='valid')
                # Pad to maintain original length
                pad_size = len(mean_values) - len(moving_avg)
                moving_avg = np.pad(moving_avg, (pad_size//2, pad_size - pad_size//2), mode='edge')
            else:
                moving_avg = mean_values
            
            # Calculate derivative (rate of change)
            if len(time_points) > 1:
                dt = np.diff(time_points)
                dy = np.diff(mean_values)
                derivative = np.append(dy / dt, dy[-1] / dt[-1])  # Extend last value
            else:
                derivative = np.zeros_like(mean_values)
            
            trends[group_name] = {
                'raw': mean_values,
                'smoothed': moving_avg,
                'derivative': derivative,
                'time_points': time_points
            }
        
        return trends
    
    def find_equilibrium_points(self, 
                              copy_data: Dict[str, Any], 
                              species_groups: List[List[str]], 
                              threshold: float = 0.01, 
                              min_duration: float = 10.0) -> Dict[str, List[Tuple[float, float]]]:
        """Find time periods where species reach equilibrium (low rate of change)."""
        trends = self.calculate_trends(copy_data, species_groups)
        
        equilibrium_periods = {}
        for group_name, trend_data in trends.items():
            time_points = trend_data['time_points']
            derivative = np.abs(trend_data['derivative'])
            
            # Find periods where derivative is below threshold
            below_threshold = derivative < threshold
            
            # Find continuous periods
            periods = []
            start_idx = None
            
            for i, is_below in enumerate(below_threshold):
                if is_below and start_idx is None:
                    start_idx = i
                elif not is_below and start_idx is not None:
                    duration = time_points[i-1] - time_points[start_idx]
                    if duration >= min_duration:
                        periods.append((time_points[start_idx], time_points[i-1]))
                    start_idx = None
            
            # Handle case where equilibrium extends to end
            if start_idx is not None:
                duration = time_points[-1] - time_points[start_idx]
                if duration >= min_duration:
                    periods.append((time_points[start_idx], time_points[-1]))
            
            equilibrium_periods[group_name] = periods
        
        return equilibrium_periods
    
    def calculate_time_to_equilibrium(self, 
                                    copy_data: Dict[str, Any], 
                                    species_groups: List[List[str]], 
                                    target_fraction: float = 0.95) -> Dict[str, float]:
        """Calculate time to reach target fraction of final value."""
        stats = self.compute_statistics(copy_data, species_groups)
        
        if 'common_time_points' in copy_data:
             time_points = copy_data['common_time_points']
        elif 'aligned_data' in copy_data:
             time_points = copy_data['aligned_data']
        else:
             time_points = determine_target_time_points(copy_data['dataframes'])
        
        equilibrium_times = {}
        for group_name, group_stats in stats.items():
            mean_values = group_stats['mean']
            final_value = mean_values[-1]
            target_value = final_value * target_fraction
            
            # Find first time point where target is reached
            reached_indices = np.where(mean_values >= target_value)[0]
            if len(reached_indices) > 0:
                equilibrium_times[group_name] = time_points[reached_indices[0]]
            else:
                equilibrium_times[group_name] = float('inf')  # Never reached
        
        return equilibrium_times
    
    def export_processed_data(self, 
                            copy_data: Dict[str, Any], 
                            species_groups: List[List[str]], 
                            output_path: str):
        """Export processed data to CSV files."""
        stats = self.compute_statistics(copy_data, species_groups)
        
        if 'common_time_points' in copy_data:
             time_points = copy_data['common_time_points']
        elif 'aligned_data' in copy_data:
             time_points = copy_data['aligned_data']
        else:
             time_points = determine_target_time_points(copy_data['dataframes'])
        
        for group_name, group_stats in stats.items():
            df = pd.DataFrame({
                'Time (s)': time_points,
                'Mean': group_stats['mean'],
                'Std': group_stats['std'],
                'Min': group_stats['min'],
                'Max': group_stats['max'],
                'Median': group_stats['median']
            })
            
            filename = f"{output_path}_{group_name}.csv"
            df.to_csv(filename, index=False)
    
    def clear_cache(self):
        """Clear processor cache."""
        self._cache.clear()

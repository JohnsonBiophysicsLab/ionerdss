"""
Histogram data processor for complex analysis.
"""
import os
import numpy as np
import re
import pandas as pd
from typing import List, Dict, Any, Tuple, Optional

# Helper functions
from .utils import parse_histogram_complex, parse_histogram_line, filter_by_time_frame, determine_target_time_points


# Configure logging, 
import logging
# inherit from the global level (should be setup in main)
logger = logging.getLogger(__name__)



class HistogramProcessor:
    """
    Specialized processor for histogram complex data.
    
    Handles complex size calculations, time series processing,
    and statistical analysis of complex distributions.
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
        Read histogram complex data from a simulation directory.
        
        Parameters:
            sim_dir (str): Path to the simulation directory
            
        Returns:
            Dict[str, Any]: Dictionary containing time series and complex data
        """
        data_file = os.path.join(sim_dir, "DATA", "histogram_complexes_time.dat")
        
        if not os.path.exists(data_file):
            logger.warning(f"Histogram complexes file not found: {data_file}")
            return {"time_series": [], "complexes": []}
        
        time_series = []
        all_complexes = []
        
        try:
            with open(data_file, "r") as f:
                lines = f.readlines()
            
            current_time = None
            current_complexes = []
            
            for line_num, line in enumerate(lines, 1):
                line = line.strip()
                if not line:
                    continue
                    
                time_match = re.match(r"Time \(s\):\s+([\d.]+(?:[eE][+-]?\d+)?)", line)
                if time_match:
                    if current_time is not None:
                        time_series.append(current_time)
                        all_complexes.append(current_complexes)
                        current_complexes = []
                    
                    current_time = float(time_match.group(1))
                else:
                    count, species_dict = parse_histogram_line(line)
                    if species_dict:
                        current_complexes.append((count, species_dict))
                    elif line.strip() and not line.startswith('#'):
                        logger.debug(f"Could not parse line {line_num}: {line}")
            
            # Add the last time point
            if current_time is not None and current_complexes:
                time_series.append(current_time)
                all_complexes.append(current_complexes)
            
            logger.debug(f"Successfully read histogram complexes from {data_file}")
            
        except Exception as e:
            logger.error(f"Error reading histogram complexes from {data_file}: {e}")
            return {"time_series": [], "complexes": []}
        
        return {
            "Time (s)": time_series,
            "complexes": all_complexes
        }

    def read_multiple(
            self, 
            selected_dirs: Optional[List[str]] = None, 
            config:Dict[str,Any] = {"time_frame":None}
        ) -> List[Dict[str, Any]]:
        """read multiple files"""

        # check cache
        cache_key = "all_data"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        if not selected_dirs:
            if not self._selected_dirs:
                raise FileNotFoundError("No directory selected for reading.")
            selected_dirs = self._selected_dirs

        # Load raw data
        all_data = []
        for sim_dir in selected_dirs:
            data = self.read_single(sim_dir)
            if data["Time (s)"]:
                if config['time_frame']:
                    data = filter_by_time_frame(data, config['time_frame'])
                all_data.append(data)

        self._cache[cache_key] = all_data

        return all_data

    def calculate_complex_sizes(self, 
                              histogram_data: Dict[str, Any], 
                              legend: List[str]) -> List[List[int]]:
        """Calculate complex sizes for specified species across all simulations."""
        cache_key = f"sizes_{hash(tuple(legend))}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        all_sizes = []
        for data in histogram_data['raw_data']:
            sim_sizes = []
        for complexes in data["complexes"]:
            for count, species_dict in complexes:
                complex_size = sum(species_dict[species] for species in legend if species in species_dict)
                sim_sizes.extend([complex_size] * count)
            all_sizes.append(sim_sizes)
        
        self._cache[cache_key] = all_sizes
        return all_sizes
    
    
    def get_size_distribution_stats(self, 
                                  histogram_data: Dict[str, Any], 
                                  legend: List[str]) -> Dict[str, float]:
        """Calculate statistical measures of size distribution."""
        all_sizes = self.calculate_complex_sizes(histogram_data, legend)
        combined_sizes = [size for sim_sizes in all_sizes for size in sim_sizes]
        
        if not combined_sizes:
            return {'mean': 0, 'std': 0, 'median': 0, 'max': 0, 'min': 0}
        
        sizes_array = np.array(combined_sizes)
        return {
            'mean': float(np.mean(sizes_array)),
            'std': float(np.std(sizes_array)),
            'median': float(np.median(sizes_array)),
            'max': int(np.max(sizes_array)),
            'min': int(np.min(sizes_array)),
            'total_complexes': len(combined_sizes)
        }
    
    def get_time_series(self, 
                        legends: List[List[str]],
                        legend_names: Optional[List[str]] = None
        ) -> pd.DataFrame:

        """
        Calculate time series statistics (mean, std, median) for different legends.
        
        This method processes histogram data across multiple simulations and time points
        to generate statistical measures for each legend over time.
        
        Parameters:
            histogram_data (Dict[str, Any]): Processed histogram data from get_histogram_data()
            legends (List[List[str]]): List of legend definitions, where each legend is a list of species.
                Example: ["A: 1.", "B: 2.", "A: 3. B: 4."] for separate A, B, and combined A+B analysis
            legend_names (Optional[List[str]]): Custom names for each legend. If None, use legends.
                Example: ["A1", "B2", "A3B4"]
        
        Returns:
            pd.DataFrame: DataFrame with columns:
                - Time (s): Time points
                - For each legend: Each trajectory's data for aligned time points
                
        Example:
            # Get time series statistics for individual and combined species
            legends = ["A: 1.", "B: 2.", "A: 3. B: 4."]
            legend_names = ["A1", "B2", "A3B4"]
            stats_df = processor.get_time_series(data, legends, legend_names)
        """
        
        cache_key = f"all_data"
        if cache_key in self._cache:
            all_data = self._cache[cache_key]
        else:
            logger.info("No data been read. Start reading...")
            try:
                all_data = self.read_multiple()
            except Exception as e:
                logger.error(f"No data provided: {e}")


        cache_key = f"time_series_{hash(tuple(sorted(legends)))}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Generate legend names if not provided
        if legend_names is None:
            legend_names = legends
        
        if len(legend_names) != len(legends):
            raise ValueError("Number of legend_names must match number of legends")
        
        legend_dicts = {}
        for lname, l in zip(legend_names, legends):
            legend_dicts[lname] = parse_histogram_complex(l)

        logger.debug('Dictionary of legends: ' + str(legend_dicts))
        
        common_time_points = determine_target_time_points(all_data)
        min_length = len(common_time_points)
        
        # Initialize result dictionary
        result_data = {'Time (s)': common_time_points}
        for lname in legend_dicts:
            result_data[lname] = [] # all data points
        
        # Process over all data
            
        for data in all_data:
            if not data['Time (s)'] or len(data['Time (s)']) < min_length:
                logger.critical(f"Illegal time_series. This should not happen. Min length calculated is {min_length}. \n" +
                                f"This time series has length {len(data['time_series'])}")
                continue
                
            # initialize
            time_series = {}
            for lname in legend_dicts:
                time_series[lname] = [0 for i in range(min_length)]
            
            logger.debug('Time series initialized with all 0: ' + str(time_series))

            for time_idx in range(min_length):
                complexes = data['complexes'][time_idx]
                for count, species_dict in complexes:
                    # find the corresponding legend name
                    for lname in legend_dicts:
                        if legend_dicts[lname] == species_dict:
                            time_series[lname][time_idx] = count

            logger.debug('Time series after reading: ' + str(time_series))
            
            for lname in legend_dicts:
                result_data[lname].append(time_series[lname])

        self._cache[cache_key] = result_data

        return result_data
        

    def calculate_time_series_statistics(
            self, 
            legends: List[List[str]],
            legend_names: Optional[List[str]] = None
        ) -> Dict[Dict,Any]:
        """
        Calculate time series statistics (mean, std, median) for different legends.
        
        This method processes histogram data across multiple simulations and time points
        to generate statistical measures for each legend over time.
        
        Parameters:
            histogram_data (Dict[str, Any]): Processed histogram data from get_histogram_data()
            legends (List[List[str]]): List of legend definitions, where each legend is a list of species.
                Example: ["A: 1.", "B: 2.", "A: 3. B: 4."] for separate A, B, and combined A+B analysis
            legend_names (Optional[List[str]]): Custom names for each legend. If None, use legends.
                Example: ["A1", "B2", "A3B4"]
        
        Returns:
            pd.DataFrame: DataFrame with columns:
                - Time (s): Time points
                - For each legend: {legend_name}_Mean, {legend_name}_Std, {legend_name}_Median
                
        Example:
            # Get time series statistics for individual and combined species
            legends = ["A: 1.", "B: 2.", "A: 3. B: 4."]
            legend_names = ["A1", "B2", "A3B4"]
            stats_df = processor.calculate_time_series_statistics(data, legends, legend_names)
            
            # Results in DataFrame with columns:
            # Time (s), Species_A_Mean, Species_A_Std, Species_A_Median, 
            # Species_B_Mean, Species_B_Std, Species_B_Median,
            # Combined_AB_Mean, Combined_AB_Std, Combined_AB_Median
        """
        cache_key = f"time_series_stats_{hash(tuple(sorted(legends)))}"
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        result_data = self.get_time_series(legends=legends, legend_names=legend_names)
            
        
        # Calculate statistics across simulations for each time point
        result_stat = {'Time (s)': result_data['Time (s)']}
        legend_names = [key for key in result_data if key != 'Time (s)']
        try:
            for species in legend_names:
                result_stat[species] = {
                    'mean':np.mean(result_data[species], axis=0), 
                    'std':np.std(result_data[species], axis=0), 
                    'median':np.median(result_data[species], axis=0),
                }
        except Exception as e:
            logger.error(f"Failed calculating stats: {e}")
        
        # Cache the result
        self._cache[cache_key] = result_stat
        
        return result_stat

    def clear_cache(self):
        """Clear processor cache."""
        self._cache.clear()



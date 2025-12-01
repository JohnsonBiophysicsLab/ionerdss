from typing import List, Dict, Any, Tuple, Optional
import re
import numpy as np
import pandas as pd

# Configure logging, 
import logging
# inherit from the global level (should be setup in main)
logger = logging.getLogger(__name__)

# ========================================================================
# Helper functions
# ========================================================================

def parse_histogram_complex(species_data_str: str) -> Optional[Dict[str, int]]:
    """
    Parse species data string and return a dictionary with species and counts.
    
    Parameters:
        species_data_str (str): Species data like "A: 2. B: 1."
        
    Returns:
        Optional[Dict[str, int]]: Dictionary mapping species to counts, or None if parsing fails
        
    Example:
        >>> parse_histogram_complex("A: 2. B: 1.")
        {'A': 2, 'B': 1}
        
        >>> parse_histogram_complex("A: 3. A: 1.")  # Duplicate species
        {'A': 4}
    """
    species_data = species_data_str.split()
    species_dict = {}

    try:
        for i in range(0, len(species_data), 2):
            species_name = species_data[i].strip(":")
            species_count = int(species_data[i + 1].strip("."))
            species_dict[species_name] = species_dict.get(species_name, 0) + species_count
    except (IndexError, ValueError) as e:
        logger.warning(f"Error parsing species data '{species_data_str}': {e}")
        return None

    return species_dict
     
def parse_histogram_line(line: str) -> Tuple[Optional[int], Optional[Dict[str, int]]]:
    """
    Parse a single complex line and return a dictionary with species and counts.
    
    Parameters:
        line (str): A line of text containing complex data
        
    Returns:
        Tuple[Optional[int], Optional[Dict[str, int]]]: 
            A tuple containing the count of complexes and a dictionary mapping species to counts
            
    Example:
        >>> parse_complex_line("5 A: 2. B: 1.")
        (5, {'A': 2, 'B': 1})
    """
    match = re.match(r"(\d+)\s+([\w\.\s:]+)", line)
    if not match:
        return None, None

    count = int(match.group(1))
    species_data_str = match.group(2)
    
    # Use the new parsing function
    species_dict = parse_histogram_complex(species_data_str)
    
    if species_dict is None:
        logger.warning(f"Error parsing complex line '{line}': failed to parse species data")
        return None, None

    return count, species_dict

def filter_by_time_frame(data: Dict[str, Any], time_frame: Tuple[float, float]) -> Dict[str, Any]:
        """Filter data by time frame."""
        start, end = time_frame
        filtered_indices = [
            i for i, t in enumerate(data["Time (s)"]) 
            if start <= t <= end
        ]
        
        return {
            "time_series": [data["Time (s)"][i] for i in filtered_indices],
            "complexes": [data["complexes"][i] for i in filtered_indices]
        }

def determine_target_time_points(all_data: List[Dict[str, Any]]) -> np.ndarray:
    """
    Find common time points across simulations and align datasets.
    
    This function determines a common set of time points (based on the shortest simulation
    by earliest end time) and returns these target time points.
    
    Args:
        all_data: List of data dictionaries or DataFrames. 
                 Each element must contain time information.
                 For histograms: {'Time (s)': [...]}
                 For copy numbers: DataFrame with 'Time' column
                 
    Returns:
        np.ndarray: Common time points array from the simulation with earliest end time.
    """
    if not all_data:
        return np.array([])

    # Extract time series from input data
    time_series_list = []
    
    for item in all_data:
        if isinstance(item, pd.DataFrame):
            # Handle Copy Number Data (DataFrame)
            if 'Time' in item.columns:
                t = item['Time'].values
                time_series_list.append(t)
            else:
                continue
        elif isinstance(item, dict) and 'Time (s)' in item:
            # Handle Histogram Data (Dict)
            time_series_list.append(np.array(item['Time (s)']))
        else:
            continue

    if not time_series_list:
        return np.array([])

    # Determine target time points (series with earliest end time)
    valid_series = [ts for ts in time_series_list if len(ts) > 0]
    if not valid_series:
        return np.array([])
        
    # Sort by end time (primary) and length (secondary) to ensure stability
    valid_series.sort(key=lambda x: (x[-1], len(x)))
    
    target_time_points = valid_series[0]
    
    return target_time_points

def align_NERDSS(t_points: np.ndarray, y_points: np.ndarray, err_points: Optional[np.ndarray], target_time_points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Aligns multiple NERDSS simulation datasets to a common set of time points.

    This function takes a list of simulation datasets, where each dataset is a tuple
    containing a time array and a corresponding data array. It uses linear interpolation
    to resample each dataset at the specified target time points.

    Args:
        t_points (array): time point arrays
        y_points (array): y arrays in format of [Y1, Y2, Y3, ...]
        err_points (array): err arrays in format of [E1, E2, E3, ...] (can be None)
        target_time_points (np.ndarray): An array of the desired time points to
                                         which all simulations should be aligned.

    Returns:
        aligned_y, aligned_err: corresponding to the target time series.
    """
    # copy number
    y_array = np.array(y_points)
    
    aligned_y = np.interp(
        target_time_points,
        t_points,
        y_array,
        left=y_array[0],
        right=y_array[-1]
    )

    # error
    if err_points is not None:
        err_array = np.array(err_points)
        aligned_err = np.interp(
            target_time_points,
            t_points,
            err_array,
            left=err_array[0],
            right=err_array[-1]
        )
    else:
        aligned_err = np.zeros_like(aligned_y)

    return aligned_y, aligned_err

def compute_average_assembly_size(complexes: List[Tuple[int, Dict[str, float]]], conditions: List[str]) -> Dict[str, float]:
    """
    Compute the average assembly size for given conditions.

    Parameters:
        complexes (list): List of tuples (count, species_dict) representing each complex.
        conditions (list): List of conditions, e.g., ["A>=2", "A+B>=4"].

    Returns:
        dict: Condition -> average assembly size mapping.
    """
    results = {}

    for condition in conditions:
        species_conditions = condition.split(", ")
        numerator, denominator = 0, 0

        for count, species_dict in complexes:
            valid = True
            total_size = 0

            for cond in species_conditions:
                species_match = re.match(r"(\w+)([>=<]=?|==)(\d+)", cond)
                if not species_match:
                    continue

                species, operator, threshold = species_match.groups()
                threshold = int(threshold)
                species_count = species_dict.get(species, 0)

                if operator == ">=" and species_count < threshold:
                    valid = False
                elif operator == ">" and species_count <= threshold:
                    valid = False
                elif operator == "<=" and species_count > threshold:
                    valid = False
                elif operator == "<" and species_count >= threshold:
                    valid = False
                elif operator == "==" and species_count != threshold:
                    valid = False

                total_size += species_count

            if valid:
                numerator += count * total_size
                denominator += count

        results[condition] = numerator / denominator if denominator > 0 else 0

    return results


def eval_condition(species_dict: Dict[str, float], condition: str) -> Tuple[bool, str]:
    """
    Evaluates whether a complex meets a condition based on species count.
    
    Parameters:
        species_dict (dict): Dictionary containing species counts in one complex.
        condition (str): A condition string like "B>=3".
    
    Returns:
        Tuple[bool, str]: (True if the complex satisfies the condition, species name)
    """
    species_match = re.match(r"(\w+)([>=<]=?|==)(\d+)", condition)
    if not species_match:
        return False, ""

    species, operator, threshold = species_match.groups()
    threshold = int(threshold)
    
    species_count = species_dict.get(species, 0)
    
    try:
        result = eval(f"{species_count} {operator} {threshold}")
        return result, species
    except:
        return False, species




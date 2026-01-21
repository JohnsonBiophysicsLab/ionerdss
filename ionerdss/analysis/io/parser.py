"""
File parsing logic for ionerdss.analysis.

This module handles reading legacy data files with robust regex-based parsing.
It replaces the brittle string splitting and eval() calls from the legacy code.
"""

import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.sparse import csc_array, lil_array

from ..core.types import TransitionData, LifetimeData

logger = logging.getLogger(__name__)

# Compiled regex patterns for performance
# Matches: "time: 0.123" or "Time (s): 0.123"
TIME_PATTERN = re.compile(r"(?:time|Time\s*\(s\)):\s*([\d\.]+)")
# Matches: "transition matrix for each mol type:"
TRANSITION_HEADER_PATTERN = re.compile(r"transition\s+matrix\s+for\s+each\s+mol\s+type:", re.IGNORECASE)
# Matches: "lifetime for each mol type:"
LIFETIME_HEADER_PATTERN = re.compile(r"lifetime\s+for\s+each\s+mol\s+type:", re.IGNORECASE)
# Matches: "size of the cluster: 5"
SIZE_HEADER_PATTERN = re.compile(r"size\s+of\s+the\s+cluster:\s*(\d+)", re.IGNORECASE)


def parse_transition_file(file_path: Path) -> Tuple[List[TransitionData], List[LifetimeData]]:
    """
    Parses the transition_matrix_time.dat file.

    Args:
        file_path (Path): Path to the .dat file.

    Returns:
        Tuple containing:
            - List[TransitionData]: Time-series of transition matrices.
            - List[LifetimeData]: Time-series of lifetime dictionaries.
    """
    if not file_path.exists():
        logger.warning(f"File not found: {file_path}")
        return [], []

    try:
        with open(file_path, 'r') as f:
            content = f.read()
    except Exception as e:
        logger.error(f"Failed to read {file_path}: {e}")
        return [], []

    # Split content by "time:" markers, but keep the delimiter or reconstruct
    # Easier strategy: Split by the time pattern, capturing the time value
    # The result will be [preamble, time1, block1, time2, block2, ...]
    parts = TIME_PATTERN.split(content)
    
    # parts[0] is usually empty or preamble
    # parts[1] is time1, parts[2] is block1, etc.
    
    transitions: List[TransitionData] = []
    lifetimes_list: List[LifetimeData] = []

    # Iterate in pairs of (time, block)
    for i in range(1, len(parts), 2):
        try:
            time_val = float(parts[i])
            block = parts[i+1]
            
            matrix, lifetime_dict = _parse_time_block(block)
            
            # Only add if we have valid data
            if matrix.size > 0:
                transitions.append({"time": time_val, "matrix": matrix})
            
            if lifetime_dict:
                lifetimes_list.append({"time": time_val, "lifetimes": lifetime_dict})
                
        except (ValueError, IndexError) as e:
            logger.debug(f"Skipping malformed block at index {i}: {e}")
            continue

    return transitions, lifetimes_list


def _parse_time_block(block: str) -> Tuple[np.ndarray, Dict[int, List[float]]]:
    """
    Parses a single time block for matrix and lifetimes.
    """
    lines = block.strip().splitlines()
    
    matrix_lines = []
    lifetime_dict: Dict[int, List[float]] = {}
    
    mode = "unknown" # Modes: "matrix", "lifetime"
    current_cluster_size = None
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
            
        # Check for headers
        if TRANSITION_HEADER_PATTERN.search(line):
            mode = "matrix"
            continue
        if LIFETIME_HEADER_PATTERN.search(line):
            mode = "lifetime"
            continue
            
        if mode == "matrix":
            # Parse matrix row: "0 1 5 2"
            try:
                # Filter out non-numeric lines (sometimes headers or noise appear)
                # The old code checked `not lines[i].startswith(('A', 'B', 'C'))`
                if line[0].isdigit() or line[0] == '-':
                    row = [int(x) for x in line.split()]
                    matrix_lines.append(row)
            except ValueError:
                continue
                
        elif mode == "lifetime":
            # Check for size header: "size of the cluster: 1"
            size_match = SIZE_HEADER_PATTERN.search(line)
            if size_match:
                current_cluster_size = int(size_match.group(1))
            elif current_cluster_size is not None:
                # Parse lifetimes: "0.1 0.2 0.5"
                try:
                    vals = [float(x) for x in line.split()]
                    if current_cluster_size not in lifetime_dict:
                        lifetime_dict[current_cluster_size] = []
                    lifetime_dict[current_cluster_size].extend(vals)
                except ValueError:
                    continue

    # Convert matrix to numpy array
    if matrix_lines:
        # Ensure rectangular shape
        max_len = max(len(row) for row in matrix_lines)
        padded_matrix = np.zeros((len(matrix_lines), max_len), dtype=int)
        for idx, row in enumerate(matrix_lines):
            padded_matrix[idx, :len(row)] = row
        matrix = padded_matrix
    else:
        matrix = np.array([])

    return matrix, lifetime_dict


def parse_copy_numbers(file_path: Path) -> pd.DataFrame:
    """
    Parses copy_numbers_time.dat (CSV format).
    """
    if not file_path.exists():
        logger.warning(f"File not found: {file_path}")
        return pd.DataFrame()
    
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        logger.error(f"Error reading {file_path}: {e}")
        return pd.DataFrame()


def parse_complex_histogram(file_path: Path) -> tuple[npt.NDArray[np.float64], list[dict], csc_array]:
    """
    Parses histogram_complexes_time.dat into a column-efficient sparse matrix.
    
    Returns:
        Tuple (times, compositions, histogram_matrix) of types (np.ndarray, list[dict], scipy.csc_array)
    """
    if not file_path.exists():
        return []

    with open(file_path, 'r') as f:
        content = f.read()
        
    parts = TIME_PATTERN.split(content)
    data = []
    all_comps = [] # for keeping track of every distinct composition seen so far
    
    for i in range(1, len(parts), 2):
        try:
            time_val = float(parts[i])
            block = parts[i+1]
            
            complexes = []
            for line in block.strip().splitlines():
                if not line.strip():
                    continue
                
                # Line format: "4\tB: 5."
                # Or "10\tC: 1. A: 1."
                try:
                    parts_line = line.split('\t', 1)
                    if len(parts_line) < 2:
                        # Sometimes split by space?
                        parts_line = line.split(None, 1)
                        if len(parts_line) < 2:
                             continue
                    
                    count = int(parts_line[0])
                    comp_str = parts_line[1]
                    
                    # Parse composition: "C: 1. A: 1."
                    # Regex finds "Key: Value" pairs
                    comp_dict = {}
                    for match in re.finditer(r"([A-Za-z0-9]+):\s*(\d+)", comp_str):
                        comp_dict[match.group(1)] = int(match.group(2))
                        
                    complexes.append({'count': count, 'composition': comp_dict})
                    if comp_dict not in all_comps:
                        all_comps.append(comp_dict)
                except ValueError:
                    continue
            
            if complexes:
                data.append({'time': time_val, 'complexes': complexes})
                
        except (ValueError, IndexError):
            continue
    
    time_dim = len(data)
    comp_dim = len(all_comps)

    time_values = np.zeros(time_dim)

    # build matrix using LIL, then convert to CSC for better column slicing efficiency
    hist_matrix = lil_array((time_dim,comp_dim))
    for i,step in enumerate(data):
        time_values[i] = step['time']
        row = np.zeros(comp_dim)
        for comp in step['complexes']:
            row[all_comps.index(comp['composition'])] = comp['count']
        hist_matrix[i,:] = row

    hist_matrix = hist_matrix.tocsc()
    
    return time_values, all_comps, hist_matrix



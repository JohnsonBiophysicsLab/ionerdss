import pytest
from pathlib import Path
import pandas as pd
import numpy as np
from io import StringIO

@pytest.fixture
def sample_transition_file(tmp_path):
    """Creates a dummy transition_matrix_time.dat file."""
    content = """time: 0.0
transion matrix for each mol type: 
0 0
0 0
lifetime for each mol type: 
time: 0.1
transion matrix for each mol type: 
0 2
1 0
lifetime for each mol type: 
size of the cluster: 2
0.5 0.5
"""
    p = tmp_path / "transition_matrix_time.dat"
    p.write_text(content)
    return p

@pytest.fixture
def sample_copy_numbers_file(tmp_path):
    """Creates a dummy copy_numbers_time.dat file."""
    content = """Time (s),A,B,Complex
0.0,10,10,0
0.1,8,8,2
0.2,5,5,5
"""
    p = tmp_path / "copy_numbers_time.dat"
    p.write_text(content)
    return p

@pytest.fixture
def sample_complex_histogram_file(tmp_path):
    """Creates a dummy histogram_complexes_time.dat file."""
    content = """Time (s): 0.0
10\tA: 1.
Time (s): 0.1
5\tA: 1. B: 1.
2\tC: 3.
"""
    p = tmp_path / "histogram_complexes_time.dat"
    p.write_text(content)
    return p

@pytest.fixture
def mock_simulation_dir(tmp_path, sample_transition_file, sample_copy_numbers_file, sample_complex_histogram_file):
    """Creates a mock simulation directory structure."""
    # Structure: root/1/DATA/files...
    sim_dir = tmp_path / "1"
    data_dir = sim_dir / "DATA"
    data_dir.mkdir(parents=True)
    
    # Copy/Link files
    (data_dir / "transition_matrix_time.dat").write_text(sample_transition_file.read_text())
    (data_dir / "copy_numbers_time.dat").write_text(sample_copy_numbers_file.read_text())
    (data_dir / "histogram_complexes_time.dat").write_text(sample_complex_histogram_file.read_text())
    
    return tmp_path



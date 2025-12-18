import pytest
import numpy as np
import pandas as pd
from ionerdss.analysis.io import parser

def test_parse_transition_file(sample_transition_file):
    transitions, lifetimes = parser.parse_transition_file(sample_transition_file)
    
    # Check transitions
    assert len(transitions) == 2
    assert transitions[0]['time'] == 0.0
    assert np.all(transitions[0]['matrix'] == np.zeros((2,2)))
    
    assert transitions[1]['time'] == 0.1
    expected_mat = np.array([[0, 2], [1, 0]])
    assert np.all(transitions[1]['matrix'] == expected_mat)
    
    # Check lifetimes
    assert len(lifetimes) == 1 # Only one time point had lifetimes
    assert lifetimes[0]['time'] == 0.1
    assert lifetimes[0]['lifetimes'][2] == [0.5, 0.5]

def test_parse_copy_numbers(sample_copy_numbers_file):
    df = parser.parse_copy_numbers(sample_copy_numbers_file)
    assert not df.empty
    assert len(df) == 3
    assert 'Complex' in df.columns
    assert df.iloc[1]['A'] == 8

def test_parse_complex_histogram(sample_complex_histogram_file):
    data = parser.parse_complex_histogram(sample_complex_histogram_file)
    assert len(data) == 2
    
    t0 = data[0]
    assert t0['time'] == 0.0
    assert len(t0['complexes']) == 1
    assert t0['complexes'][0]['count'] == 10
    assert t0['complexes'][0]['composition'] == {'A': 1}
    
    t1 = data[1]
    assert t1['time'] == 0.1
    assert len(t1['complexes']) == 2
    # Check "5 A: 1. B: 1."
    c1 = t1['complexes'][0]
    assert c1['count'] == 5
    assert c1['composition'] == {'A': 1, 'B': 1}


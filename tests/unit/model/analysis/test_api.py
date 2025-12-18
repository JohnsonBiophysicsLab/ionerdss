import pytest
from ionerdss.analysis import Analyzer
import pandas as pd

def test_analyzer_loading(mock_simulation_dir):
    analyzer = Analyzer(mock_simulation_dir)
    
    assert len(analyzer.simulations) == 1
    sim = analyzer.get_simulation(0)
    
    assert sim.id == "1"
    
    # Test lazy loading
    assert sim._data is None
    sim.load()
    assert sim._data is not None
    assert len(sim.data.transitions) == 2
    assert sim.data.copy_numbers is not None

def test_analyzer_integration_compute(mock_simulation_dir):
    analyzer = Analyzer(mock_simulation_dir)
    sim = analyzer.get_simulation(0)
    
    # Compute Free Energy (should trigger load)
    df_fe = analyzer.compute_free_energy(sim)
    
    assert not df_fe.empty
    assert 'free_energy' in df_fe.columns
    
    # Check caching
    assert sim.data.df_free_energy is not None
    assert sim.data.df_free_energy is df_fe


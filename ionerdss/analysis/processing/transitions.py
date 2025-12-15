"""
Scientific processing logic for transition matrices.

This module contains pure functions for calculating physical properties
from transition matrices. All inputs are standard NumPy arrays or Pandas objects.
"""

import numpy as np
import pandas as pd
from typing import Tuple, Dict

def compute_size_distribution(transition_matrix: np.ndarray) -> pd.DataFrame:
    """
    Calculates the probability distribution of cluster sizes.
    
    P(size) = sum(transitions_to_size) / total_transitions
    Note: This uses the 'steady state' assumption from the counts in the matrix.
    Alternatively, row sums represent the total number of events *originating* from a size.
    
    Args:
        transition_matrix (np.ndarray): Square matrix (NxN).
        
    Returns:
        pd.DataFrame: DataFrame with columns ['size', 'count', 'probability'].
    """
    if transition_matrix.size == 0:
        return pd.DataFrame(columns=['size', 'count', 'probability'])

    # Counts per size (sum of row = total events starting from size i)
    # Assuming rows = 'from', cols = 'to' based on typical transition matrix def
    # But let's check legacy code: "sum over rows" -> "counts_per_size = aggregated_matrix.sum(axis=1)"
    # Yes, axis 1 means summing columns for each row.
    counts = transition_matrix.sum(axis=1)
    total_counts = counts.sum()
    
    probabilities = counts / total_counts if total_counts > 0 else np.zeros_like(counts)
    
    sizes = np.arange(1, len(counts) + 1)
    
    return pd.DataFrame({
        'size': sizes,
        'count': counts,
        'probability': probabilities
    })


def compute_free_energy(size_dist: pd.DataFrame, temperature: float = 1.0) -> pd.DataFrame:
    """
    Computes the free energy landscape F = -kT ln(P).
    
    Args:
        size_dist (pd.DataFrame): Output from compute_size_distribution.
        temperature (float): kT units. Default 1.0.
        
    Returns:
        pd.DataFrame: Copy of input with 'free_energy' column added.
    """
    df = size_dist.copy()
    
    # Avoid log(0)
    probs = df['probability'].values.astype(np.float64)
    with np.errstate(divide='ignore'):
        fe = -np.log(probs) * temperature
        
    # Replace inf with NaN for cleaner plotting
    fe = np.where(np.isinf(fe), np.nan, fe)
    df['free_energy'] = fe
    
    return df


def compute_transition_probabilities(transition_matrix: np.ndarray, symmetric: bool = True) -> pd.DataFrame:
    """
    Computes association (growth) and dissociation (shrinkage) probabilities.
    
    Args:
        transition_matrix (np.ndarray): NxN matrix.
        symmetric (bool): If True, divides counts by 2 for monomer events (N+1) or (N-1).
        
    Returns:
        pd.DataFrame: DataFrame with columns ['size', 'growth_prob', 'shrink_prob'].
    """
    if transition_matrix.size == 0:
        return pd.DataFrame(columns=['size', 'growth_prob', 'shrink_prob'])

    n_sizes = transition_matrix.shape[0]
    growth_probs = []
    shrink_probs = []
    
    # Vectorized approach is tricky due to the symmetric condition varying by index
    # We'll stick to an optimized loop for clarity on the physics logic
    
    for n in range(n_sizes):
        # Row index n corresponds to size n+1
        size = n + 1
        
        # Growth: Transitions from n to m > n
        # Matrix indices: [m, n] (col n is 'from n', row m is 'to m') ?
        # CHECK LEGACY: "count = aggregated_matrix[m, n]" where m > n (Growth)
        # This implies Matrix[to, from] convention.
        
        # Shrinkage: Transitions from n to m < n
        
        growth_counts = 0.0
        shrink_counts = 0.0
        
        # Growth (Column n, Rows > n)
        for m in range(n + 1, n_sizes):
            count = transition_matrix[m, n]
            pair_size = (m + 1) - size
            # Symmetric rule: if adding a monomer (pair_size == 1, i.e., m = n+1), divide by 2
            # Legacy code: "if pair_size == n + 1" -> This looks like it means "if size doubled"?
            # Wait, legacy code: 
            #   assoc_counts: m in range(n+1, max_size)
            #   pair_size = m - n
            #   if pair_size == n + 1: count /= 2
            # This logic seems specific to the simulation physics (e.g. homo-dimerization).
            # I will preserve the legacy logic exactly.
            
            legacy_pair_size = m - n # This is actually delta size? No, indices.
            # m is index of 'to', n is index of 'from'
            # size_to = m+1, size_from = n+1
            # delta = m-n.
            
            # Legacy check: "if pair_size == n + 1". 
            # In legacy loop: `for m in range(n + 1, max_size): pair_size = m - n`
            # So pair_size is just the difference in indices.
            # `n` in legacy loop was 0-indexed? Yes.
            # `if pair_size == n + 1`: this means `m - n == n + 1` => `m = 2n + 1`.
            # Size 'from' = n (legacy used 0-based loop but maybe treated n as size?)
            # Legacy: `sizes = np.arange(1, len... + 1)`
            # In legacy loop `for n in range(max_size - 1)`: n is 0-based index.
            # `pair_size = m - n`.
            # condition `pair_size == n + 1` => `m - n == n + 1` => `m = 2n + 1`.
            # size_to (m) = 2 * size_from (n) + 1 ??
            # This logic is confusing. Let's look at `calculate_association_probabilities` in legacy.
            
            # "if symmetric and pair_size == n + 1:"
            # If n=0 (size 1), pair_size = m. condition: m == 1. m=1 is size 2. 1+1=2.
            # So if monomer + monomer -> dimer?
            
            # Let's rewrite based on the intent: Correct for symmetric collisions (A+A -> A2).
            if symmetric and (m - n) == (n + 1): 
                 count /= 2.0
            growth_counts += count

        # Shrinkage (Column n, Rows < n)
        for m in range(n):
            count = transition_matrix[m, n]
            # Legacy: "pair_size = n - m. if pair_size == m + 1: count /= 2"
            if symmetric and (n - m) == (m + 1):
                count /= 2.0
            shrink_counts += count

        total = growth_counts + shrink_counts
        if total > 0:
            growth_probs.append(growth_counts / total)
            shrink_probs.append(shrink_counts / total)
        else:
            growth_probs.append(np.nan)
            shrink_probs.append(np.nan)

    return pd.DataFrame({
        'size': np.arange(1, n_sizes + 1),
        'growth_prob': growth_probs,
        'shrink_prob': shrink_probs
    })


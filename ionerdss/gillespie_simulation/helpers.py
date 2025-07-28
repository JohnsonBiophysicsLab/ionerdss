"""
Helper functions for Gillespie simulations

@author: MSANG
@email: msang2@jh.edu
"""

import numpy as np

def align_gillespie_simulations(t_points_all, y_points_all, target_time_points):
    """
    Aligns multiple Gillespie simulation datasets to a common set of time points.

    This function takes a list of simulation datasets, where each dataset is a tuple
    containing a time array and a corresponding data array. It uses linear interpolation
    to resample each dataset at the specified target time points.

    Args:
        t_points_all (list of array): A list of all time point arrays
        y_points_all (list of array): A list of all y arrays. A y array is in format of 
                                      [[N1, N2, N3, ...], ...]
        target_time_points (np.ndarray): An array of the desired time points to
                                         which all simulations should be aligned.

    Returns:
        np.ndarray: A 2D numpy array where each row corresponds to a simulation
                    and each column corresponds to a target time point.
    """
    aligned_data = []
    for time_points, data_values in zip(t_points_all, y_points_all):
        # Use numpy's interpolation function.
        # It requires that the new time points are within the range of the original time points.
        # We handle cases where target_time_points might extend beyond the simulation time
        # by specifying left and right fill values (using the first and last data values).
        if len(np.shape(data_values)) == 1:
            aligned_values = np.interp(
                target_time_points,
                time_points,
                data_values,
                left=data_values[0],
                right=data_values[-1]
            )
            aligned_data.append(aligned_values)
        elif len(np.shape(data_values)) == 2:
            aligned_values = []
            for data_values_species in np.transpose(data_values):
                aligned_values_species = np.interp(
                    target_time_points,
                    time_points,
                    data_values_species,
                    left=data_values_species[0],
                    right=data_values_species[-1]
                )
                aligned_values.append(aligned_values_species)
            aligned_data.append(np.transpose(aligned_values))
    return np.array(aligned_data)
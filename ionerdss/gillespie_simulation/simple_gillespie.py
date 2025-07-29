#!/usr/bin/env python3.8.10
"""
This script contains functions for simulating simple
gillespie trajectories.

Usage: from <dir>.reaction_gillespie import *

Docstrings and comments are helped written by GPT3.5 and GPT API

@author: MYING
@email: yying7@jh.edu

@Modified by MSANG
@email: msang2@jh.edu
"""

# >>>>>> For autocompletion >>>>>>
__all__ = [
    "run_Gillespie",
    "convert_to_microscopic_rate_constants",
    "rate_constants_volume_correction",
    "calculate_propensity",
    "update_rates",
    "gillespie_simulation",
    "align_gillespie_simulations",
    "run_Gillespie_repeats",
]
def __dir__():
    return __all__
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



import numpy as np
import math
from .helpers import align_gillespie_simulations

def run_Gillespie_repeats(
        N_repeats: int = 1,
        show_progress: bool = False,
        target_time_points: np.ndarray = None,
        *args, **kwargs,
    ):
    """
    Args:
        N_repeats (int): The number of simulation repeats to run.
        show_progress (bool): True for showing progressbar
        target_time_points (numpy.ndarray): Time points to which the
            simulation data are aligned. If None, the time points from the
            first trajectory are used as the target. Defaults to None.
        *args: Positional arguments to be passed directly to `run_Gillespie`.
        **kwargs: Keyword arguments to be passed directly to `run_Gillespie`.
                  See the `run_Gillespie` function's docstring for a detailed
                  list of available arguments. The 'seed' argument is handled
                  specially by this wrapper to ensure different trajectories
                  are generated if a seed is provided.

    Returns:
        tuple: A tuple containing two numpy arrays:
               - aligned_t: The common time points.
               - aligned_y: A list of the state arrays, each aligned to `aligned_t`.

    Example:
            
    Note:
        
    """
    # Make a copy of kwargs to avoid modifying the original dictionary passed by the user.
    gillespie_kwargs = kwargs.copy()

    # Handle the random seed. We set the seed once here, then remove it from
    # kwargs so that each call to run_Gillespie within the loop produces a
    # different random trajectory.
    if 'seed' in gillespie_kwargs:
        if gillespie_kwargs['seed'] is not None:
            np.random.seed(gillespie_kwargs.pop('seed'))

    y_all = []
    t_all = []
    
    if show_progress:
        from progressbar import progressbar
        iterator = progressbar(range(N_repeats)) 
    else: 
        range(N_repeats)

    for _ in iterator:
        # Pass the arguments down to the core simulation function.
        # Note that the 'seed' kwarg has already been removed.
        y_record, t_record = run_Gillespie(
            *args, **gillespie_kwargs
        )
        t_all.append(t_record)
        y_all.append(y_record)

    # If no target time points are provided, use the first simulation's times.
    if target_time_points is None:
        if not t_all:
             # Handle case with zero repeats
            return np.array([]), []
        target_time_points = t_all[0]

    aligned_y = align_gillespie_simulations(t_all, y_all, target_time_points)

    return target_time_points, aligned_y

def run_Gillespie(
        max_time:float, 
        y_init: np.ndarray, 
        reactant_matrix: np.ndarray, 
        product_matrix: np.ndarray, 
        rate_constants: np.ndarray,
        volume:float = None,
        macroscopic:bool = None,
        volume_corrected:bool = None,
        record_interval:float = 1,
        full_update_scheme:bool = False,
        rate_update_rules:callable = None,
        avogadro=6.02214e23,
        seed:int=None,
        excluded_volume:float=None,
        excluded_species_id:int=None,
    ):
    """
    Args:
        max_time (float): Maximum simulation time.
        y_init (numpy.ndarray): Initial state of the system (species counts).
        reactant_matrix (numpy.ndarray): Matrix representing reactants in each reaction.
        product_matrix (numpy.ndarray): Matrix representing products in each reaction.
        rate_constants (np.ndarray): Rate constants for each reaction.
        volume (float): The volume of the system. Make the units consistent with rates.
        macroscopic (bool): True for input rates are macroscopic rates.
        volume_corrected (bool): True for input rates are volume corrected rates.
        record_interval (int): Step interval for saving data. Default `record_interval=1`.
        full_update_scheme (bool): controls if update every propensity entry in each iteration. Defualt `full_update_scheme=False`.
        rate_update_rules (callable): How to update kinetic rates. 
            It takes inputs (rate_constants, y_curr, reactant_matrix, volume)
        avogadro=6.02214e23, set `avogadro=1` if Molar is not used in units.
        seed (int): Set random seed if given. 
        excluded_volume (float): Consider excluded volume in 1D or other reasonable systems.
        excluded_species_id (int): Index of the species that creates excluded volume.

    Returns:
        tuple: A tuple containing arrays for recorded time points (t_record) and
            corresponding system states (y_record).

    Example:
        Suppose you have the following input matrices and arrays:

        max_time = 100.0
        y_init = np.array([10, 5, 3])  # Initial state (species counts)
        
        reactant_matrix = np.array([[2, 1, 0],  # Example reactant matrix
                                    [0, 1, 1]])

        product_matrix = np.array([[0, 1, 0],  # Example product matrix
                                [1, 0, 1]])

        volume_corrected_rate_constants = np.array([0.1, 0.05])  # Example rate constants

        y_record, t_record = run_Gillespie(
            max_time, y_init, reactant_matrix, product_matrix, 
            volume_corrected_rate_constants, volume_corrected=True
        )
        print(y_record)
        print(t_record)
            
    Note:
        This function performs a Gillespie simulation for a chemical reaction system.
        It records the system state and corresponding time points during the simulation.
    """

    # Validate rate_update_rules
    need_update_rates = False
    if (rate_update_rules is not None) and (not macroscopic):
        raise ValueError('To update rates, macroscopic rates must be provided.')
    elif (rate_update_rules is not None) and (bool(macroscopic)):
        need_update_rates = True
    # print('[DEBUG] rate_update_rules:', rate_update_rules)
    # print('[DEBUG] need_update_rates:', need_update_rates)
        
    # process excluded volume
    if excluded_volume is not None:
        if excluded_species_id is None:
            raise ValueError('Must give a species that excludes others.')
        def reduced_volume(y_curr, excluded_volume, volume): 
            new_vol = volume - y_curr[excluded_species_id]*excluded_volume
            if new_vol == 0:
                # when there is no place to move, consider the available 
                # space is infinitesmal to the excluded volume
                return excluded_volume / 1000
            else:
                return new_vol
    else:
        reduced_volume = lambda y_curr, excluded_volume, volume: volume

    # calculate volume_corrected_rates based on inputs
    if (macroscopic is None) and (volume_corrected is None):
        print('By default, rate constants are volume corrected (with unit 1/s)')
        volume_corrected_rates = rate_constants 
    elif bool(macroscopic) ^ bool(volume_corrected):
        if macroscopic:
            new_volume = reduced_volume(y_init, excluded_volume, volume)
            if need_update_rates:
                volume_corrected_rates = update_rates(
                    rate_constants, y_init, new_volume, rate_update_rules, reactant_matrix, avogadro
                )
            else:
                volume_corrected_rates = rate_constants_volume_correction(
                    rate_constants, reactant_matrix, new_volume, avogadro
                )
        else:
            volume_corrected_rates = rate_constants
    else:
        raise ValueError("The rates must be either macroscopic or volume corrected, but not both. By default, they are volume corrected (with unit 1/s)")
    
    # Validate "reactant matrix", it should contain only integers
    if not np.all(np.mod(reactant_matrix, 1) == 0):
        raise ValueError("For gillespie, all entries in the reactant matrix must be mathematically integers.")
    # Validate "product matrix", it should contain only integers
    if not np.all(np.mod(product_matrix, 1) == 0):
        raise ValueError("For gillespie, all entries in the product matrix must be mathematically integers.")
    
    # prepare simulation
    time = 0.0 # Simulation time elapsed
    y = np.array([x for x in y_init])  # Initial copy numbers
    is_propensity_update_needed = np.zeros(len(reactant_matrix)) # 1 for updated needed
    delta_y = product_matrix - reactant_matrix  # Yield matrix
    index = np.array(range(0, len(reactant_matrix)))  # np.random.choice must be 1-d array; use indexing instead
    y_record = [np.copy(y)]  # Record array for copy numbers
    t_record = [time]  # Record array for time
    n_steps = 0 # Record every record_interval step(s)
    if seed is not None: np.random.seed(seed) # set random seed if given

    while time < max_time:  # Control simulation time scale

        # update volume if needed
        if excluded_volume is not None:
            new_volume = reduced_volume(y, excluded_volume, volume)
            volume_corrected_rates = rate_constants_volume_correction(
                rate_constants, reactant_matrix, new_volume, avogadro
            )
        else:
            new_volume = volume
        
        # update rates if needed (adaptive rates)
        if need_update_rates: volume_corrected_rates = update_rates(
                rate_constants, y, new_volume, rate_update_rules, reactant_matrix, avogadro
            )
        # print('[DEBUG] volume corrected rates:', volume_corrected_rates)
        if full_update_scheme:
            # Calculate propensity
            propensities = calculate_propensity(y, reactant_matrix, volume_corrected_rates)
        
        else:
            # Calculate propensity
            propensities = calculate_propensity(y, reactant_matrix, volume_corrected_rates,
                                            propensities, is_propensity_update_needed)
        # print(volume_corrected_rates)
        # print(propensities)
        # Calculate r_tot and sojourn time
        r_tot = np.sum(propensities)
        tau = - (1.0 / r_tot) * np.log(np.random.rand())

        # Choose reaction and add to species
        reaction_index_chose = np.random.choice(index, p=propensities / r_tot)
        y += delta_y[reaction_index_chose]
        
        
        # Update which propensities need to be updated in next iteration
        if not full_update_scheme:
            ##@deprecated: entries_changed = np.abs(np.transpose(delta_y[[reaction_index_chose]]))
            ##@deprecated: is_propensity_update_needed = np.squeeze(np.matmul(reactant_matrix, entries_changed))
            entries_changed = np.abs(delta_y[reaction_index_chose].T)
            is_propensity_update_needed = np.dot(reactant_matrix, entries_changed.squeeze())

        # Progress time
        time += tau

        # Record
        if n_steps % record_interval == 0:
            y_record.append(np.copy(y))
            t_record.append(time)
        n_steps += 1

    return np.array(y_record), np.array(t_record)
    
def convert_to_microscopic_rate_constants(macroscopic_rate_constants, reactant_matrix, volume, avogadro=6.02214e23):
    """
    Legacy API to rate_constants_volume_correction
    """
    return rate_constants_volume_correction( macroscopic_rate_constants,  reactant_matrix, volume, avogadro)

def rate_constants_volume_correction(
        macroscopic_rate_constants, 
        reactant_matrix, 
        volume,
        avogadro=6.02214e23
        ):
    """
    Note:
        All volumes are assumed to be in liters, and concentrations are assumed to be in mol/L!
    
    Convert macroscopic rate constants to microscopic rate constants for Gillespie algorithm.

    Args:
        macroscopic_rate_constants (numpy.ndarray): Array of macroscopic rate constants.
        reactant_matrix (numpy.ndarray): Matrix representing reactants in each reaction.
        volume (float): Volume of the system (assumed to be in liters).
        avogadro (float, optional): Avogadro's number (default: 6.02214e-23).

    Raises:
        ValueError: If any entry in the reactant matrix is not a mathematical integer.

    Returns:
        numpy.ndarray: Array of microscopic rate constants.
    """
    #### ALL volumes are assumed to be in liters, and concentrations are assumed to be in mol/L!
    # Check if all entries in the reactant matrix are mathematical integers
    if not np.all(np.mod(reactant_matrix, 1) == 0):
        raise ValueError("For Gillespie, all entries in the matrix must be mathematical integers.")

    # Initialize an array for microscopic rate constants
    microscopic_rate_constants = np.zeros(len(macroscopic_rate_constants))

    # Calculate microscopic rate constants
    for reaction_index, reaction in enumerate(reactant_matrix):
        scalar = 1
        power = 1
        for species_index, species_count in enumerate(reaction):
            scalar *= math.factorial(int(species_count))
            power -= species_count
        microscopic_rate_constants[reaction_index] = (
            scalar * macroscopic_rate_constants[reaction_index] * np.power((volume * avogadro), power)
        )

    return microscopic_rate_constants

def calculate_propensity(
        y, 
        reactant_matrix, 
        microscopic_rate_constants,
        previous_propensities = None, 
        is_propensity_update_needed = None
        ):
    """
    Calculate propensities for Gillespie algorithm.

    Args:
        y (numpy.ndarray): Current state of the system (species counts).
        reactant_matrix (numpy.ndarray): Matrix representing reactants in each reaction.
        microscopic_rate_constants (numpy.ndarray): Rate constants for each reaction.

    Returns:
        numpy.ndarray: Array of propensities for each reaction.

    Example:
        Suppose you have the following input matrices and arrays:

        y = np.array([10, 5, 3])  # Current state (species counts)
        
        reactant_matrix = np.array([[2, 1, 0],  # Example reactant matrix
                                    [0, 1, 1]])

        microscopic_rate_constants = np.array([0.1, 0.05])  # Example rate constants

        propensities = calculate_propensity(y, reactant_matrix, microscopic_rate_constants)
        print(propensities)
        # Output: [0.1 * comb(10, 2) * comb(5, 1), 0.05 * comb(5, 1) * comb(3, 1)]

    Note:
        The function calculates the propensity of each reaction in a Gillespie algorithm.
        Propensity is the product of the microscopic rate constant and combinatorial terms
        based on the reactant matrix and current state (y) of the system.
    """
    
    if previous_propensities is None: # simple case
        
        propensities = np.zeros(len(reactant_matrix))

        # Loop over each reaction
        for reaction_index, reaction in enumerate(reactant_matrix):
            propensity = microscopic_rate_constants[reaction_index]

            # Multiply by the combinatorial term for each reactant
            for species_index, species_count in enumerate(reaction):
                propensity *= math.comb(y[species_index], species_count)

            propensities[reaction_index] = propensity
        
        return propensities

    else: # optimization case
        # Loop over each reaction
        for reaction_index, reaction in enumerate(reactant_matrix):
            if is_propensity_update_needed[reaction_index] != 0:
                # Update propensity if entry is not zero
                propensity = microscopic_rate_constants[reaction_index]

                # Multiply by the combinatorial term for each reactant
                for species_index, species_count in enumerate(reaction):
                    propensity *= math.comb(y[species_index], species_count)
                
                previous_propensities[reaction_index] = propensity
            
        return previous_propensities
    
def update_rates(
        rate_constants, 
        y_curr, 
        volume, 
        rate_update_rules, 
        reactant_matrix,
        avogadro=6.02214e23
        ):
    """
    For systems that kinetic rates change as concentration changes (e.g. 1D and 2D), 
    update according to macroscopic rates and recalculate propensities.

    Args:
        rate_constants (numpy.ndarray): Array of rate constants.
            These are the intrinsic rates, i.e. the independent variable is at a limit that 
            the kinetic rates do not depend on them. For example, for 2D and 1D systems
            that diffusion affects rates, diffusion constants are infinite.
        y_curr (numpy.ndarray): Current state of the system (species counts).
        volume (float): Volume of the system (assumed to be in liters).
        rate_update_rules (callable): How to update kinetic rates. 
            It takes inputs (macroscopic_rate_constants, y_curr, reactant_matrix, volume)
        reactant_matrix
        avogadro (float, optional): Avogadro's number (default: 6.02214e-23).
    
    Returns:
        volume_corrected_rates
    """
    new_macro_rates = rate_update_rules(rate_constants, y_curr, reactant_matrix)
    return rate_constants_volume_correction(new_macro_rates, reactant_matrix, volume, avogadro)

def gillespie_simulation(max_time, y_init,reactant_matrix, product_matrix, microscopic_rate_constants,record_interval = 1,full_update_scheme = False):
    """
    Legacy API to gillespie_simulation_volume_corrected_rates
    """
    return run_Gillespie(
        max_time, y_init, reactant_matrix, product_matrix, microscopic_rate_constants, 
        volume_corrected=True, record_interval=record_interval, full_update_scheme=full_update_scheme,
    )


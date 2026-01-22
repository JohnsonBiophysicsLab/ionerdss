"""
Author MYING
Using Python 3.9.7
"""
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
# figure size in inches
# rcParams['figure.figsize'] = 11,8

from scipy.integrate import solve_ivp

def calculate_macroscopic_reaction_rates(y, reactant_matrix, k):
    """
    Calculate macroscopic reaction rates for reactions based on the equation sum(k * product(conc^stoichiometry)).

    Parameters:
        y (array-like): The current concentrations of species.
        reactant_matrix (array-like): The matrix representing the stoichiometry of reactants in each reaction.
        k (array-like): The rate constants for each reaction.

    Returns:
        array-like: An array containing the macroscopic reaction rates for each reaction.
    """
    num_reactions = len(reactant_matrix)
    num_species = len(reactant_matrix[0])

    # Initialize rates
    reaction_rates = np.zeros(num_reactions)

    # Calculate rates for each reaction
    for reaction_index, reactant_array in enumerate(reactant_matrix):

        # Calculate rate as k[i] * Product of y[j]^reactant_stoichiometry
        reaction_rates[reaction_index] = k[reaction_index]  # Initialize rates with k
        for species in np.arange(num_species):
            reaction_rates[reaction_index] *= y[species] ** reactant_array[species]

    return reaction_rates

def dydt(t, y, reactant_matrix, product_matrix, k):
    """
    Define the function required by scipy.integrate.solve_ivp for solving a system of ODEs.

    Parameters:
        t (float): The current time.
        y (array-like): The current concentrations of species.
        reactant_matrix (array-like): The matrix representing the stoichiometry of reactants in each reaction.
        product_matrix (array-like): The matrix representing the stoichiometry of products in each reaction.
        k (array-like): The rate constants for each reaction.

    Returns:
        array-like: The rate of change of concentrations for each species at the given time.
    """
    net_change_matrix = product_matrix - reactant_matrix

    # Get the number of species (note that it is the secondary axis)
    num_species = len(reactant_matrix[0])

    # Initialize the rate of change vector
    dydt = np.zeros(num_species)

    # Calculate reaction rates
    reaction_rates = calculate_macroscopic_reaction_rates(y, reactant_matrix, k)

    # Iterate over species
    for species in np.arange(num_species):

        # Iterate over reactions and sum up net_change * rate
        for reaction_index, net_change_array in enumerate(net_change_matrix):
            dydt[species] += (net_change_array[species] * reaction_rates[reaction_index])

    return dydt

def solve_reaction_ode(dydt, t_span, y_initial, reactant_matrix, product_matrix, k,
                    plotting=True, dense_output=True, method="BDF", atol = 1e-4, plotting_sample_points=1000, species_names = None):
    """
    Solve a system of ordinary differential equations (ODEs) for a chemical reaction and optionally plot the results.

    Parameters:
        dydt (function): the target function representing the system of reactions
        t_span (tuple): A tuple specifying the time span (initial and final times) for integration.
        y_initial (array-like): The initial concentrations of species.
        reactant_matrix (array-like): The matrix representing the stoichiometry of reactants in each reaction.
        product_matrix (array-like): The matrix representing the stoichiometry of products in each reaction.
        k (array-like): The rate constants for each reaction.
        plotting (bool, optional): Whether to plot the results. Defaults to True.
        plotting_sample_points (int, optional): Number of points for plotting. Defaults to 1000.
        species_names (list, optional): Names of species for legend. If None, species are labeled as y0, y1, etc.

    Returns:
        time
        concentration
    """
    sol = solve_ivp(dydt, t_span, y_initial, args=(reactant_matrix, product_matrix, k), dense_output=dense_output, atol=atol, method=method)
    t = np.linspace(min(t_span), max(t_span), plotting_sample_points)
    y = sol.sol(t)
    
    if plotting:
        plt.plot(t, y.T)
        plt.xlabel('time')
        plt.ylabel('concentration')
        if species_names is None:
            plt.legend([f"y{i}" for i in range(len(y_initial))])
        else:
            plt.legend(species_names)
        plt.show()
    
    return t, y.T, species_names

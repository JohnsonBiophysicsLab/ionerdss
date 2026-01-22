#!/usr/bin/env python3.8.10
"""
This script contains a class for parsing and
manipulating chemical reaction strings.

Usage: from <dir>.reaction_string_parser import *

Docstrings and comments are helped written by GPT3.5 and GPT API

@author: MYING
@email: yying7@jh.edu
"""

import re
import numpy as np

class ReactionStringParser:
    """
    A class for parsing and manipulating chemical reaction strings.

    Note:
        The string processing gets rid of any whitespaces in the given in the
        input. Please consider changing whitespaces into underscores '_' in you
        would like to reserve naming.

    Note:
        If you would like to customize regex for reaction symbols, note that 
        the level of intepretation is reversible > right > left!

    Note:
        By default, all rate constants are required (e.g. a reversible reaction
        requires)

    Attributes (__init__):
        rightward_reaction_symbol (str, regex): rightward reaction symbol.
        leftward_reaction_symbol (str, regex): leftward reaction symbol.
        reversible_reaction_symbol (str, regex): reversible reaction symbol.
        reaction_rate_separator (str, regex): reaction rate separator.
        species_separator (str, regex): species separator.
        reaction_rate_value_assigner (str, regex): reaction rate value assigner.
        stoich_species_regex (str, regex): regex pattern for stoichiometry.
        is_rate_constant_required (bool): when set to False, automatically generate
                naming for missing rate constant names. Default to True.
        auto_rate_constant_name (str): automatic generated rate constant name string when
                naming is missing for rate constants.

    Methods:
        parse_reaction_string(reaction_string): Parse a reaction string into components.
        parse_stoichiometry_string(reactants_or_products_string): Parse stoichiometry strings.
        extract_species_dictionaries_from_reaction_strings(reaction_strings):
            Extract dictionaries and rate constants.
        parse_reaction_strings(reaction_strings, dtype=int,
                               sort_reactions_by=None, sort_species_by=None, VERBOSE_MODE=False):
            Parse and sort reaction strings.
        sort_by_rate_constants(reactant_matrix, product_matrix, rate_constant_names, sort_order ):
            Sort matrices based on rate constants.
        sort_by_species_names(reactant_matrix, product_matrix, species_names, sort_order ):
            Sort matrices based on species names.
    """

    # Regular expression of all parts
    __DEFAULT_RIGHTWARD_REACTION_SYMBOL = r"-+>"
    __DEFAULT_LEFTWARD_REACTION_SYMBOL = r"<-+"
    __DEFAULT_REVERSIBLE_REACTION_SYMBOL = r"<-+>"
    __DEFAULT_REACTION_RATE_SEPARATOR = r"[,;]"
    __DEFAULT_SPECIES_SEPARATOR = r'\+'
    __DEFAULT_REACTION_RATE_VALUE_ASSIGNER = "="  # Placeholder
    __DEFAULT_STOICH_SPECIES_REGEX = r"([\d.]+|\d+\s*\/\s*\d+)?\s*(\w+)"
    __DEFAULT_IS_RATE_CONSTANT_REQUIRED = False
    __DEFAULT_AUTO_RATE_CONSTANT_NAME = "k"
    
    # Controls autoassigning of rate consant name
    __RATE_CONSTANT_NAME_INDEX = 0

    # Debug mode
    DEBUG_MODE = False

    def __init__(self, **kwargs):
        """
        Note the level of intepretation is reversible > right > left!
        """
        self.__rightward_reaction_symbol = \
            kwargs.get('rightward_reaction_symbol',
                       self.__DEFAULT_RIGHTWARD_REACTION_SYMBOL)
        self.__leftward_reaction_symbol = \
            kwargs.get('leftward_reaction_symbol',
                       self.__DEFAULT_LEFTWARD_REACTION_SYMBOL)
        self.__reversible_reaction_symbol = \
            kwargs.get('reversible_reaction_symbol',
                       self.__DEFAULT_REVERSIBLE_REACTION_SYMBOL)
        self.__reaction_rate_separator = \
            kwargs.get('reaction_rate_separator',
                       self.__DEFAULT_REACTION_RATE_SEPARATOR)
        self.__species_separator = \
            kwargs.get('species_separator', self.__DEFAULT_SPECIES_SEPARATOR)
        self.__reaction_rate_value_assigner = \
            kwargs.get('reaction_rate_value_assigner',
                       self.__DEFAULT_REACTION_RATE_VALUE_ASSIGNER)
        self.__stoich_species_regex = \
            kwargs.get('stoich_species_regex',
                       self.__DEFAULT_STOICH_SPECIES_REGEX)
        self.__is_rate_constant_required = \
            kwargs.get('is_rate_constant_required',
                       self.__DEFAULT_IS_RATE_CONSTANT_REQUIRED)
        self.__auto_rate_constant_name = \
            kwargs.get('auto_rate_constant_name',
                       self.__DEFAULT_AUTO_RATE_CONSTANT_NAME)
        self.DEBUG_MODE = \
            kwargs.get("DEBUG_MODE", False)


    def __find_direction_and_split_reaction(self, reaction_string):
        """
        Search for the symbol indicating reaction direction and split reaction
        string to two parts by the direction indicating symbol.
        """
        # Check if the REVERSIBLE regular expression is contained in the string
        if re.search(self.__reversible_reaction_symbol, reaction_string):
            # Split the reaction string by the regular expression
            reaction_parts = re.split(
                self.__reversible_reaction_symbol, reaction_string)
            reaction_direction = 0  # 0 signifies reversible reaction

        # Check if the RIGHTWARD regular expression is contained in the string
        elif re.search(self.__rightward_reaction_symbol, reaction_string):
            # Split the reaction string by the regular expression
            reaction_parts = re.split(
                self.__rightward_reaction_symbol, reaction_string)
            reaction_direction = 1  # 1 signifies rightward reaction

        # Check if the LEFTWARD regular expression is contained in the string
        elif re.search(self.__leftward_reaction_symbol, reaction_string):
            reaction_parts = re.split(
                self.__leftward_reaction_symbol, reaction_string)
            reaction_direction = -1  # -1 signifies leftward reaction

        else:
            # Handle the case where the regular expression is not found in the string
            raise ValueError(
                "The regular expression indicating direction is not in the string." + 
                  reaction_string
            )
        return reaction_direction, reaction_parts

    def parse_reaction_string(self, reaction_string):
        """Parse a reaction string into left-hand-side species, right-hand-side species,
        rate constant, and reaction direction.

        Note that the rate constant output type can be a float or a tuple of
        two floats, depending on the direction of the reaction. The reaction direction is
        defines as (+1) for rightward reaction, (-1) for leftward reaction, and (0) for
        reversible reaction.

        Reactants and products are not strictly defined by the sequence however, as the
        output is given by (lefthandside_species, righthandside_species, ...). This is to
        accomadate reversible reactions that takes up two reactions where the definition
        of reactants and products is ambiguous.

        Args:
            reaction_string (str): The input reaction string.

        Returns:
            tuple: A tuple containing for elements -
                - left_stoich_species (str),
                - right_stoich_species (str),
                - rate constant (str or tuple(str) in case of reversible reaction),
                - reaction_direction (-1, 0, or 1).
                    -1 denotes leftward reaction
                    0 denotes reversible reaction
                    1 denotes rightward reaction

        Raises:
            ValueError: If the input reaction string does not have the expected format.

        Example usage:
            input_string = "A + 2B -> C, k1"
            left_stoich_species, right_stoich_species, rate_constant, direction = \
                    parse_reaction_string(input_string)
            print(left_stoich_species)
            print(right_stoich_species)  # Output: C
        """
        reaction_direction, reaction_parts = \
            self.__find_direction_and_split_reaction(reaction_string)

        if len(reaction_parts) != 2:
            raise ValueError(
                "Invalid reaction string format: " + reaction_string)

        # Extract species on LHS and the second half of the reaction string
        # .strip() gets rid of whitespaces!
        left_stoich_species = reaction_parts[0].strip() # stoich-species string on LHS
        second_half_string = reaction_parts[1].strip() # RHS + rate constants

        # Split the second half of the reaction string by reaction_rate_separator
        # (default to ',') to separate products and rate constant
        second_half_parts = re.split(
            self.__reaction_rate_separator, second_half_string)
        
        # If there is no rate constant name and rate constant required set to False
        # generate rate constant name
        # Reversible reactions are treated as two separated reactions:
        # one forward and one backward
        # this means we will have two rate constants
        if reaction_direction == 0:  # reversible reaction
            if len(second_half_parts) != 3:  # reverisble reactions have two rate constants
                if self.__is_rate_constant_required:
                    raise ValueError(
                        "Invalid reversible reaction string format [psa0]: " + reaction_string)
                else: # if rate constant not required, come up with unique naming
                    right_stoich_species = second_half_parts[0]
                    rate_constant_forward = self.__auto_rate_constant_name + str(self.__RATE_CONSTANT_NAME_INDEX)
                    rate_constant_backward = self.__auto_rate_constant_name + str(self.__RATE_CONSTANT_NAME_INDEX + 1)
                    self.__RATE_CONSTANT_NAME_INDEX += 2
            else: # regular case where there are two rate constants
                right_stoich_species, rate_constant_forward, rate_constant_backward =\
                    [part.strip() for part in second_half_parts[:3]]
            # put the forward and backward constant into a tuple
            rate_constant = (rate_constant_forward, rate_constant_backward)
        # Leftward or rightward reactions
        else:
            if len(second_half_parts) < 2:  # unidirection reactions have one rate constant
                if self.__is_rate_constant_required:
                    raise ValueError(
                        "Invalid unidirection reaction string format [psa1]: " + reaction_string)
                else: # if rate constant not required, come up with unique naming
                    right_stoich_species = second_half_parts[0]
                    rate_constant = self.__auto_rate_constant_name + str(self.__RATE_CONSTANT_NAME_INDEX)
                    self.__RATE_CONSTANT_NAME_INDEX += 1
            else: # regular case where there is one rate constant
                # Extract products and rate constant
                right_stoich_species, rate_constant = \
                    [part.strip() for part in second_half_parts[:2]]

        return left_stoich_species, right_stoich_species, rate_constant, reaction_direction

    # import re

    def parse_stoichiometry_string(self, reactants_or_products_string):
        """Parse a stoichiometry string into a stoichiometry dictionary.

        Note:
            - species with the same naming will be merged, such as A + 2A -> B will
              be merged into {'A' : 3.0}
            - Although while output into string, the value may show up as integers,
              ALL the stochiometry values are converted to float.
            - All the fractions are calculated and converted to float, such as
              1/2 A -> B will be converted to {'A' : 0.5}

        Args:
            reactants_or_products_string (str): The input stoichiometry string.

        Returns:
            dict: A dictionary representing the stoichiometry, where keys are elements and values are coefficients.

        Raises:
            ValueError: If the input stoichiometry string has an invalid format.

        Example:
            input_string = "A + 1.42857  B   +1/2C+  A"
            stoichiometry = parse_stoichiometry_string(input_string)
            print(stoichiometry)  # Output: {'A': 2, 'B': 1.42857  B, 'C': 0.5}
        """
        # Split the input string into individual reactants by the '+' symbol
        reactants = re.split(self.__species_separator,
                             reactants_or_products_string)
        if self.DEBUG_MODE:
            print(reactants)

        # Initialize an empty dictionary to store the stoichiometry
        stoichiometry_dict = {}

        for reactant in reactants:
            # Remove leading and trailing white spaces
            reactant = reactant.strip()

            # Use regular expression to match coefficients and elements in the reactant
            # The default regex supports decimals and fractions in stoichiometry
            match = re.match(self.__stoich_species_regex, reactant)

            # If match is found for regex
            if match:
                coefficient, element = match.groups()
                if coefficient:  # if coefficient is present
                    if '/' in coefficient:  # If coefficient contains fraction
                        # Separate the fraction into nominator and denominator by '/'
                        fraction_parts = coefficient.split('/')
                        # Only supports simple fraction A/B, output error otherwise
                        if len(fraction_parts) != 2:
                            raise ValueError(
                                "Invalid stoichiometry: " + coefficient)
                        # Calculate float value from fraction
                        coefficient = float(
                            fraction_parts[0])/float(fraction_parts[1])
                    else:  # If coefficient does not contain fraction
                        coefficient = float(coefficient)
                else:  # if there is not coefficient, stochiometry is 1
                    coefficient = 1.0

                # Update the stoichiometry dictionary with the coefficient and element
                if element in stoichiometry_dict:  # case of repeated species in string
                    stoichiometry_dict[element] += coefficient
                else:  # case of new species
                    stoichiometry_dict[element] = coefficient
            else:
                if re.search(r'\S', reactant):
                    raise ValueError(
                        f"Invalid reactant or product format: {reactant}")

        return stoichiometry_dict

    def extract_species_dictionaries_from_reaction_strings(self, reaction_strings):
        """
        Extract species dictionaries, rate constants, and names from a list of reaction strings.
        For reversible reaction, there will be two rate constants provided. Accordingly,
        the forward and backward reaction is treated as two separate reactions in
        the system and will be added as two reactions to the reactant and product list
        for the reactions.

        Warning:
            A warning message will be output in commandline in case of duplicate rate
            constant naming.

        Args:
            reaction_strings (list): List of reaction strings.

        Returns:
            tuple: A tuple containing species_names_set (set), rate_constant_names (list),
            reactant_dictionaries (list), and product_dictionaries (list).

        Example:
            reaction_strings = ["A + B -> C, kon", "2X -> Y, kf", "C -> A + B, koff", "Y + A -> X + C, ki"]
            extract_species_dictionaries_from_reaction_strings(reaction_strings)
        """
        species_names_set = set()
        reactant_dictionaries = []
        product_dictionaries = []
        # for future use, I think rate constant not from string is better for now?
        rate_constant_values = []
        rate_constant_names = []

        for reaction_string in reaction_strings:
            # split reaction string into four entries
            left_stoich_species, right_stoich_species, rate_constant, reaction_direction\
                = self.parse_reaction_string(reaction_string)

            # convert to species-stoichiometry dict for reactant and product
            left_dict = self.parse_stoichiometry_string(left_stoich_species)
            right_dict = self.parse_stoichiometry_string(right_stoich_species)

            # add new names to set
            species_names_set.update(left_dict.keys())
            species_names_set.update(right_dict.keys())

            # deal with different direction of reaction
            # forward reaction (going from left to right)
            if reaction_direction >= 0:

                # append dictionaries to the respective lists
                reactant_dictionaries.append(left_dict)
                product_dictionaries.append(right_dict)

                # append rate constant to the list
                # in case of reversible reaction, rate_constant will be a tuple
                # and we are taking the first one as forward reaction rate constant
                rate_constant_names.append(rate_constant[0]
                                           if isinstance(rate_constant, (list, tuple)) and rate_constant else rate_constant)

            # backward reaction (going from right to left)
            if reaction_direction <= 0:

                # append dictionaries to the respective lists
                reactant_dictionaries.append(right_dict)
                product_dictionaries.append(left_dict)

                # append rate constant to the list
                # in case of reversible reaction, rate_constant will be a tuple
                # and we are taking the first one as forward reaction rate constant
                rate_constant_names.append(rate_constant[1]
                                           if isinstance(rate_constant, (list, tuple)) and rate_constant else rate_constant)

        # duplicate rate constant naming check
        duplicates = [
            name for name in rate_constant_names if rate_constant_names.count(name) > 1]
        if duplicates:
            print("WARNING: Repeated names found in rate_constant_names:",
                  set(duplicates), "\n This can lead to ambiguity in ODE definition.",
                  "\n Please consider change the naming.")

        return species_names_set, rate_constant_names, reactant_dictionaries, product_dictionaries

    def parse_reaction_strings(self, reaction_strings, dtype=int,
                               sort_reactions_by=None,
                               sort_species_by="alphabetical",
                               VERBOSE_MODE=False):
        """
        Parse reaction strings into species names, rate constant names, reactant matrices
        and product matrices.

        Avoid using dtype = float as long as there is no decimals or fractions
        in stoichiometry is advised.

        Args:
            reaction_strings (list): List of reaction strings.
            dtype (type, optional): Data type for matrix values (default: int).
            sort_reactions_by (list, optional): List of rate constant names to sort reactions.
            sort_species_by (list, optional): List of species names to sort species.
            VERBOSE_MODE (bool, optional): If True, print additional information for debugging.

        Returns:
            tuple: A tuple containing species_names (list), rate_constant_names (list),
                reactant_matrix_array (numpy.ndarray), and product_matrix_array (numpy.ndarray).

        Examples:
            # Example usage:
            reaction_strings = ["A + B -> C, kon", "2X -> Y, kf",
                                "C -> A + B, koff",
                                "Y + A -> X + C, ki"]
            species_names, reactant_matrix, product_matrix, rate_constant_names =\
                    parse_reaction_strings(reaction_strings)
        """

        rate_constant_values = []  # Placeholder

        # Extract species dictionaries, names, and count
        species_names_set, rate_constant_names, reactant_dictionaries, product_dictionaries = \
            self.extract_species_dictionaries_from_reaction_strings(
                reaction_strings)

        species_names = list(species_names_set)
        num_species = len(species_names)
        num_reactions = len(reaction_strings)

        if VERBOSE_MODE:
            print(f"species_names : {species_names}")
            print(f"reactant_dictionaries : {reactant_dictionaries}")
            print(f"product_dictionaries : {product_dictionaries}")

        # Initialize reactant and product matrices
        reactant_matrix = np.zeros(
            (num_reactions, num_species), dtype=dtype)
        product_matrix = np.zeros(
            (num_reactions, num_species), dtype=dtype)

        # Create a dictionary to map species names to indices
        species_to_index = {species: i for i,
                            species in enumerate(species_names)}

        # Fill in reactant and product matrices
        for i in range(num_reactions):
            for species, stoichiometry in reactant_dictionaries[i].items():
                reactant_matrix[i,
                                species_to_index[species]] = stoichiometry
            for species, stoichiometry in product_dictionaries[i].items():
                product_matrix[i,
                               species_to_index[species]] = stoichiometry

        # Sort reactoins by given list of rate constants
        if sort_reactions_by is not None:
            reactant_matrix, product_matrix, rate_constant_names =\
                self.sort_by_rate_constants(reactant_matrix,
                                            product_matrix,
                                            rate_constant_names,
                                            sort_reactions_by)

        # Sort species by given list of species
        if sort_species_by is not None:
            reactant_matrix, product_matrix, species_names =\
                self.sort_by_species_names(reactant_matrix,
                                           product_matrix,
                                           species_names,
                                           sort_species_by)

        return species_names, rate_constant_names, reactant_matrix, product_matrix

    def sort_by_rate_constants(self, reactant_matrix, product_matrix,
                               rate_constant_names, sort_order):
        """
        Sorts reactant and product matrices along with rate constants based on a specified sort order.
        If the given sort_order has entry that is not in rate_constant_names, an error will
        be raised. If the given sort_order does not contain all the entries in rate_constant_names,
        the output will be truncated based on the given sort_order.

        Warning:
            The output type after sorting however will be numpy array instead of
            python lists. Numpy array is used to enhance the efficiency of sorting.

        Args:
            reactant_matrix (numpy.ndarray): The matrix representing reactants.
            product_matrix (numpy.ndarray): The matrix representing products.
            rate_constant_names (list): List of rate constant names.
            sort_order (list): List of rate constant names in the desired sorting order.

        Returns:
            numpy.ndarray: Sorted reactant matrix.
            numpy.ndarray: Sorted product matrix.
            numpy.ndarray: Sorted rate constants.

        Example:
            Suppose you have the following input matrices and rate constant names:

            reactant_matrix = np.array([[1, 2, 3],
                                       [4, 5, 6]])

            product_matrix = np.array([[7, 8, 9],
                                       [10, 11, 12]])

            rate_constant_names = ['ki', 'kon', 'koff', 'kf']

            sort_order = ['kon', 'koff', 'kf', 'ki']

            # Create an instance of YourClassName (assuming it's a class method)
            instance = YourClassName()

            # Call the sort_by_rate_constants method
            sorted_reactant, sorted_product, sorted_constants = instance.sort_by_rate_constants(
                reactant_matrix, product_matrix, rate_constant_names, sort_order)

            # The function will return the sorted matrices and constants as follows:
            # sorted_reactant = np.array([[2, 3, 1],
            #                            [5, 6, 4]])

            # sorted_product = np.array([[8, 9, 7],
            #                           [11, 12, 10]])

            # sorted_constants = np.array(['kon', 'koff', 'kf', 'ki'])

        Note:
            This function assumes that the input matrices and rate constant names
            are compatible and correctly ordered.
        """
        # Convert rate_constant_names to a NumPy array
        rate_constant_names = np.array(rate_constant_names)

        # Get the indices to sort the rate constants
        sorted_indices = [list(rate_constant_names).index(name)
                          for name in sort_order]

        # Commandline print if debug mode enabled
        if self.DEBUG_MODE:
            print(sorted_indices)

        # Sort the rate constants, reactant_matrix, and product_matrix based on the sorted indices
        sorted_rate_constants = rate_constant_names[sorted_indices]
        sorted_reactant_matrix = reactant_matrix[sorted_indices, :]
        sorted_product_matrix = product_matrix[sorted_indices, :]

        return sorted_reactant_matrix, sorted_product_matrix, sorted_rate_constants

    def sort_by_species_names(self, reactant_matrix, product_matrix,
                              species_names, sort_order):
        """
        Sorts reactant and product matrices along with species names based on a specified sort order.

        Args:
            reactant_matrix (numpy.ndarray): The matrix representing reactants.
            product_matrix (numpy.ndarray): The matrix representing products.
            species_names (list): List of species names.
            sort_order (string, list, or callable):
                - If a string, options include "alphabetical," "increasing," "decreasing,"
                  "alphabetical" and "increasing" both use ASCII sequence and "decreasing"
                  uses reverse ASCII sequence.
                - If a list, provide a custom list of species names for sorting.
                - If a callable, sorter_order specifies a key function for custom sorting.
                  see python sorted() command for more information.


        Returns:
            numpy.ndarray: Sorted reactant matrix.
            numpy.ndarray: Sorted product matrix.
            numpy.ndarray: Sorted species names.

        Example:
            Suppose you have the following input matrices and species names:

            reactant_matrix = np.array([[1, 2, 3],
                                    [4, 5, 6]])

            product_matrix = np.array([[7, 8, 9],
                                    [10, 11, 12]])

            species_names = ['A', 'B', 'C', 'D']

            sort_order = ['B', 'A', 'D', 'C']

            # Create an instance of YourClassName (assuming it's a class method)
            instance = YourClassName()

            # Call the sort_by_species_names method
            sorted_reactant, sorted_product, sorted_species = instance.sort_by_species_names(
                reactant_matrix, product_matrix, species_names, sort_order)

            # The function will return the sorted matrices and species names as follows:
            # sorted_reactant = np.array([[2, 1, 3],
            #                            [5, 4, 6]])

            # sorted_product = np.array([[8, 7, 9],
            #                           [11, 10, 12]])

            # sorted_species = np.array(['B', 'A', 'D', 'C'])

        Note:
            This function assumes that the input matrices and species names are compatible and correctly ordered.
        """
        # Convert rate_constant_names to a NumPy array
        species_names = np.array(species_names)

        # If the given sort order is a string, generate a sort order list
        # based on the type of sequence given
        if isinstance(sort_order, str):
            if sort_order.casefold() == "alphabetical" or sort_order.casefold() == "increasing":
                sort_order = sorted(species_names)
            elif sort_order.casefold() == "decreasing":
                sort_order = sorted(species_names, reverse=True)
        elif callable(sort_order):
            sort_order = sorted(species_names, key=sort_order)

        # Get the indices to sort the rate constants
        sorted_indices = [list(species_names).index(name)
                          for name in sort_order]

        # Commandline print if debug mode enabled
        if self.DEBUG_MODE:
            print(sorted_indices)

        # Sort the rate constants, reactant_matrix, and product_matrix based on the sorted indices
        sorted_species_names = species_names[sorted_indices]
        sorted_reactant_matrix = reactant_matrix[:, sorted_indices]
        sorted_product_matrix = product_matrix[:, sorted_indices]

        return sorted_reactant_matrix, sorted_product_matrix, sorted_species_names
    
    def reset_rate_constant_autonumbering(self):
        """Reset the autonumbering index for rate constant names.

        This method resets the internal index used for autonumbering rate constant names
        in the ReactionStringParser class. The autonumbering index is used to generate
        unique names for rate constants when they are not provided in the input.

        Example:
            parser = ReactionStringParser()
            parser.reset_rate_constant_autonumbering()

        Note:
            After calling this method, the next autonumbered rate constant name generated
            will start from 0.

        """
        self.__RATE_CONSTANT_NAME_INDEX = 0

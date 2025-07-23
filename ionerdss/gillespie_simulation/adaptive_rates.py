import numpy as np
def adaptive_bimolecular_rate_1D(
        ratelist, y_curr, reactant_matrix:np.ndarray, Length,
        diffusion_constants, sigmalist, reverse_reaction_pairs:dict[int,int],
    ):
    """
    Update the bimolecular macroscopic rate.
        k_on = ( 1/k_a + (b-sigma)/3D )^-1
    Here ka is the microscopic rate; b is the mean-field space taken by the reactant 
    with more copy numbers; sigma is the bond length; and D is the total diffusion 
    constant = D1 + D2
    
    Args:
        ratelist: it contains both ka in nm/s and kb in s^-1
        y_curr: current copy numbers
        reactant_matrix (numpy.ndarray): Matrix representing reactants in each reaction.
        Length: length in nm
        diffusion_constants (array like): diffusion constant of each species
        sigmas (array like): sigma of each reaction. Dissociation may take sigma=0. 
        reverse_reaction_pairs (dict): pair forward reaction to its reverse reaction in
            the form of {forward_index:reverse_index}

    Returns:
        kon: Macroscopic rate constants

    Usage:
    ```python
        import numpy as np

        # For dimerization
        ratelist = np.array([10, 0.2])
        y_curr = np.array([10, 3])
        reactant_matrix = np.array([[2, 0], [0,1]])
        Length = 50
        diffusion_constants = np.array([0.1, 0.05])
        sigmalist = np.array([1, 0])
        reverse_reaction_pairs = {0:1}
        print(adaptive_bimolecular_rate_1D(
            ratelist, y_curr, reactant_matrix, Length,
            diffusion_constants, sigmalist, reverse_reaction_pairs,
        ))

        # For other bimolecular reactions
        ratelist = np.array([10, 0.2])
        y_curr = np.array([5, 5, 3])
        reactant_matrix = np.array([[1, 1, 0], [0, 0, 1]])
        Length = 50
        diffusion_constants = np.array([0.1, 0.1, 0.05])
        sigmalist = np.array([1, 0])
        reverse_reaction_pairs = {0:1}
        print(adaptive_bimolecular_rate_1D(
            ratelist, y_curr, reactant_matrix, Length,
            diffusion_constants, sigmalist, reverse_reaction_pairs,
        ))
    ```
    """

    # san check:
    # Species and diffusion constants are one-by-one mapped
    if len(diffusion_constants) != reactant_matrix.shape[1]:
        raise ValueError(f'diffusion_constants size ({len(diffusion_constants)}) does not match the number of species ({reactant_matrix.shape[1]})')
    # Reactions and sigmas are one-by-one mapped
    if len(sigmalist) != reactant_matrix.shape[0]:
        raise ValueError(f'sigmas size ({len(sigmalist)}) does not match the number of reactions ({reactant_matrix.shape[0]})')
    assert len(ratelist)==reactant_matrix.shape[0], f"ratelist size ({len(ratelist)}) does not match the number of reactions ({reactant_matrix.shape[0]})"
    assert len(y_curr)==reactant_matrix.shape[1], f"y_curr size ({len(y_curr)}) does not match the number of species ({reactant_matrix.shape[1]})"
    
    # initiate kon list
    new_ratelist = np.zeros_like(ratelist)
    # calculate kon
    for reactionid, row in enumerate(reactant_matrix):
        # Find positions where the value is 1
        reactant_ids = np.where(row == 1)[0]
        dimerization_id = np.where(row == 2)[0]
        # Find bimolecular reactions
        if len(reactant_ids) == 2:
            # find the spicies with more copy numbers
            species_id_more_counts = reactant_ids[np.argmax(y_curr[reactant_ids])]
            b = Length / y_curr[species_id_more_counts]
            D_tot = np.sum([diffusion_constants[i] for i in reactant_ids])
            sigma = sigmalist[reactionid]
            ka = ratelist[reactionid]
            new_ratelist[reactionid] = ( 1/ka + (b-sigma)/3/D_tot )**(-1)
            reverse_reaction_id = reverse_reaction_pairs[reactionid]
            kb = ratelist[reverse_reaction_id]
            new_ratelist[reverse_reaction_id] = kb * new_ratelist[reactionid] / ka
        # dimerization is a special case
        elif len(dimerization_id) == 1:
            b = Length / y_curr[dimerization_id[0]]
            D_tot = 2 * diffusion_constants[dimerization_id[0]]
            sigma = sigmalist[reactionid]
            ka = ratelist[reactionid]
            new_ratelist[reactionid] = ( 1/ka + (b-sigma)/3/D_tot )**(-1)
            reverse_reaction_id = reverse_reaction_pairs[reactionid]
            kb = ratelist[reverse_reaction_id]
            new_ratelist[reverse_reaction_id] = kb * new_ratelist[reactionid] / ka
        else:
            if new_ratelist[reverse_reaction_id] == 0 and ratelist[reactionid] != 0:
                new_ratelist[reverse_reaction_id] = ratelist[reactionid]
    # validate whether all rates are updated
    missed_rates_id = []
    info = 'reaction_id, old_rate, new_rate\n'
    for reactionid, (new_rate, old_rate) in enumerate(zip(new_ratelist, ratelist)):
        if new_rate == 0 and old_rate != 0:
            missed_rates_id.append(reactionid)
            info += f'{reactionid:>11}, {old_rate:>3.2e}, {new_rate:>3.2e}\n'
    if len(missed_rates_id) > 0:
        raise ValueError(info)
    return new_ratelist

if __name__ == '__main__':
    import numpy as np
    ratelist = np.array([10, 0.2])
    y_curr = np.array([10, 3])
    reactant_matrix = np.array([[2, 0], [0,1]])
    Length = 50
    diffusion_constants = np.array([0.1, 0.05])
    sigmalist = np.array([1, 0])
    reverse_reaction_pairs = {0:1}
    print(adaptive_bimolecular_rate_1D(
        ratelist, y_curr, reactant_matrix, Length,
        diffusion_constants, sigmalist, reverse_reaction_pairs,
    ))

    ratelist = np.array([10, 0.2])
    y_curr = np.array([5, 5, 3])
    reactant_matrix = np.array([[1, 1, 0], [0, 0, 1]])
    Length = 50
    diffusion_constants = np.array([0.1, 0.1, 0.05])
    sigmalist = np.array([1, 0])
    reverse_reaction_pairs = {0:1}
    print(adaptive_bimolecular_rate_1D(
        ratelist, y_curr, reactant_matrix, Length,
        diffusion_constants, sigmalist, reverse_reaction_pairs,
    ))

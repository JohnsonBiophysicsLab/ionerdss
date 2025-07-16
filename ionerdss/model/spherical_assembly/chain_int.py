# This function will go over every atom between two chains to determine whether they are interacting (distance smaller
# than 3.5A)
# remember to import math package when use the function
# Input variables:
# return variables: a tuple includes


def chain_int(unique_chain, split_position, split_resi_count, split_atom_count, split_resi_type, split_atom_type, split_resi_position):
    import math
    distance = 0
    reaction_chain = [] # list of lists (each sublist will include two letters indicating these two chains have
    # interaction) eg: in this protein, only chain A&B, A&D and C&D are interacting, then the list will look like
    # [[A,B],[A,D],[C,D]]
    reaction_resi_type = [] # list of lists of lists(each sub-sublist will include a bunch of lists of residue pairs
    # (without repeats)) eg: [[[resia,resib],[resic,resid]],[[resie,resif],[resig,resih]],[[resii,resij],[resik,resil]]]
    # ----reaction residues of chain-------- A&B------------------------A&D-------------------------C&D -------------
    reaction_resi_count = []
    reaction_atom = []
    reaction_atom_position = []
    reaction_atom_distance = []
    reaction_atom_type = []
    reaction_resi_position = []

    for i in range(len(unique_chain) - 1):
        for j in range(i+1, len(unique_chain)):
            inner_atom_position = []
            inner_atom_distance = []
            inner_atom = []
            inner_reaction_resi_count = []
            inner_reaction_resi_type = []
            inner_reaction_atom_type = []
            inner_reaction_resi_position = []

            for m in range(len(split_position[i])):
                for n in range(len(split_position[j])):
                    distance = math.sqrt((split_position[i][m][0]-split_position[j][n][0])**2
                                         + (split_position[i][m][1]-split_position[j][n][1])**2
                                         + (split_position[i][m][2]-split_position[j][n][2])**2)
                    if distance <= 0.85:
                        inner_atom.append([split_atom_count[i][m], split_atom_count[j][n]])
                        inner_atom_distance.append(distance)
                        inner_atom_position.append([split_position[i][m], split_position[j][n]])
                        inner_reaction_atom_type.append([split_atom_type[i][m], split_atom_type[j][n]])
                        if [split_resi_count[i][m], split_resi_count[j][n]] not in inner_reaction_resi_count:
                            inner_reaction_resi_count.append([split_resi_count[i][m], split_resi_count[j][n]])
                            inner_reaction_resi_position.append([split_resi_position[i][m], split_resi_position[j][n]])
                            inner_reaction_resi_type.append([split_resi_type[i][m], split_resi_type[j][n]])

            if len(inner_reaction_resi_count) > 0:
                reaction_chain.append([unique_chain[i], unique_chain[j]])
                reaction_resi_count.append(inner_reaction_resi_count)
                reaction_resi_type.append(inner_reaction_resi_type)
                reaction_atom.append(inner_atom)
                reaction_atom_position.append(inner_atom_position)
                reaction_atom_distance.append(inner_atom_distance)
                reaction_atom_type.append(inner_reaction_atom_type)
                reaction_resi_position.append(inner_reaction_resi_position)
    return reaction_chain, reaction_atom, reaction_atom_position, reaction_atom_distance, reaction_resi_count, \
           reaction_resi_type, reaction_atom_type, reaction_resi_position







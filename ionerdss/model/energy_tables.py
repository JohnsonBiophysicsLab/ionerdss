"""
Energy tables and related utilities for PDB model processing.

This module contains energy tables and functions for calculating interaction
energies between residues and molecular components.
"""


def get_default_energy_table():
    """
    Return energy table for residue-residue interactions.

    Reference:
        Miyazawa, S., & Jernigan, R. L. (1996). Residue-residue potentials 
        with a favorable contact pair term and an unfavorable high packing density term,
        for simulation and threading. J Mol Biol, 256(3), 623–644.

    Returns:
        dict: A symmetric dictionary with residue pair tuples as keys and contact energies (in RT units) as values.
    """
    residues = [
        'CYS', 'MET', 'PHE', 'ILE', 'LEU', 'VAL', 'TRP', 'TYR', 'ALA', 'GLY',
        'THR', 'SER', 'ASN', 'GLN', 'ASP', 'GLU', 'HIS', 'ARG', 'LYS', 'PRO'
    ]

    # Extracted from the upper triangle of the table (manually transcribed)
    energy_matrix = [
        [-5.44],
        [-4.99, -5.46],
        [-5.80, -6.56, -7.26],
        [-5.50, -6.02, -6.84, -6.54],
        [-5.83, -6.41, -7.28, -7.04, -7.37],
        [-4.96, -5.32, -6.29, -6.05, -6.48, -5.52],
        [-4.95, -5.55, -6.16, -5.78, -6.14, -5.18, -5.06],
        [-4.16, -4.91, -5.66, -5.25, -5.67, -4.62, -4.66, -4.17],
        [-3.57, -3.94, -4.81, -4.58, -4.91, -4.04, -3.82, -3.36, -2.72],
        [-3.16, -3.39, -4.13, -3.78, -4.16, -3.38, -3.42, -3.01, -2.31, -2.24],
        [-3.11, -3.51, -4.28, -4.03, -4.34, -3.46, -3.22, -3.01, -2.32, -2.08, -2.12],
        [-2.86, -3.03, -4.02, -3.52, -3.92, -3.05, -2.99, -2.78, -2.01, -1.82, -1.96, -1.67],
        [-2.59, -2.95, -3.75, -3.24, -3.74, -2.83, -3.07, -2.76, -1.84, -1.74, -1.88, -1.58, -1.68],
        [-2.85, -3.30, -4.10, -3.67, -4.04, -3.07, -3.11, -2.97, -1.89, -1.66, -1.90, -1.49, -1.71, -1.54],
        [-2.41, -2.57, -3.48, -3.17, -3.40, -2.48, -2.84, -2.76, -1.70, -1.59, -1.80, -1.63, -1.68, -1.46, -1.21],
        [-2.27, -2.89, -3.56, -3.27, -3.59, -2.67, -2.99, -2.79, -1.51, -1.22, -1.74, -1.48, -1.51, -1.42, -1.02, -0.91],
        [-3.60, -3.98, -4.77, -4.14, -4.54, -3.58, -3.98, -3.52, -2.41, -2.15, -2.42, -2.11, -2.08, -1.98, -2.32, -2.15, -3.05],
        [-2.57, -3.12, -3.98, -3.63, -4.03, -3.07, -3.41, -3.16, -1.83, -1.72, -1.90, -1.62, -1.64, -1.80, -2.29, -2.27, -2.16, -1.55],
        [-1.95, -2.48, -3.36, -3.01, -3.37, -2.49, -2.69, -2.60, -1.31, -1.15, -1.31, -1.05, -1.21, -1.29, -1.68, -1.80, -1.35, -0.59, -0.12],
        [-3.07, -3.45, -4.25, -3.76, -4.20, -3.32, -3.73, -3.19, -2.03, -1.87, -1.90, -1.57, -1.53, -1.73, -1.33, -1.26, -2.25, -1.70, -0.97, -1.75]
    ]

    energy_table = {}

    for i, res_i in enumerate(residues):
        for j, res_j in enumerate(residues[:i+1]):
            energy = energy_matrix[i][j] + 2.27  # Adjusted energy value
            energy_table[(res_i, res_j)] = energy
            energy_table[(res_j, res_i)] = energy  # symmetry

    return energy_table
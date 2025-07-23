"""
energy_table.py

@deprecated This method does not work well in estimateing the binding free energy.

Provides a residue-residue interaction energy lookup table used in the coarse-graining
pipeline for identifying and scoring protein-protein interfaces.

This table serves as a simplified, empirical scoring function that estimates the 
binding affinity of residue pairs based on their types. It is primarily used during 
interface detection to assign a total "interaction energy" score to each potential 
interface between chains.

The energy values are heuristic and intended for qualitative comparison, not physical 
accuracy. Typical sources of such scoring schemes include knowledge-based potentials 
(e.g., Miyazawa–Jernigan matrix) or parameter tuning based on known structures.

Functions
---------
- get_default_energy_table(): Returns a dictionary of (res1, res2) → energy (float) pairs.
                              The table is symmetric: E(A, B) == E(B, A)
"""

def get_default_energy_table():
    """
    Returns a dictionary mapping residue pairs to interaction energies (kcal/mol).

    The returned dictionary is symmetric: both ('ARG', 'GLU') and ('GLU', 'ARG') will be present.

    Example:
        {
            ('ARG', 'GLU'): -3.5,
            ('GLU', 'ARG'): -3.5,
            ...
        }

    Returns
    -------
    dict of tuple(str, str) → float
        A lookup table of interaction energies between residue types.
    """
    raw_table = {
        ('ARG', 'GLU'): -3.5,
        ('ARG', 'ASP'): -3.2,
        ('LYS', 'GLU'): -3.0,
        ('LYS', 'ASP'): -2.8,
        ('HIS', 'GLU'): -2.0,
        ('HIS', 'ASP'): -1.8,
        ('ARG', 'ARG'): 0.5,
        ('LYS', 'LYS'): 0.5,
        ('GLU', 'GLU'): 0.5,
        ('ASP', 'ASP'): 0.5,
        ('LEU', 'ILE'): -2.5,
        ('LEU', 'VAL'): -2.4,
        ('ILE', 'VAL'): -2.3,
        ('PHE', 'TYR'): -2.2,
        ('PHE', 'PHE'): -2.1,
        ('TYR', 'TYR'): -2.0,
        ('TRP', 'TRP'): -1.8,
        ('SER', 'THR'): -1.2,
        ('ASN', 'GLN'): -1.0,
        ('ALA', 'ALA'): -0.5,
    }

    # Symmetrize the table
    symmetric_table = {}
    for (res1, res2), energy in raw_table.items():
        symmetric_table[(res1, res2)] = energy
        symmetric_table[(res2, res1)] = energy

    return symmetric_table

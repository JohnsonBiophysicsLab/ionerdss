"""
repeated_interface_signature
------------------------------------

This module provides utility functions to compute and compare geometric and sequence-based 
signatures for identifying homologous or repeated molecular binding interfaces. These 
signatures are used to detect and match structurally or sequentially equivalent interfaces 
across multiple protein chains or coarse-grained molecules.

Key functionalities include:

- `identify_interface_signature`: Computes a geometry-based signature of a binding interface
  between two chains using inter-residue distances and angles formed by centers of mass and
  interaction points.
  
- `identify_interface_sequence_signature`: Generates a sequence-based signature using the
  sorted residue names (e.g., amino acid types) involved in the interface.
  
- `signature_are_similar`: Compares two geometric interface signatures based on user-defined 
  thresholds for distances and angles.
  
- `signature_difference`: Calculates the relative difference between two geometric interface
  signatures.
  
- `build_signature`: Constructs a geometric signature dictionary given the COM and interface
  coordinates.
  
- `invert_signature`: Returns the conjugated form of a signature by swapping A ↔ B values.

These tools are intended for detecting repeated interfaces in symmetric molecular assemblies 
and for consolidating equivalent interfaces in coarse-grained models used in NERDSS.
"""

import numpy as np

from ionerdss.math.angles import angles_from_points
from ionerdss.math.coords import Coords

def identify_interface_signature(mol_a, mol_b,
                                 cg_model,
                                 interface_index_on_a = -1):
    """
    Computes a geometric interaction signature between two chains/molecules.

    Args:
        mol_a (str): Chain ID of molecule A.
        mol_b (str): Chain ID of molecule B.
        cg_model (dict): Coarse-grained structure data including COMs and interfaces.
        interface_index_on_a (int, optional): Index of the interface on molecule A. Defaults to -1 (last).

    Returns:
        tuple:
            - signature (dict): Geometric descriptor from A to B.
            - signature_conjugated (dict): Inverted signature from B to A.
            - I_B (Coords): Interface coordinate on molecule B.
            - k (int): Index of the matching interface on molecule B.
    """
    index_a = cg_model["chains"].index([chain for chain in cg_model["chains"] if chain.id == mol_a][0])
    index_b = cg_model["chains"].index([chain for chain in cg_model["chains"] if chain.id == mol_b][0])

    com_a = cg_model["COMs"][index_a]
    com_b = cg_model["COMs"][index_b]

    # if the interface id for the interactions are provided, use the input
    # to avoid enumerative recalculation of the interacting interfaces

    I_A = cg_model["interface_coords"][index_a][interface_index_on_a]

    # Find the matching interface on B that binds to A
    for k, partner_interface_id in enumerate(cg_model["interfaces"][index_b]):
        if partner_interface_id == mol_a:
            I_B = cg_model["interface_coords"][index_b][k]
            break
    else:
        raise ValueError(f"No matching interface found on {mol_b} for partner {mol_a}")

    signature = {
        "dA": np.linalg.norm([(com_a - I_A).x, (com_a - I_A).y, (com_a - I_A).z]),
        "dB": np.linalg.norm([(com_b - I_B).x, (com_b - I_B).y, (com_b - I_B).z]),
        "dAB": np.linalg.norm([(I_A - I_B).x, (I_A - I_B).y, (I_A - I_B).z]),
        "thetaA": angles_from_points(com_a, I_A, I_B),
        "thetaB": angles_from_points(com_b, I_B, I_A),
    }

    signature_conjugated = {
        "dA": signature["dB"],
        "dB": signature["dA"],
        "dAB": signature["dAB"],
        "thetaA": signature["thetaB"],
        "thetaB": signature["thetaA"],
    }

    return signature, signature_conjugated, I_B, k

def identify_interface_sequence_signature(A, B, i, cg_model):
    """
    Computes a sequence-based signature for an interaction using sorted residue types.

    Args:
        A (str): Chain ID of molecule A.
        B (str): Chain ID of molecule B.
        i (int): Index of the interface on molecule A.
        cg_model (dict): Coarse-grained structure data.

    Returns:
        tuple:
            - signature (dict): {'seqA': str, 'seqB': str}
            - signature_conjugated (dict): Reversed version of the sequence pair
            - residues_B (list): Residues on the matching interface on B
            - k (int): Index of the matching interface on B
    """
    index_a = cg_model["chains"].index([chain for chain in cg_model["chains"] if chain.id == A][0])
    index_b = cg_model["chains"].index([chain for chain in cg_model["chains"] if chain.id == B][0])

    residues_a = cg_model["interface_residues"][index_a][i]

    # Find matching interface on B
    for k, partner_interface_id in enumerate(cg_model["interfaces"][index_b]):
        if partner_interface_id == A:
            residues_b = cg_model["interface_residues"][index_b][k]
            break
    else:
        raise ValueError(f"No matching interface found on {B} for partner {A}")

    # Build signature based on sorted residue identities (amino acid types only)
    seq_a = ''.join(sorted([residue.get_resname().strip() for residue in residues_a]))
    seq_b = ''.join(sorted([residue.get_resname().strip() for residue in residues_b]))

    signature = {
        "seqA": seq_a,
        "seqB": seq_b,
    }

    signature_conjugated = {
        "seqA": seq_b,
        "seqB": seq_a,
    }

    return signature, signature_conjugated, residues_b, k

def signature_are_similar(sig1, sig2,
                     dist_threshold_intra,
                     dist_threshold_inter,
                     angle_threshold):
    """
    Compares two groups of interface interaction geometry signatures.

    Args:
        sig1 (dict): The first interface signature.
        sig2 (dict): The second interface signature.
        dist_threshold_intra (float): Distance threshold for intra-molecular comparisons.
        dist_threshold_inter (float): Distance threshold for inter-molecular comparisons.
        angle_threshold (float): Angle threshold for comparisons.

    Returns:
        bool: True if the signatures are similar within the given thresholds, False otherwise.
    """
    for key in ("dA", "dB"):
        if abs(sig1[key] - sig2[key]) > dist_threshold_intra:
            return False
    for key in ("dAB",):
        if abs(sig1[key] - sig2[key]) > dist_threshold_inter:
            return False
    for key in ("thetaA", "thetaB"):
        if abs(sig1[key] - sig2[key]) > angle_threshold:
            return False
    return True

def signature_difference(sig1, sig2):
    """
    Computes a relative distance between two signatures as a normalized sum of component-wise differences.

    Args:
        sig1 (dict): First signature.
        sig2 (dict): Second signature.

    Returns:
        float: Total normalized difference across all signature keys.
    """
    total_diff = 0.0
    for key in ("dA", "dB", "dAB", "thetaA", "thetaB"):
        val1 = sig1[key]
        val2 = sig2[key]

        if np.isnan(val1) or np.isnan(val2):
            return float("nan")

        denom = abs(val1) if abs(val1) > 1e-6 else 1.0
        total_diff += abs(val1 - val2) / denom
    return total_diff

def build_signature(COM_A, I_A, COM_B, I_B):
    """
    Builds a geometric interaction signature between two COM-interface pairs.

    Args:
        COM_A (Coords): Center of mass for molecule A.
        I_A (Coords): Interface point on molecule A.
        COM_B (Coords): Center of mass for molecule B.
        I_B (Coords): Interface point on molecule B.

    Returns:
        dict: {'dA', 'dB', 'dAB', 'thetaA', 'thetaB'} signature of interaction.
    """
    return {
        "dA": np.linalg.norm([(COM_A - I_A).x, (COM_A - I_A).y, (COM_A - I_A).z]),
        "dB": np.linalg.norm([(COM_B - I_B).x, (COM_B - I_B).y, (COM_B - I_B).z]),
        "dAB": np.linalg.norm([(I_A - I_B).x, (I_A - I_B).y, (I_A - I_B).z]),
        "thetaA": angles_from_points(COM_A, I_A, I_B),
        "thetaB": angles_from_points(COM_B, I_B, I_A)
    }

def invert_signature(sig):
    """
    Produces the conjugated version of an interaction signature (swap A <-> B).

    Args:
        sig (dict): Input signature.

    Returns:
        dict: Inverted version of the signature.
    """
    return {
        "dA": sig["dB"],
        "dB": sig["dA"],
        "dAB": sig["dAB"],
        "thetaA": sig["thetaB"],
        "thetaB": sig["thetaA"]
    }
    
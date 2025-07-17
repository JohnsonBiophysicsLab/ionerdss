"""
Data classes for PDB model processing.

This module contains the core data structures used in the PDB model processing
pipeline, including templates, molecules, interfaces, and reactions.
"""

import numpy as np
from Bio.PDB.Polypeptide import is_aa
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1
from scipy.spatial import KDTree
from sklearn.cluster import KMeans
import math


class MoleculeTemplate:
    """
    Template for a particular molecule type that can be reused across multiple complexes.

    Attributes:
        name (str): Identifier for the molecule type.
        interface_template_list (list): A list of BindingInterfaceTemplate objects that 
            describe the molecule's binding sites.
        normal_point (list): Default normal vector direction.
    """

    def __init__(self, name: str):
        """
        Initialize a MoleculeTemplate.

        Args:
            name (str): Name/identifier of the molecule type.
        """
        self.name = name
        self.interface_template_list = []
        self.normal_point = [0, 0, 1]
        self.diffusion_translation = None
        self.diffusion_rotation = None
        self.radius = None

    def __str__(self):
        interfaces = "\n  ".join(str(it) for it in self.interface_template_list)
        return f"Molecule Template: {self.name}\n  Interfaces:\n  {interfaces}"
    
    def __eq__(self, other):
        if not isinstance(other, MoleculeTemplate):
            return False
        return self.name == other.name


class BindingInterfaceTemplate:
    """
    Template for a binding interface between molecules.

    Attributes:
        name (str): Identifier of the interface template.
        coord (Coords): Relative coordinates of the interface.
        my_residues (list): Residues forming this interface.
        required_free_list (list): Other interface templates that must remain unbound 
            for this interface to bind.
        signature (dict): Stores interface geometry information.
    """

    def __init__(self, name: str):
        """
        Initialize a BindingInterfaceTemplate.

        Args:
            name (str): Identifier for the interface template.
        """
        self.name = name
        self.coord = None
        self.my_residues = []
        self.required_free_list = []  # The list of interface templates that need to be free to bind to this interface template
        self.signature = {}
        self.energy = None

    def __str__(self):
        residues = ", ".join(self.my_residues)
        required_free = ", ".join(self.required_free_list)
        return (f"Interface Template: {self.name}\n"
                f"  Coordinates: {self.coord}\n"
                f"  Residues: {residues}\n"
                f"  Required Free: {required_free}")
    
    def __eq__(self, other):
        if not isinstance(other, BindingInterfaceTemplate):
            return False
        return self.name == other.name


class CoarseGrainedMolecule:
    """
    Coarse-grained molecule in NERDSS, potentially derived from a PDB chain.

    Attributes:
        name (str): Identifier of the molecule.
        my_template (MoleculeTemplate): Reference to the associated molecule template.
        coord (Coords): Center-of-mass coordinates.
        interface_list (list): List of binding interfaces.
        normal_point (list): Normal vector direction.
    """

    def __init__(self, name: str):
        """
        Initialize a CoarseGrainedMolecule.

        Args:
            name (str): Name/identifier of the molecule.
        """
        self.name = name
        self.my_template = None
        self.coord = None
        self.interface_list = []
        self.normal_point = None
        self.diffusion_translation = None
        self.diffusion_rotation = None
        self.radius = None

    def __str__(self):
        interfaces = "\n  ".join(str(interface) for interface in self.interface_list)
        return (f"CoarseGrainedMolecule: {self.name}\n"
                f"  Template: {self.my_template}\n"
                f"  Coordinates: {self.coord}\n"
                f"  Interfaces:\n  {interfaces}")
    
    def __repr__(self):
        return self.name
    
    def __eq__(self, other):
        if not isinstance(other, CoarseGrainedMolecule):
            return False
        return self.name == other.name
    
    def __hash__(self):
        return hash(self.name)


class BindingInterface:
    """
    Binding interface between molecules.

    Attributes:
        name (str): Identifier of the binding interface.
        coord (Coords): Position of the interface.
        my_template (BindingInterfaceTemplate): Reference to the associated interface template.
        my_residues (list): Residues included in the interface.
        signature (dict): Stores interface geometry information.
    """

    def __init__(self, name: str):
        """
        Initialize a BindingInterface.

        Args:
            name (str): Identifier for the binding interface.
        """
        self.name = name
        self.coord = None
        self.my_template = None
        self.my_residues = []
        self.signature = {}
        self.energy = None

    def __str__(self):
        return (f"BindingInterface: {self.name}\n"
                f"  Template: {self.my_template}\n"
                f"  Coordinates: {self.coord}\n"
                f"  Residue Count: {len(self.my_residues)}\n"
                f"  Residues: {self.my_residues}")
    
    def __eq__(self, other):
        if not isinstance(other, BindingInterface):
            return False
        return self.my_template == other.my_template


class ReactionTemplate:
    """
    Template for a reaction between two MoleculeTemplates.

    Attributes:
        expression (str): Textual representation of the reaction.
        reactants (list): List of reactant molecule/interface templates.
        products (list): List of product molecule/interface templates.
        binding_angles (tuple): Tuple describing binding angles (theta1, theta2, phi1, phi2, omega).
        binding_radius (float): Distance between binding interfaces.
        norm1 (list): Normal vector of the first reactant.
        norm2 (list): Normal vector of the second reactant.
    """

    def __init__(self):
        """
        Initialize a ReactionTemplate with default values.
        """
        self.expression = None
        self.reactants = None
        self.products = None
        self.binding_angles = None
        self.binding_radius = None
        self.norm1 = None
        self.norm2 = None
        self.kd = None
        self.ka = None
        self.kb = None
        self.energy = None

    def __str__(self):
        return (f"Reaction Template: {self.expression}\n"
                f"  Reactants: {self.reactants}\n"
                f"  Products: {self.products}\n"
                f"  Binding Angles: {self.binding_angles}\n"
                f"  Binding Radius: {self.binding_radius / 10:.6f} nm\n"
                f"  norm1: {self.norm1}\n"
                f"  norm2: {self.norm2}")
    
    def __eq__(self, other):
        if not isinstance(other, ReactionTemplate):
            return False
        return self.expression == other.expression


class Reaction:
    """
    Actual reaction between two Molecule objects.

    Attributes:
        expression (str): Textual representation of the reaction.
        reactants (list): List of reactant molecules/interfaces.
        products (list): List of product molecules/interfaces.
        binding_angles (tuple): Tuple describing binding angles (theta1, theta2, phi1, phi2, omega).
        binding_radius (float): Distance between binding interfaces.
        norm1 (list): Normal vector of the first reactant.
        norm2 (list): Normal vector of the second reactant.
    """

    def __init__(self):
        """
        Initialize a Reaction with default values.
        """
        self.expression = None
        self.reactants = None
        self.products = None
        self.binding_angles = None
        self.binding_radius = None
        self.norm1 = None
        self.norm2 = None
        self.my_template = None
        self.kd = None
        self.ka = None
        self.kb = None
        self.energy = None

    def __str__(self):
        return (f"Reaction: {self.expression}\n"
                f"  Reactants: {self.reactants}\n"
                f"  Products: {self.products}\n"
                f"  Binding Angles: {self.binding_angles}\n"
                f"  Binding Radius: {self.binding_radius / 10:.6f} nm")

    def __repr__(self):
        return f"Reaction({self.expression})"

    def __eq__(self, other):
        if not isinstance(other, Reaction):
            return False
        return self.expression == other.expression
    
    def __hash__(self):
        return hash(self.expression)


# -------------------------------------------------------------------------
# Helper functions - geometry transformation
# -------------------------------------------------------------------------

def rigid_transform_3d(points_a: np.ndarray, points_b: np.ndarray):
    """
    Compute a rigid transformation (rotation + translation) that aligns 
    `points_a` to `points_b` using Singular Value Decomposition (SVD).

    Args:
        points_a (np.ndarray): Shape (N, 3), first set of 3D points.
        points_b (np.ndarray): Shape (N, 3), second set of 3D points.

    Returns:
        tuple:
            - np.ndarray: 3x3 rotation matrix `R`
            - np.ndarray: 3-element translation vector `t`
    """
    assert len(points_a) == len(points_b), "Point sets must be same length."
    centroid_a = points_a[0]
    centroid_b = points_b[0]
    pa = points_a[1:] - centroid_a
    pb = points_b[1:] - centroid_b
    h = pa.T @ pb
    u, s, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    t = centroid_b - r @ centroid_a
    return r, t


def apply_rigid_transform(r: np.ndarray, t: np.ndarray, point: np.ndarray):
    """
    Apply a rigid transformation (rotation + translation) to a point.

    Args:
        r (np.ndarray): A 3x3 rotation matrix.
        t (np.ndarray): A 3-element translation vector.
        point (np.ndarray): A shape (3,) array representing a point.

    Returns:
        np.ndarray: Transformed point.
    """
    return (r @ point.T).T + t


def rigid_transform_chains(chain1, chain2):
    """
    Align chain1 to chain2 by:
    1. Extracting amino acid sequences.
    2. Performing sequence alignment.
    3. Identifying matching residues.
    4. Computing a coarse-grained set of representative points.
    5. Computing a rigid transformation.

    Args:
        chain1 (Bio.PDB.Chain.Chain): First molecular chain.
        chain2 (Bio.PDB.Chain.Chain): Second molecular chain.

    Returns:
        tuple:
            - np.ndarray: 3x3 rotation matrix `R`
            - np.ndarray: 3-element translation vector `t`
    """

    # Step 1: Extract sequences from both chains
    def extract_sequence(chain):
        """Extract the amino acid sequence from a chain."""
        return "".join(seq1(residue.resname) for residue in chain.get_residues() if is_aa(residue))

    sequence1 = extract_sequence(chain1)
    sequence2 = extract_sequence(chain2)

    # Step 2: Find the best overlap between the two sequences using PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5
    
    alignments = aligner.align(sequence1, sequence2)
    alignment = alignments[0]  # Get the best alignment

    aligned_seq1 = alignment[0]
    aligned_seq2 = alignment[1]

    # Step 3: Identify matching residue pairs in the aligned sequences
    residue_pairs = []
    idx1, idx2 = 0, 0
    residues1 = [res for res in chain1 if is_aa(res)]
    residues2 = [res for res in chain2 if is_aa(res)]

    for a1, a2 in zip(aligned_seq1, aligned_seq2):
        if a1 == '-' or a2 == '-':
            if a1 != '-':
                idx1 += 1
            if a2 != '-':
                idx2 += 1
            continue
        residue_pairs.append((residues1[idx1]['CA'].coord, residues2[idx2]['CA'].coord))
        idx1 += 1
        idx2 += 1

    # Step 4: Group residues into four spatially groups
    def group_residues(residues, n_groups=4):
        """Group residues into n_groups based on their spatial proximity."""
        coords = np.array([res for res, _ in residues])
        kmeans = KMeans(n_clusters=n_groups).fit(coords)
        groups = [[] for _ in range(n_groups)]
        for i, label in enumerate(kmeans.labels_):
            groups[label].append(residues[i])
        return groups

    groups = group_residues(residue_pairs)

    # Step 5: Compute the average position of each group and COM
    P = [np.mean([res[0] for res in group], axis=0) for group in groups]
    Q = [np.mean([res[1] for res in group], axis=0) for group in groups]
    P.insert(0, np.mean([res[0] for res in residue_pairs], axis=0))
    Q.insert(0, np.mean([res[1] for res in residue_pairs], axis=0))

    P = np.array(P)
    Q = np.array(Q)

    # Step 6: Apply rigid transformation
    R, t = rigid_transform_3d(P, Q)

    return R, t


def check_steric_clashes(points_1, points_2, cutoff: float = 3.5, number_threshold: int = 2):
    """
    Detect steric clashes between two sets of molecular points.

    Args:
        points_1 (np.ndarray): N x 3 coordinates for the first molecule.
        points_2 (np.ndarray): M x 3 coordinates for the second molecule.
        cutoff (float, optional): Distance threshold (default: 3.5 Å).
        number_threshold (int, optional): Minimum number of close contacts to flag a clash (default: 2).

    Returns:
        bool: True if a steric clash is detected, False otherwise.
    """
    tree = KDTree(points_2)
    clashes = tree.query_ball_point(points_1, r=cutoff)
    return any(len(clash) >= number_threshold for clash in clashes)
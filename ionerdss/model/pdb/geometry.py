"""
geometry.py

Provides geometric transformation utilities used throughout the NERDSS coarse-graining
pipeline, including 3D rigid alignment of repeated chains, vector transformations,
steric clash detection, and angle measurements.

Key Features
------------
- Rigid-body alignment using singular value decomposition (SVD), ensuring optimal
  rotation and translation that minimizes RMSD between two point clouds.
- Transformation application for atomic or coarse-grained coordinates.
- Utility functions for angular analysis between unit vectors or planes.
- Steric clash detection between coarse-grained spheres.

Functions
---------
- rigid_transform_chains(chain1_coords, chain2_coords): Convenience wrapper for rigid_transform_3d.
- check_steric_clashes(pos1, pos2, r1, r2, buffer=0.0): Returns True if two spheres overlap.
- angle_between(v1, v2): Computes angle (in radians) between two 3D vectors.
- dihedral(p1, p2, p3, p4): Returns dihedral angle (in radians) defined by 4 points.
"""

import numpy as np
from Bio.PDB import is_aa
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1
from scipy.cluster.vq import kmeans, vq

from ionerdss.utils.rigid_transform import rigid_transform_3d

def rigid_transform_chains(chain1, chain2, n_cluster_groups = 4):
    """
    Aligns chain1 to chain2 by:
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
        """Extracts the amino acid sequence from a chain."""
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
    def group_residues(residues, n_cluster_groups=4):
        """Groups residues into n_cluster_groups based on their spatial proximity."""

        # get coordinates of residues
        coords = np.array([res for res, _ in residues])

        # Adjust number of groups if not enough data
        n_from_points = len(coords) - 1
        if n_from_points < 1:
            raise ValueError("At least two residues must be present per chain for alignment.")
        n_clusters = min(n_cluster_groups, n_from_points - 1)

        if n_clusters < 1:
            raise ValueError("Number of cluster groups cannot be fewer than 1 for KMeans clustering.")

        # `coords` should be a NumPy array of shape (N, D)
        centroids, _ = kmeans(coords, n_cluster_groups)
        labels, _ = vq(coords, centroids)

        groups = [[] for _ in range(n_cluster_groups)]
        for i, label in enumerate(labels):
            groups[label].append(residues[i])
        return groups

    groups = group_residues(residue_pairs, n_cluster_groups=n_cluster_groups)

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


def check_steric_clashes(pos1, pos2, r1, r2, buffer=0.0):
    """
    Returns True if two coarse-grained spheres overlap beyond allowed buffer.

    Parameters
    ----------
    pos1 : array-like of shape (3,)
        Center of first molecule.
    pos2 : array-like of shape (3,)
        Center of second molecule.
    r1 : float
        Radius of first molecule.
    r2 : float
        Radius of second molecule.
    buffer : float, optional
        Extra distance allowed without clash. Default is 0.0 nm.

    Returns
    -------
    bool
        True if spheres are clashing.
    """
    d = np.linalg.norm(np.asarray(pos1) - np.asarray(pos2))
    return d < (r1 + r2 - buffer)

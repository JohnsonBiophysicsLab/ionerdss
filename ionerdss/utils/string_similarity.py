"""
ionerdss.math.string_similarity

This module provides utility functions to compute:
- the Levenshtein distance and normalized similarity between two strings.
- the Jaccard similarity

Levenshtein distance:

Purpose:
--------
The Levenshtein distance is a measure of how many single-character edits
(insertions, deletions, or substitutions) are required to change one string
into another. It is useful in sequence comparison tasks such as:
- Interface residue sequence matching
- Fuzzy identity clustering
- Sequence-based template assignment

Method:
-------
We implement the classical dynamic programming algorithm using a
(len(s1)+1) × (len(s2)+1) matrix, where each entry represents the edit distance
between substrings.

Dependencies:
-------------
- numpy (used for efficient array operations)

Usage:
------
>>> levenshtein_distance("ABC", "ADC")
1

>>> levenshtein_similarity("ABC", "ADC")
0.666...

>>> levenshtein_similarity("AAAA", "AAAA")
1.0

>>> levenshtein_similarity("AAAA", "")
0.0

Jaccard Similarity:

Purpose:
--------
Jaccard similarity is a simple and efficient measure of set overlap, defined as:

    J(A, B) = |A ∩ B| / |A ∪ B|

This is useful in molecular modeling when:
- Comparing interface residues where residue **types** matter more than order
- Handling **deletions, insertions**, or **partial interfaces**
- Detecting approximate conservation of **binding site chemistry**

Unlike Levenshtein distance, which penalizes for position-dependent edits,
Jaccard similarity is **order-insensitive** and tolerant to gaps.

Dependencies:
-------------
- numpy (only used for type checking, can be removed if undesired)

Usage:
------
>>> jaccard_similarity("ADEFG", "ADEFG")
1.0

>>> jaccard_similarity("ADEFG", "ADCFG")
0.6

>>> jaccard_similarity("ABCDE", "ACDF")
0.5

>>> jaccard_similarity("", "ABC")
0.0

>>> jaccard_similarity("ABC", "")
0.0

>>> jaccard_similarity("", "")
1.0

"""

import numpy as np

def levenshtein_distance(s1, s2):
    """
    Compute the Levenshtein (edit) distance between two strings.

    The distance represents the minimum number of single-character edits
    (insertions, deletions, or substitutions) required to transform s1 into s2.

    Args:
        s1 (str): First string.
        s2 (str): Second string.

    Returns:
        int: Levenshtein edit distance.
    """
    len_s1, len_s2 = len(s1), len(s2)

    # Initialize (len_s1+1) × (len_s2+1) matrix
    dp = np.zeros((len_s1 + 1, len_s2 + 1), dtype=int)

    # Base cases: distance from empty string
    dp[:, 0] = np.arange(len_s1 + 1)  # Cost of deleting all characters from s1
    dp[0, :] = np.arange(len_s2 + 1)  # Cost of inserting all characters into s1

    # Fill dynamic programming matrix
    for i in range(1, len_s1 + 1):
        for j in range(1, len_s2 + 1):
            if s1[i - 1] == s2[j - 1]:
                cost = 0  # No edit needed
            else:
                cost = 1  # Substitution required

            dp[i, j] = min(
                dp[i - 1, j] + 1,      # Deletion
                dp[i, j - 1] + 1,      # Insertion
                dp[i - 1, j - 1] + cost  # Substitution
            )

    return dp[len_s1, len_s2]

def levenshtein_similarity(s1, s2):
    """
    Compute the normalized Levenshtein similarity between two strings.

    Similarity is defined as:
        1 - (edit_distance / max(len(s1), len(s2)))

    This yields a value between 0.0 (completely different) and 1.0 (identical).

    Args:
        s1 (str): First string.
        s2 (str): Second string.

    Returns:
        float: Similarity score in the range [0.0, 1.0].
    """
    if not s1 and not s2:
        return 1.0  # Both strings empty → identical
    max_len = max(len(s1), len(s2))
    distance = levenshtein_distance(s1, s2)
    return 1.0 - (distance / max_len)

def jaccard_similarity(seq1, seq2):
    """
    Compute the Jaccard similarity between two sequences by treating them as sets
    of amino acid residue types (characters).

    Args:
        seq1 (str or iterable): First residue sequence.
        seq2 (str or iterable): Second residue sequence.

    Returns:
        float: Similarity score in the range [0.0, 1.0].
    """
    set1 = set(seq1)
    set2 = set(seq2)

    if not set1 and not set2:
        return 1.0  # Both empty → identical

    intersection = set1 & set2
    union = set1 | set2

    return len(intersection) / len(union) if union else 1.0

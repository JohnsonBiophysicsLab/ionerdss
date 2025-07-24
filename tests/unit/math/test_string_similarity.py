"""
tests.unit.math.test_string_similarity

Unit tests for ionerdss.math.string_similarity module.

This test suite verifies the correctness of:

1. Levenshtein distance:
   - Correct edit distance between common string pairs
   - Behavior with empty strings and identical strings

2. Normalized Levenshtein similarity:
   - Values between [0, 1]
   - Exact matches give 1.0
   - Disjoint strings give 0.0

3. Jaccard similarity:
   - Based on set overlap of character tokens
   - Order-insensitive and tolerant to insertions/deletions
   - Correct treatment of empty sets

These tests ensure robust behavior in interface residue comparison and fuzzy sequence matching workflows.

"""

import unittest
from ionerdss.math.string_similarity import (
    levenshtein_distance,
    levenshtein_similarity,
    jaccard_similarity,
)

class TestStringSimilarity(unittest.TestCase):

    def test_levenshtein_distance(self):
        self.assertEqual(levenshtein_distance("kitten", "sitting"), 3)
        self.assertEqual(levenshtein_distance("flaw", "lawn"), 2)
        self.assertEqual(levenshtein_distance("abc", "abc"), 0)
        self.assertEqual(levenshtein_distance("", "abc"), 3)
        self.assertEqual(levenshtein_distance("abc", ""), 3)
        self.assertEqual(levenshtein_distance("", ""), 0)

    def test_levenshtein_similarity(self):
        self.assertAlmostEqual(levenshtein_similarity("kitten", "sitting"), 1 - 3/7)
        self.assertAlmostEqual(levenshtein_similarity("abc", "abc"), 1.0)
        self.assertAlmostEqual(levenshtein_similarity("abc", ""), 0.0)
        self.assertAlmostEqual(levenshtein_similarity("", "abc"), 0.0)
        self.assertAlmostEqual(levenshtein_similarity("", ""), 1.0)

    def test_jaccard_similarity(self):
        self.assertAlmostEqual(jaccard_similarity("abc", "abc"), 1.0)
        self.assertAlmostEqual(jaccard_similarity("abc", "abd"), 2/4)
        self.assertAlmostEqual(jaccard_similarity("abc", "xyz"), 0.0)
        self.assertAlmostEqual(jaccard_similarity("", "abc"), 0.0)
        self.assertAlmostEqual(jaccard_similarity("abc", ""), 0.0)
        self.assertAlmostEqual(jaccard_similarity("", ""), 1.0)
        self.assertAlmostEqual(jaccard_similarity("aabcc", "abccc"), 1.0)  # sets: {a,b,c}

if __name__ == '__main__':
    unittest.main()

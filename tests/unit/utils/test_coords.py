import unittest
import numpy as np
from ionerdss.utils.coords import Coords, get_perpendicular_vector

class TestCoords(unittest.TestCase):

    def test_addition_and_subtraction(self):
        a = Coords(1, 2, 3)
        b = Coords(4, 5, 6)
        self.assertEqual(a + b, Coords(5, 7, 9))
        self.assertEqual(b - a, Coords(3, 3, 3))

    def test_distance_and_distance_squared(self):
        a = Coords(0, 0, 0)
        b = Coords(3, 4, 0)
        self.assertAlmostEqual(a.distance(b), 5.0)
        self.assertEqual(a.distance_squared(b), 25)

    def test_str_and_repr(self):
        a = Coords(1.2345, 6.7890, -3.210)
        # Python’s default rounding for "{:.3f}".format(1.2345) 
        # is not rounding up; it gives "1.234".
        self.assertEqual(str(a), "(1.234, 6.789, -3.210)")
        self.assertEqual(repr(a), "Coords(x=1.234, y=6.789, z=-3.210)")

    def test_equality(self):
        a = Coords(1.0000001, 2.0, 3.0)
        b = Coords(1.0, 2.0, 3.0)
        self.assertTrue(a == b)

    def test_iteration(self):
        a = Coords(1, 2, 3)
        self.assertEqual(list(a), [1, 2, 3])

    def test_to_tuple_and_numpy(self):
        a = Coords(1, 2, 3)
        self.assertEqual(a.to_tuple(), (1.0, 2.0, 3.0))
        np.testing.assert_array_equal(a.to_numpy(), np.array([1.0, 2.0, 3.0]))

    def test_from_numpy_valid(self):
        arr = [1, 2.5, -3]
        c = Coords.from_numpy(arr)
        self.assertEqual(c, Coords(1.0, 2.5, -3.0))

    def test_from_numpy_invalid_type(self):
        with self.assertRaises(TypeError):
            Coords.from_numpy("invalid")

    def test_from_numpy_wrong_length(self):
        with self.assertRaises(ValueError):
            Coords.from_numpy([1, 2])

    def test_from_numpy_invalid_values(self):
        with self.assertRaises(TypeError):
            Coords.from_numpy([1, "bad", 3])


class TestGetPerpendicularVector(unittest.TestCase):

    def test_3d_perpendicular_vector(self):
        v = np.array([1.0, 0.0, 0.0])
        perp = get_perpendicular_vector(v)
        np.testing.assert_allclose(np.dot(perp, v), 0.0, atol=1e-8)
        self.assertAlmostEqual(np.linalg.norm(perp), 1.0, places=6)

    def test_coords_input(self):
        v = Coords(0, 1, 0)
        perp = get_perpendicular_vector(v)
        np.testing.assert_allclose(np.dot(perp, v.to_numpy()), 0.0, atol=1e-8)
        self.assertAlmostEqual(np.linalg.norm(perp), 1.0, places=6)

    def test_vector_not_normalized(self):
        v = [1.0, 2.0, 3.0]
        perp = get_perpendicular_vector(v, normalize=False)
        self.assertAlmostEqual(np.dot(perp, v), 0.0, places=6)

    def test_zero_vector(self):
        with self.assertRaises(ValueError):
            get_perpendicular_vector([0, 0, 0])

    def test_2d_vector(self):
        v = [1, 0]
        perp = get_perpendicular_vector(v)
        self.assertEqual(len(perp), 2)
        self.assertAlmostEqual(np.dot(perp, v), 0.0, places=6)


if __name__ == "__main__":
    unittest.main()

import unittest 
import itertools
from bsplinemi import *
from pymibspline.bspline import *

class TestUtils(unittest.TestCase):
    def setUp(self):
        self.x = [1, 3, 4, 6, 3, 6, 7, 3, 2]
        self.y = [3, 2, 4, 4, 3, 2, 4, 4, 2]
        self.bs = 5
        self.so = 3

    def test_basis_function(self):
        z = x2z(self.x)
        knots = knot_vector(self.bs, self.so)
        for (b, n) in itertools.product(range(5), range(len(z))):
            self.assertTrue(abs(basis_function(b, self.so, z[n], knots, self.bs) - basisFunction(b, self.so, z[n], knots, self.bs)) < 1E-10)

    def test_find_weights(self):
        w = find_weights(self.x, 5, 3)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUtils)
    unittest.TextTestRunner(verbosity=2).run(suite)

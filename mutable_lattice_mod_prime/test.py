import unittest

from . import (
    is_prime8, is_prime32, inverse32
)

from .pylattice import LatticeModPrime as PyLatticeModPrime

class TestPyLatticeModPrime(unittest.TestCase):
    def test_construction(self):
        for p in [-10, -1, 0, 1, 4, 6, 8, 9, 10]:
            with self.assertRaisesRegex(ValueError, "is not prime"):
                PyLatticeModPrime(p, 10)
        for p in [2, 3, 5, 7, 11, 13]:
            PyLatticeModPrime(p, 10)

    def test_basic_1(self):
        L = PyLatticeModPrime(7, 3)
        self.assertEqual(L.basis, [])
        L.add_vector([1, 0, 0])
        self.assertEqual(L.basis, [[1, 0, 0]])
        L.add_vector([0, 2, 0])
        self.assertEqual(L.basis, [[1, 0, 0], [0, 1, 0]])
        L.add_vector([0, 0, 3])
        self.assertEqual(L.basis, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    def test_basic_2(self):
        L = PyLatticeModPrime(7, 3)
        L.add_vector([5, 5, 5])
        self.assertEqual(L.basis, [[1, 1, 1]])
        L.add_vector([0, 2, 3])
        self.assertEqual(L.basis, [[1, 1, 1], [0, 1, 5]])

SMALL_PRIME_LIST = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251
]


class Test_IsPrime(unittest.TestCase):
    def test_is_prime8(self):
        self.assertEqual([x for x in range(256) if is_prime8(x)], SMALL_PRIME_LIST)
    
    def test_is_prime8_overflow(self):
        with self.assertRaises(OverflowError):
            is_prime8(256)

    def test_is_prime32_small(self):
        self.assertTrue(is_prime32(23))
        self.assertTrue(is_prime32(61))
        self.assertEqual([x for x in range(256) if is_prime32(x)], SMALL_PRIME_LIST)

    def test_is_prime32_large(self):
        self.assertEqual([x for x in range(2**32-1000, 2**32) if is_prime32(x)],
                         [4294966297, 4294966337, 4294966367, 4294966373, 4294966427, 4294966441,
                          4294966447, 4294966477, 4294966553, 4294966583, 4294966591, 4294966619,
                          4294966639, 4294966651, 4294966657, 4294966661, 4294966667, 4294966769,
                          4294966813, 4294966829, 4294966877, 4294966909, 4294966927, 4294966943,
                          4294966981, 4294966997, 4294967029, 4294967087, 4294967111, 4294967143,
                          4294967161, 4294967189, 4294967197, 4294967231, 4294967279, 4294967291])

class Test_ModInverse(unittest.TestCase):
    def test_inverse32(self):
        for p in SMALL_PRIME_LIST:
            for a in range(1, p):
                b = inverse32(a, p)
                self.assertIn(b, range(p))
                self.assertEqual((a*b) % p, 1)
                self.assertEqual(b, pow(a, -1, p))
    
    def test_inverse32_random(self):
        import random
        random.seed(123)
        for _ in range(1000):
            a, p = sorted([random.getrandbits(32), random.getrandbits(32)])
            if not is_prime32(p):
                continue
            b = inverse32(a, p)
            self.assertIn(b, range(p))
            self.assertEqual((a*b) % p, 1)

if __name__ == "__main__":
    unittest.main()
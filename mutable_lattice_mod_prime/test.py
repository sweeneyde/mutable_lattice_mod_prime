import unittest

from . import (
    is_prime_u8, is_prime_u32, inverse_u32, SMALL_PRIMES,
    Prime_u8, Prime_u32,
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

class Test_IsPrime(unittest.TestCase):
    def test_is_prime_u8(self):
        self.assertEqual([x for x in range(256) if is_prime_u8(x)],
                         [x for x in SMALL_PRIMES if x < 256])

    def test_is_prime_u8_overflow(self):
        with self.assertRaises(OverflowError):
            is_prime_u8(256)

    def test_is_prime_u32_small(self):
        self.assertTrue(is_prime_u32(23))
        self.assertTrue(is_prime_u32(61))
        self.assertEqual([x for x in range(max(SMALL_PRIMES) + 1) if is_prime_u32(x)],
                         [x for x in SMALL_PRIMES])

    def test_is_prime_u32_large(self):
        self.assertEqual([x for x in range(2**32-1000, 2**32) if is_prime_u32(x)],
                         [4294966297, 4294966337, 4294966367, 4294966373, 4294966427, 4294966441,
                          4294966447, 4294966477, 4294966553, 4294966583, 4294966591, 4294966619,
                          4294966639, 4294966651, 4294966657, 4294966661, 4294966667, 4294966769,
                          4294966813, 4294966829, 4294966877, 4294966909, 4294966927, 4294966943,
                          4294966981, 4294966997, 4294967029, 4294967087, 4294967111, 4294967143,
                          4294967161, 4294967189, 4294967197, 4294967231, 4294967279, 4294967291])

class Test_ModInverse(unittest.TestCase):
    def test_inverse_u32(self):
        for p in SMALL_PRIMES:
            for a in range(1, p):
                b = inverse_u32(a, p)
                self.assertIn(b, range(p))
                self.assertEqual((a*b) % p, 1)
                self.assertEqual(b, pow(a, -1, p))
    
    def test_inverse_u32_random(self):
        import random
        random.seed(123)
        for _ in range(1000):
            a, p = sorted([random.getrandbits(32), random.getrandbits(32)])
            if not is_prime_u32(p):
                continue
            b = inverse_u32(a, p)
            self.assertIn(b, range(p))
            self.assertEqual((a*b) % p, 1)

class Test_Prime_u8(unittest.TestCase):
    def test_prime_u8_creation(self):
        Prime_u8(2)
        Prime_u8(3)
        Prime_u8(251)

    def test_prime_u8_round_trip(self):
        for p in [2, 3, 251]:
            self.assertEqual(Prime_u8(p).p, p)
            self.assertEqual(repr(Prime_u8(p)), f"Prime_u8({p})")
            self.assertEqual(str(Prime_u8(p)), f"Prime_u8({p})")

    def test_prime_u8_error(self):
        with self.assertRaises(OverflowError): Prime_u8(-1)
        with self.assertRaises(OverflowError): Prime_u8(256)
        with self.assertRaises(OverflowError): Prime_u8(500)
        with self.assertRaises(OverflowError): Prime_u8(1000)
        with self.assertRaisesRegex(ValueError, "not prime"): Prime_u8(0)
        with self.assertRaisesRegex(ValueError, "not prime"): Prime_u8(1)
        with self.assertRaisesRegex(ValueError, "not prime"): Prime_u8(4)
        with self.assertRaisesRegex(ValueError, "not prime"): Prime_u8(255)

    def test_inverses(self):
        self.assertEqual(list(map(Prime_u8(2).inverse, range(2))), [None, 1])
        self.assertEqual(list(map(Prime_u8(3).inverse, range(3))), [None, 1, 2])
        self.assertEqual(list(map(Prime_u8(5).inverse, range(5))), [None, 1, 3, 2, 4])
        self.assertEqual(list(map(Prime_u8(7).inverse, range(7))), [None, 1, 4, 5, 2, 3, 6])
        self.assertEqual(list(map(Prime_u8(11).inverse, range(11))), [None, 1, 6, 4, 3, 9, 2, 8, 7, 5, 10])
        self.assertEqual(list(map(Prime_u8(251).inverse, range(251))),
                         [None] + [pow(x, -1, 251) for x in range(1, 251)])

class Test_Prime_u32(unittest.TestCase):
    def test_prime_u32_creation(self):
        Prime_u32(2)
        Prime_u32(3)
        Prime_u32(251)
        Prime_u32(7919)
        Prime_u32(4294967291)

    def test_prime_u32_round_trip(self):
        for p in [2, 3, 251, 7919, 4294967291]:
            self.assertEqual(Prime_u32(p).p, p)
            self.assertEqual(repr(Prime_u32(p)), f"Prime_u32({p})")
            self.assertEqual(str(Prime_u32(p)), f"Prime_u32({p})")

    def test_prime_u32_error(self):
        with self.assertRaises(OverflowError): Prime_u32(-1)
        with self.assertRaises(OverflowError): Prime_u32(2**32)
        with self.assertRaises(OverflowError): Prime_u32(2**64)
        with self.assertRaises(OverflowError): Prime_u32(2**100)
        with self.assertRaisesRegex(ValueError, "not prime"): Prime_u32(0)
        with self.assertRaisesRegex(ValueError, "not prime"): Prime_u32(1)
        with self.assertRaisesRegex(ValueError, "not prime"): Prime_u32(4)
        with self.assertRaisesRegex(ValueError, "not prime"): Prime_u32(255)
        with self.assertRaisesRegex(ValueError, "not prime"): Prime_u32(2**32 - 3)

    def test_inverses(self):
        self.assertEqual(list(map(Prime_u32(2).inverse, range(2))), [None, 1])
        self.assertEqual(list(map(Prime_u32(3).inverse, range(3))), [None, 1, 2])
        self.assertEqual(list(map(Prime_u32(5).inverse, range(5))), [None, 1, 3, 2, 4])
        self.assertEqual(list(map(Prime_u32(7).inverse, range(7))), [None, 1, 4, 5, 2, 3, 6])
        self.assertEqual(list(map(Prime_u32(11).inverse, range(11))), [None, 1, 6, 4, 3, 9, 2, 8, 7, 5, 10])
        self.assertEqual(list(map(Prime_u32(251).inverse, range(251))),
                         [None] + [pow(x, -1, 251) for x in range(1, 251)])
        p = 4294967291
        self.assertEqual(Prime_u32(p).inverse(1234), pow(1234, -1, p))


if __name__ == "__main__":
    unittest.main()

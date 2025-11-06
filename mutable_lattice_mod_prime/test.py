import unittest
from array import array
import sys

WORD = {2**31-1: 'i', 2**63-1: 'q'}[sys.maxsize]
assert array(WORD, []).itemsize == {'i': 4, 'q': 8}[WORD]

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

class Test_VectorModPrime_u8(unittest.TestCase):
    def test_VectorModPrime_u8_creation(self):
        Prime_u8(2).Vector([0, 1, 1])

    def test_VectorModPrime_u8_round_trip(self):
        v = Prime_u8(5).Vector([4, 3, 2, 1, 0])
        self.assertEqual(v.tolist(), [4, 3, 2, 1, 0])
        self.assertEqual(v.p, 5)

    def test_VectorModPrime_u8_add(self):
        pp = Prime_u8(7)
        self.assertEqual(pp.Vector([0, 1, 2, 3, 4, 5, 6]) + pp.Vector([2, 4, 6, 1, 3, 5, 0]),
                         pp.Vector([2, 5, 1, 4, 0, 3, 6]))
        self.assertEqual(pp.Vector([]), pp.Vector([]), pp.Vector([]))
        pp = Prime_u8(251)
        self.assertEqual(pp.Vector([100, 200]) + pp.Vector([100, 200]), pp.Vector([200, 149]))

    def test_VectorModPrime_u8_subtract(self):
        pp = Prime_u8(7)
        self.assertEqual(pp.Vector([0, 1, 2, 3, 4, 5, 6]) - pp.Vector([2, 4, 6, 1, 3, 5, 0]),
                         pp.Vector([5, 4, 3, 2, 1, 0, 6]))
        self.assertEqual(pp.Vector([]), pp.Vector([]), pp.Vector([]))
        pp = Prime_u8(251)
        self.assertEqual(pp.Vector([1, 2, 3]) - pp.Vector([0, 100, 200]), pp.Vector([1, 153, 54]))

    def test_VectorModPrime_u8_inplace_add(self):
        pp = Prime_u8(7)
        v = pp.Vector([0, 1, 2, 3, 4, 5, 6])
        v0 = v
        v += pp.Vector([4, 4, 4, 4, 4, 4, 4])
        self.assertEqual(v, pp.Vector([4, 5, 6, 0, 1, 2, 3]))
        self.assertIs(v, v0)

    def test_VectorModPrime_u8_inplace_subtract(self):
        pp = Prime_u8(7)
        v = pp.Vector([4, 4, 4, 4, 4, 4, 4])
        v0 = v
        v -= pp.Vector([0, 1, 2, 3, 4, 5, 6])
        self.assertEqual(v, pp.Vector([4, 3, 2, 1, 0, 6, 5]))
        self.assertIs(v, v0)

    def test_VectorModPrime_u8_mul(self):
        pp = Prime_u8(251)
        v = pp.Vector([0, 20, 40, 60, 80, 100])
        self.assertEqual(v * 3, pp.Vector([0, 60, 120, 180, 240, 49]))
        self.assertEqual(3 * v, pp.Vector([0, 60, 120, 180, 240, 49]))
        self.assertEqual(v * (251*10**100 + 3), pp.Vector([0, 60, 120, 180, 240, 49]))

    def test_VectorModPrime_u8_imul(self):
        pp = Prime_u8(251)
        v = pp.Vector([0, 20, 40, 60, 80, 100])
        v0 = v
        v *= 3
        self.assertEqual(v, pp.Vector([0, 60, 120, 180, 240, 49]))
        self.assertIs(v, v0)

    def test_VectorModPrime_u8_negative(self):
        pp = Prime_u8(251)
        self.assertEqual(-pp.Vector([0, 1, 2, 3, 248, 249, 250]), pp.Vector([0, 250, 249, 248, 3, 2, 1]))

    def test_VectorModPrime_u8_subtract_multiple(self):
        pp = Prime_u8(251)
        v = pp.Vector([0, 10, 20, 30])
        v.subtract_multiple(5, pp.Vector([1, 2, 3, 4]))
        self.assertEqual(v, pp.Vector([246, 0, 5, 10]))

    def test_VectorModPrime_u8_repr(self):
        self.assertEqual(
            repr(Prime_u8(251).Vector([0, 1, 2])),
            "Prime_u8(251).Vector([0, 1, 2])"
        )

    def test_VectorModPrime_u8_str(self):
        self.assertEqual(str(Prime_u8(251).Vector([])), "[]")
        self.assertEqual(
            str(Prime_u8(7).Vector([1, 2, 3, 4, 5, 6])),
            "[1 2 3 4 5 6]"
        )
        self.assertEqual(
            str(Prime_u8(19).Vector([1, 2, 3, 10, 11, 12])),
            "[ 1  2  3 10 11 12]"
        )
        self.assertEqual(
            str(Prime_u8(251).Vector([123, 45, 6, 78, 9])),
            "[123  45   6  78   9]"
        )

    def test_shuffled_by_action(self):
        pp = Prime_u8(97)
        v = pp.Vector([0, 10, 20, 30, 40])
        self.assertEqual(v.shuffled_by_action(array(WORD, [4, 3, 2, 1, 0])),
                         pp.Vector([40, 30, 20, 10, 0]))
        self.assertEqual(v.shuffled_by_action(array(WORD, [4, 3, 2, 1, 0]), 5),
                         pp.Vector([40, 30, 20, 10, 0]))
        self.assertEqual(v.shuffled_by_action(array(WORD, [0, 1, 0, 1, 0])),
                         pp.Vector([60, 40, 0, 0, 0]))
        self.assertEqual(v.shuffled_by_action(array(WORD, [0, 1, 0, 1, 0]), 2),
                         pp.Vector([60, 40]))
        self.assertEqual(v.shuffled_by_action(array(WORD, [2, 2, 2, 2, 2])),
                         pp.Vector([0, 0, 100 % 97, 0, 0]))

class Test_Lattice_u8(unittest.TestCase):
    def test_creation(self):
        pp = Prime_u8(7)
        pp.Lattice(3)
        pp.Lattice(3, [])
        pp.Lattice(3, [[2, 1, 1]])

    def test_tolist(self):
        pp = Prime_u8(7)
        self.assertEqual(pp.Lattice(3).tolist(), [])
        self.assertEqual(pp.Lattice(3, []).tolist(), [])
        self.assertEqual(pp.Lattice(3, [[2, 1, 1]]).tolist(), [[1, 4, 4]])

    def test_get_basis(self):
        pp = Prime_u8(7)
        self.assertEqual(pp.Lattice(3).get_basis(), [])
        self.assertEqual(pp.Lattice(3, []).get_basis(), [])
        self.assertEqual(pp.Lattice(3, [[2, 1, 1]]).get_basis(), [pp.Vector([1, 4, 4])])

    def test_row_echelon_1(self):
        pp = Prime_u8(7)
        self.assertEqual(pp.Lattice(3, [[2, 2, 2], [3, 3, 3], [4, 4, 4]]).tolist(), [[1, 1, 1]])
        self.assertEqual(pp.Lattice(3, [[2, 2, 2], [0, 3, 3], [0, 0, 4]]).tolist(), [[1, 1, 1], [0, 1, 1], [0, 0, 1]])
        self.assertEqual(pp.Lattice(3, [[1, 0, 1], [2, 0, 3]]).tolist(), [[1, 0, 1], [0, 0, 1]])
        self.assertEqual(pp.Lattice(3, [[2, 0, 3], [1, 0, 1]]).tolist(), [[1, 0, 5], [0, 0, 1]])

    def test_row_echelon_2(self):
        pp = Prime_u8(3)
        self.assertEqual(pp.Lattice(0, [[], [], []]).tolist(), [])
        self.assertEqual(pp.Lattice(1, []).tolist(), [])
        self.assertEqual(pp.Lattice(1, [[0]]).tolist(), [])
        self.assertEqual(pp.Lattice(1, [[1]]).tolist(), [[1]])
        self.assertEqual(pp.Lattice(1, [[2]]).tolist(), [[1]])
        self.assertEqual(pp.Lattice(1, [[2], [1]]).tolist(), [[1]])
        self.assertEqual(pp.Lattice(2, []).tolist(), [])
        self.assertEqual(pp.Lattice(2, [[0, 0]]).tolist(), [])
        self.assertEqual(pp.Lattice(2, [[1, 0]]).tolist(), [[1, 0]])
        self.assertEqual(pp.Lattice(2, [[2, 0]]).tolist(), [[1, 0]])
        self.assertEqual(pp.Lattice(2, [[1, 1]]).tolist(), [[1, 1]])
        self.assertEqual(pp.Lattice(2, [[2, 2]]).tolist(), [[1, 1]])
        self.assertEqual(pp.Lattice(2, [[2, 1]]).tolist(), [[1, 2]])
        self.assertEqual(pp.Lattice(2, [[0, 0], [2, 0]]).tolist(), [[1, 0]])
        self.assertEqual(pp.Lattice(2, [[1, 2], [2, 1]]).tolist(), [[1, 2]])
        self.assertEqual(pp.Lattice(2, [[1, 1], [2, 0]]).tolist(), [[1, 1], [0, 1]])
        self.assertEqual(pp.Lattice(2, [[0, 2], [2, 2]]).tolist(), [[1, 1], [0, 1]])

    def test_contains(self):
        pp = Prime_u8(3)
        L = pp.Lattice(6, [[1, 1, 0, 0, 2, 2]])
        self.assertIn(pp.Vector([2, 2, 0, 0, 1, 1]), L)
        self.assertIn(pp.Vector([1, 1, 0, 0, 2, 2]), L)
        self.assertIn(pp.Vector([0, 0, 0, 0, 0, 0]), L)
        self.assertNotIn(pp.Vector([0, 0, 1, 0, 0, 0]), L)
        self.assertNotIn(pp.Vector([1, 1, 0, 0, 1, 1]), L)
        self.assertNotIn(pp.Vector([0, 0, 1, 0, 0, 0]), L)
        self.assertNotIn(pp.Vector([1, 1, 0, 0, 0, 0]), L)
        L = pp.Lattice(6, [[1, 1, 0, 0, 2, 2], [0, 0, 1, 1, 1, 1]])
        self.assertIn(pp.Vector([0, 0, 0, 0, 0, 0]), L)
        self.assertIn(pp.Vector([1, 1, 0, 0, 2, 2]), L)
        self.assertIn(pp.Vector([2, 2, 0, 0, 1, 1]), L)
        self.assertIn(pp.Vector([0, 0, 1, 1, 1, 1]), L)
        self.assertIn(pp.Vector([1, 1, 1, 1, 0, 0]), L)
        self.assertIn(pp.Vector([2, 2, 1, 1, 2, 2]), L)
        self.assertIn(pp.Vector([0, 0, 2, 2, 2, 2]), L)
        self.assertIn(pp.Vector([1, 1, 2, 2, 1, 1]), L)
        self.assertIn(pp.Vector([2, 2, 2, 2, 0, 0]), L)
        self.assertNotIn(pp.Vector([1, 1, 1, 1, 1, 1]), L)
        self.assertNotIn(pp.Vector([2, 2, 1, 1, 0, 0]), L)
        self.assertNotIn(pp.Vector([0, 0, 1, 1, 2, 2]), L)
        self.assertNotIn(pp.Vector([1, 1, 2, 2, 2, 2]), L)
        self.assertNotIn(pp.Vector([2, 2, 2, 2, 1, 1]), L)
        self.assertNotIn(pp.Vector([0, 0, 2, 2, 0, 0]), L)
        self.assertNotIn(pp.Vector([1, 1, 0, 0, 0, 0]), L)
        self.assertNotIn(pp.Vector([2, 2, 0, 0, 2, 2]), L)
        self.assertNotIn(pp.Vector([0, 0, 0, 0, 1, 1]), L)

if __name__ == "__main__":
    unittest.main()

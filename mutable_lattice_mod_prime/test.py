from array import array
import unittest
import itertools
import sys
import random

WORD = {2**31-1: 'i', 2**63-1: 'q'}[sys.maxsize]
assert array(WORD, []).itemsize == {'i': 4, 'q': 8}[WORD]

from . import (
    is_prime_u8, is_prime_u32, inverse_u32, PRIMES_UNDER_256, PRIMES_UNDER_1000,
    Prime_u8, Prime_u32, Prime2,
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
        self.assertEqual([x for x in range(256) if is_prime_u8(x)], list(PRIMES_UNDER_256))

    def test_is_prime_u8_overflow(self):
        with self.assertRaises(OverflowError):
            is_prime_u8(256)

    def test_is_prime_u32_small(self):
        self.assertTrue(is_prime_u32(23))
        self.assertTrue(is_prime_u32(61))
        self.assertEqual([x for x in range(1000) if is_prime_u32(x)], list(PRIMES_UNDER_1000))

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
        for p in PRIMES_UNDER_1000:
            for a in range(1, p):
                b = inverse_u32(a, p)
                self.assertIn(b, range(p))
                self.assertEqual((a*b) % p, 1)
                self.assertEqual(b, pow(a, -1, p))
    
    def test_inverse_u32_random(self):
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

    def test_length(self):
        pp = Prime_u8(71)
        self.assertEqual(len(pp.Vector([])), 0)
        self.assertEqual(len(pp.Vector([10])), 1)
        self.assertEqual(len(pp.Vector([10, 20])), 2)
        self.assertEqual(len(pp.Vector([10, 20, 30])), 3)

    def test_getitem(self):
        v = Prime_u8(71).Vector([0, 10, 20, 30, 40])
        self.assertEqual(v[0], 0)
        self.assertEqual(v[1], 10)
        self.assertEqual(v[2], 20)
        self.assertEqual(v[3], 30)
        self.assertEqual(v[4], 40)
        with self.assertRaises(IndexError): v[5]
        self.assertEqual(v[-1], 40)
        self.assertEqual(v[-2], 30)
        self.assertEqual(v[-3], 20)
        self.assertEqual(v[-4], 10)
        self.assertEqual(v[-5], 0)
        with self.assertRaises(IndexError): v[-6]

    def test_iterable(self):
        self.assertEqual(tuple(Prime_u8(5).Vector([1, 2, 3])), (1, 2, 3))

    def test_setitem(self):
        v = Prime_u8(71).Vector([0, 10, 20, 30, 40])
        self.assertEqual(v[2], 20)
        v[2] = 21
        self.assertEqual(v[2], 21)
        self.assertEqual(v.tolist(), [0, 10, 21, 30, 40])

    def test_copy(self):
        v = Prime_u8(71).Vector([0, 10, 20, 30, 40])
        w = v.copy()
        v[1], v[2], v[3] = 11, 21, 31
        w[1], w[2], w[3] = 12, 22, 32
        self.assertEqual(v.tolist(), [0, 11, 21, 31, 40])
        self.assertEqual(w.tolist(), [0, 12, 22, 32, 40])

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

    def test_VectorModPrime_u8_add_multiple(self):
        pp = Prime_u8(251)
        v = pp.Vector([0, 10, 20, 30])
        v.add_multiple(5, pp.Vector([1, 2, 3, 4]))
        self.assertEqual(v, pp.Vector([5, 20, 35, 50]))

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

    def test_vector_ops_exhaustive(self):
        for p in (2, 3, 5, 7, 11, 13):
            pp = Prime_u8(p)
            data1 = [i for i in range(p) for j in range(p)]
            data2 = [j for i in range(p) for j in range(p)]
            self.assertEqual((pp.Vector(data1) + pp.Vector(data2)).tolist(),
                             [(a + b) % p for a, b in zip(data1, data2)])
            self.assertEqual((pp.Vector(data1) - pp.Vector(data2)).tolist(),
                             [(a - b) % p for a, b in zip(data1, data2)])

class Test_Lattice_u8(unittest.TestCase):
    def test_creation(self):
        pp = Prime_u8(7)
        pp.Lattice(3)
        pp.Lattice(3, [])
        pp.Lattice(3, [[2, 1, 1]])
        pp.Lattice(3, maxrank=1)
        pp.Lattice(3, [[0, 1, 2], [1, 0, 0]], maxrank=-1)

    def test_tolist(self):
        pp = Prime_u8(7)
        self.assertEqual(pp.Lattice(3).tolist(), [])
        self.assertEqual(pp.Lattice(3, []).tolist(), [])
        self.assertEqual(pp.Lattice(3, [[2, 1, 1]]).tolist(), [[1, 4, 4]])

    def test_str(self):
        p7 = Prime_u8(7)
        self.assertEqual(str(p7.Lattice(0)), "<zero subspace of Z/7Z ^ 0>")
        self.assertEqual(str(p7.Lattice(1)), "<zero subspace of Z/7Z ^ 1>")
        self.assertEqual(str(p7.Lattice(2)), "<zero subspace of Z/7Z ^ 2>")
        self.assertEqual(str(p7.Lattice(1, [[1]])), "[1]\n")
        self.assertEqual(str(p7.Lattice(2, [[1, 2]])), "[1 2]\n")
        self.assertEqual(str(p7.Lattice(2, [[1, 2], [0, 1]])), "[1 2]\n[0 1]\n")
        p13 = Prime_u8(13)
        self.assertEqual(str(p13.Lattice(0)), "<zero subspace of Z/13Z ^ 0>")
        self.assertEqual(str(p13.Lattice(1)), "<zero subspace of Z/13Z ^ 1>")
        self.assertEqual(str(p13.Lattice(2)), "<zero subspace of Z/13Z ^ 2>")
        self.assertEqual(str(p13.Lattice(1, [[1]])), "[ 1]\n")
        self.assertEqual(str(p13.Lattice(2, [[1, 2]])), "[ 1  2]\n")
        self.assertEqual(str(p13.Lattice(3, [[1, 2, 4], [0, 1, 10]])), "[ 1  2  4]\n[ 0  1 10]\n")
        p101 = Prime_u8(101)
        self.assertEqual(str(p101.Lattice(0)), "<zero subspace of Z/101Z ^ 0>")
        self.assertEqual(str(p101.Lattice(1)), "<zero subspace of Z/101Z ^ 1>")
        self.assertEqual(str(p101.Lattice(2)), "<zero subspace of Z/101Z ^ 2>")
        self.assertEqual(str(p101.Lattice(1, [[1]])), "[  1]\n")
        self.assertEqual(str(p101.Lattice(2, [[1, 2]])), "[  1   2]\n")
        self.assertEqual(str(p101.Lattice(4, [[1, 2, 4, 8], [0, 1, 10, 100]])), "[  1   2   4   8]\n"
                                                                                "[  0   1  10 100]\n")

    def test_repr(self):
        self.assertEqual(repr(Prime_u8(7).Lattice(0, [])), "Prime_u8(7).Lattice(0, [])")
        self.assertEqual(repr(Prime_u8(7).Lattice(1, [])), "Prime_u8(7).Lattice(1, [])")
        self.assertEqual(repr(Prime_u8(7).Lattice(1, [[1]])), "Prime_u8(7).Lattice(1, [[1]])")
        self.assertEqual(repr(Prime_u8(7).Lattice(2, [])), "Prime_u8(7).Lattice(2, [])")
        self.assertEqual(repr(Prime_u8(7).Lattice(2, [[1, 1]])), "Prime_u8(7).Lattice(2, [[1, 1]])")
        self.assertEqual(repr(Prime_u8(7).Lattice(2, [[1, 1], [0, 1]])), "Prime_u8(7).Lattice(2, [[1, 1], [0, 1]])")

    def test_get_basis(self):
        pp = Prime_u8(7)
        self.assertEqual(pp.Lattice(3).get_basis(), [])
        self.assertEqual(pp.Lattice(3, []).get_basis(), [])
        self.assertEqual(pp.Lattice(3, [[2, 1, 1]]).get_basis(), [pp.Vector([1, 4, 4])])

    def test_full(self):
        pp = Prime_u8(3)
        self.assertEqual(pp.Lattice_full(0).tolist(), [])
        self.assertEqual(pp.Lattice_full(1).tolist(), [[1]])
        self.assertEqual(pp.Lattice_full(2).tolist(), [[1, 0], [0, 1]])
        self.assertEqual(pp.Lattice_full(3).tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    def test_maxrank(self):
        pp = Prime_u8(7)
        L = pp.Lattice(8, [[1]*8], maxrank=2)
        self.assertEqual(L.maxrank, 2)
        self.assertEqual(L.rank, 1)
        L.add_vector(pp.Vector([0, 0, 0, 0, 1, 2, 3, 4]))
        self.assertEqual(L.rank, 2)
        L.add_vector(pp.Vector([-1, -1, -1, -1, 0, 1, 2, 3]))
        self.assertEqual(L.rank, 2)
        with self.assertRaises(IndexError):
            L.add_vector(pp.Vector([0, 1, 0, 0, 0, 0, 0, 0]))

    def test_row_echelon_1(self):
        pp = Prime_u8(7)
        self.assertEqual(pp.Lattice(3, [[2, 2, 2], [3, 3, 3], [4, 4, 4]]).tolist(), [[1, 1, 1]])
        self.assertEqual(pp.Lattice(3, [[2, 2, 2], [0, 3, 3], [0, 0, 4]]).tolist(), [[1, 1, 1], [0, 1, 1], [0, 0, 1]])
        self.assertEqual(pp.Lattice(3, [[1, 0, 1], [2, 0, 3]]).tolist(), [[1, 0, 1], [0, 0, 1]])
        self.assertEqual(pp.Lattice(3, [[2, 0, 3], [1, 0, 1]]).tolist(), [[1, 0, 5], [0, 0, 1]])
        self.assertEqual(pp.Lattice(4, [[2, 1, 1, 3], [0, 0, 4, 5]]).tolist(), [[1, 4, 4, 5], [0, 0, 1, 3]])

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

    def test_equality(self):
        pp = Prime_u8(5)
        self.assertEqual(pp.Lattice(2, [[2, 2]]), pp.Lattice(2, [[3, 3]]))
        self.assertEqual(pp.Lattice(2, [[1, 3], [0, 1]]), pp.Lattice(2, [[1, 0], [0, 1]]))
        self.assertEqual(pp.Lattice(4, [[1, 2, 3, 4], [0, 0, 1, 3]]), pp.Lattice(4, [[1, 2, 0, 0], [0, 0, 1, 3]]))
        self.assertNotEqual(pp.Lattice(2), pp.Lattice(3))
        self.assertNotEqual(pp.Lattice(3, []), pp.Lattice(3, [[1, 1, 1]]))
        self.assertNotEqual(pp.Lattice(3, [[1, 1, 1]]), pp.Lattice(3, [[1, 2, 1]]))

    def test_comparison_exhaustive(self):
        p = 3
        pp = Prime_u8(p)
        # get all row echelon forms
        vectors = [[1,0,0],[1,0,1],[1,0,2],
                   [1,1,0],[1,1,1],[1,1,2],
                   [1,2,0],[1,2,1],[1,2,2],
                   [0,1,0],[0,1,1],[0,1,2],[0,0,1]]
        pairs = [[u, v] for u, v in itertools.combinations(vectors, 2) if u.index(1) < v.index(1)]
        triples = [[u, v, w] for u, v in pairs for w in vectors if v.index(1) < w.index(1)]
        spaces = (
            [pp.Lattice(3)]
            + [pp.Lattice(3, [v]) for v in vectors]
            + [pp.Lattice(3, pair) for pair in pairs]
            + [pp.Lattice(3, triple) for triple in triples]
        )
        for a in spaces:
            for b in spaces:
                a_in_b = all(map(b.__contains__, a.get_basis()))
                b_in_a = all(map(a.__contains__, b.get_basis()))
                self.assertEqual(a == b, a_in_b and b_in_a)
                self.assertEqual(a != b, not (a_in_b and b_in_a))
                self.assertEqual(a <= b, a_in_b)
                self.assertEqual(a < b, a_in_b and not b_in_a)
                self.assertEqual(a >= b, b_in_a)
                self.assertEqual(a > b, b_in_a and not a_in_b)

    def test_add_vector(self):
        pp = Prime_u8(11)
        L = pp.Lattice(4)
        self.assertEqual(L.tolist(), [])
        L.add_vector(pp.Vector([0, 0, 0, 4]))
        self.assertEqual(L.tolist(), [[0, 0, 0, 1]])
        L.add_vector(pp.Vector([0, 2, 2, 2]))
        self.assertEqual(L.tolist(), [[0, 1, 1, 1], [0, 0, 0, 1]])
        L.add_vector(pp.Vector([1, 1, 1, 1]))
        self.assertEqual(L.tolist(), [[1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 0, 1]])
        L.add_vector(pp.Vector([0, 0, 3, 3]))
        self.assertEqual(L.tolist(), [[1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1]])
        L.add_vector(pp.Vector([1, 2, 3, 5]))
        self.assertEqual(L.tolist(), [[1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1]])

    def test_copy(self):
        pp = Prime_u8(7)
        L = pp.Lattice(3, [[1, 2, 3], [0, 1, 2]])
        L2 = L.copy()
        L2.add_vector(pp.Vector([0, 0, 1]))
        self.assertEqual(L.tolist(), [[1, 2, 3], [0, 1, 2]])
        self.assertEqual(L2.tolist(), [[1, 2, 3], [0, 1, 2], [0, 0, 1]])

    def test_clear(self):
        pp = Prime_u8(7)
        for N in [0, 1, 2, 3, 4, 5]:
            L = pp.Lattice_full(N)
            basis = L.get_basis()
            self.assertIsNone(L.clear())
            self.assertEqual(L.get_basis(), [])
            for v in basis:
                L.add_vector(v) # make sure we can still add vectors

    def test_reduce(self):
        # reduced row echelon form
        pp = Prime_u8(7)
        L = pp.Lattice(3, [[2, 2, 2], [3, 3, 3], [4, 4, 4]])
        L.reduce()
        self.assertEqual(L.tolist(), [[1, 1, 1]])
        L = pp.Lattice(3, [[2, 2, 2], [0, 3, 3], [0, 0, 4]])
        L.reduce()
        self.assertEqual(L.tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        L = pp.Lattice(3, [[1, 0, 1], [2, 0, 3]])
        L.reduce()
        self.assertEqual(L.tolist(), [[1, 0, 0], [0, 0, 1]])
        L = pp.Lattice(3, [[2, 0, 3], [1, 0, 1]])
        L.reduce()
        self.assertEqual(L.tolist(), [[1, 0, 0], [0, 0, 1]])
        L = pp.Lattice(4, [[2, 1, 1, 3], [0, 0, 4, 5]])
        L.reduce()
        self.assertEqual(L.tolist(), [[1, 4, 0, 0], [0, 0, 1, 3]])

    def test_random(self):
        random.seed(123)
        p = 5
        pp = Prime_u8(p)
        for N in [0, 1, 2, 3, 4, 5, 10, 20]:
            for _ in range(200):
                mod = pp.Lattice(N)
                vecs = [
                    pp.Vector([random.randrange(p) for _ in range(N)])
                    for _ in range(random.randrange(N + 2))
                ]
                zero = pp.Vector([0] * N)
                vec0 = vecs[0] if vecs else zero
                for vec in vecs:
                    mod.add_vector(vec)
                    self.assertIn(vec, mod)
                    self.assertIn(zero-vec, mod)
                    self.assertIn(3*vec, mod)
                    self.assertIn(vec0 + vec, mod)
                    self.assertIn(3*vec0-4*vec, mod)

                if random.randrange(2):
                    copy = mod.copy()
                    mod.clear()
                    for v in copy.get_basis():
                        mod.add_vector(v)

                copy = mod.copy()
                copy.reduce()
                self.assertEqual(copy, mod)

                lin_combo = pp.Vector([0] * N)
                for vec in vecs:
                    lin_combo += random.randrange(p) * vec
                self.assertIn(lin_combo, mod)

                basis = mod.get_basis()
                for to_remove in basis:
                    submod = pp.Lattice(N)
                    for vec in basis:
                        if vec != to_remove:
                            submod.add_vector(vec)
                    self.assertNotIn(to_remove, submod)
                    also_outside = to_remove
                    for vec in submod.get_basis():
                        also_outside += random.randint(-2, 2) * vec
                    self.assertNotIn(also_outside, submod)

    def test_against_pylattice(self):
        random.seed(123)
        p = 2
        pp = Prime_u8(p)
        for N in [0, 1, 2, 3, 4, 5, 10, 20]:
            for _ in range(200):
                mod = pp.Lattice(N)
                vecs = [
                    pp.Vector([random.randrange(p) for _ in range(N)])
                    for _ in range(random.randrange(N + 2))
                ]
                pymod = PyLatticeModPrime(p, N)
                for vec in vecs:
                    mod.add_vector(vec)
                    pymod.add_vector(vec.tolist())
                    self.assertEqual(mod.tolist(), pymod.basis)

    def test_contains_does_not_mutate_vector(self):
        pp = Prime_u8(3)
        L = pp.Lattice(3, [[1, 0, 0], [0, 1, 1]])
        v = pp.Vector([2, 2, 2])
        self.assertEqual(v.tolist(), [2, 2, 2])
        v in L
        self.assertEqual(v.tolist(), [2, 2, 2])

    def test_decompose_simple(self):
        pp = Prime_u8(59)
        self.assertEqual(
            pp.Lattice(2, [[1, 0], [0, 1]]).decompose(),
            ([[0], [1]], [pp.Lattice(1, [[1]]), pp.Lattice(1, [[1]])])
        )
        self.assertEqual(
            pp.Lattice(2, [[1, 2], [0, 3]]).decompose(),
            ([[0], [1]], [pp.Lattice(1, [[1]]), pp.Lattice(1, [[1]])])
        )
        self.assertEqual(
            pp.Lattice(2, [[1, 2]]).decompose(),
            ([[0, 1]],
             [pp.Lattice(2, [[1, 2]])])
        )
        self.assertEqual(
            pp.Lattice(4, [[1, 0, 0, 2], [0, 1, 0, 0], [0, 0, 1, 1]]).decompose(),
            ([[0, 2, 3], [1]],
             [pp.Lattice(3, [[1, 0, 2], [0, 1, 1]]), pp.Lattice(1, [[1]])])
        )
        self.assertEqual(pp.Lattice(0).decompose(), ([], []))
        self.assertEqual(pp.Lattice(1).decompose(), ([[0]], [pp.Lattice(1)]))
        self.assertEqual(pp.Lattice(1, [[1]]).decompose(), ([[0]], [pp.Lattice(1, [[1]])]))
        self.assertEqual(
            pp.Lattice(10).decompose(),
            ([[i] for i in range(10)],
             [pp.Lattice(1, maxrank=0) for _ in range(10)])
        )
        self.assertEqual(
            pp.Lattice(5, [[0, 1, 0, 0, 0], [0, 0, 0, 1, 7]]).decompose(),
            ([[0], [1], [2], [3, 4]],
             [pp.Lattice(1), pp.Lattice(1, [[1]]), pp.Lattice(1), pp.Lattice(2, [[1, 7]])],
             )
        )
        self.assertEqual(
            pp.Lattice(5, [[0, 0, 1, 20, 30]]).decompose(),
            ([[0], [1], [2, 3, 4]],
             [pp.Lattice(1), pp.Lattice(1), pp.Lattice(3, [[1, 20, 30]])])
        )
        self.assertEqual(
            pp.Lattice(10, [[0,1,0,0,0,0, 0, 0,0,10],
                            [0,0,1,0,0,0, 0, 0,0, 0],
                            [0,0,0,1,0,0, 0, 0,0,30],
                            [0,0,0,0,1,0,40, 0,0, 0],
                            [0,0,0,0,0,1, 0, 0,0, 0],
                            [0,0,0,0,0,0, 0, 1,0, 0]]).decompose(),
            ([[0], [1, 3, 9], [2], [4, 6], [5], [7], [8]],
             [pp.Lattice(1),
              pp.Lattice(3, [[1, 0, 10], [0, 1, 30]]),
              pp.Lattice(1, [[1]]),
              pp.Lattice(2, [[1, 40]]),
              pp.Lattice(1, [[1]]),
              pp.Lattice(1, [[1]]),
              pp.Lattice(1),])
        )

    def test_decompose_exhaustive(self):
        pp = Prime_u8(71)
        vectors = [pp.Vector(list(v)) for v in itertools.product([0, 1, 2], repeat=4)]
        lattices = (
            [pp.Lattice(4, maxrank=0), pp.Lattice_full(4)]
            + [pp.Lattice(4, [v], maxrank=1) for v in vectors]
            + [pp.Lattice(4, [v, w], maxrank=2) for v, w in itertools.combinations(vectors, 2)]
        )
        for L in lattices:
            indexes, summands = L.decompose()
            self.assertEqual(len(indexes), len(summands))
            self.assertEqual(sum(map(len, indexes)), 4)
            self.assertEqual({x for arr in indexes for x in arr}, set(range(4)))
            self.assertEqual(sum(summand.rank for summand in summands), L.rank)
            self.assertEqual(sum(summand.ambient_dimension for summand in summands), L.ambient_dimension)
            for summand in summands:
                self.assertEqual(summand.maxrank, summand.rank)
                self.assertEqual(summand.decompose(), ([list(range(summand.ambient_dimension))], [summand]))
            reconstruction = pp.Lattice(4)
            for arr, summand in zip(indexes, summands):
                self.assertEqual(summand.ambient_dimension, len(arr))
                for vec0 in summand.get_basis():
                    vec = pp.Vector([0] * 4)
                    for i, x in enumerate(arr):
                        vec[x] = vec0[i]
                    reconstruction.add_vector(vec)
            self.assertEqual(L, reconstruction)

    def test_decompose_keeping_together(self):
        pp = Prime_u8(103)
        L = pp.Lattice(7, [[1,2,0,0,0,0,0],
                           [0,0,1,0,0,7,0],
                           [0,0,0,1,3,0,0]])
        self.assertEqual(L.decompose([]),
            ([[0, 1], [2, 5], [3, 4], [6]],
             [pp.Lattice(2, [[1, 2]]),
              pp.Lattice(2, [[1, 7]]),
              pp.Lattice(2, [[1, 3]]),
              pp.Lattice(1)])
        )
        self.assertEqual(L.decompose([[1, 5, 4, 6]]), ([[0,1,2,3,4,5,6]], [L]))
        self.assertEqual(L.decompose([[1, 5]]),
            ([[0, 1, 2, 5], [3, 4], [6]],
             [pp.Lattice(4, [[1, 2, 0, 0], [0, 0, 1, 7]]),
              pp.Lattice(2, [[1, 3]]),
              pp.Lattice(1)])
        )
        self.assertEqual(L.decompose([[0, 4]]),
            ([[0, 1, 3, 4], [2, 5], [6]],
             [pp.Lattice(4, [[1, 2, 0, 0], [0, 0, 1, 3]]),
              pp.Lattice(2, [[1, 7]]),
              pp.Lattice(1)])
        )
        self.assertEqual(L.decompose([[1, 6]]),
            ([[0, 1, 6], [2, 5], [3, 4]],
             [pp.Lattice(3, [[1, 2, 0]]),
              pp.Lattice(2, [[1, 7]]),
              pp.Lattice(2, [[1, 3]])])
        )
        self.assertEqual(L.decompose([[0, 6], [5, 4]]),
            ([[0, 1, 6], [2, 3, 4, 5]],
             [pp.Lattice(3, [[1, 2, 0]]),
              pp.Lattice(4, [[1, 0, 0, 7], [0, 1, 3, 0]])])
        )
        self.assertEqual(L.decompose([[0, 3], [2, 6]]),
            ([[0, 1, 3, 4], [2, 5, 6]],
             [pp.Lattice(4, [[1, 2, 0, 0], [0, 0, 1, 3]]),
              pp.Lattice(3, [[1, 7, 0]])])
        )
        self.assertEqual(L.decompose([[0, 2], [3, 6]]),
            ([[0, 1, 2, 5], [3, 4, 6]],
             [pp.Lattice(4, [[1, 2, 0, 0], [0, 0, 1, 7]]),
              pp.Lattice(3, [[1, 3, 0]])])
        )
        self.assertEqual(L.decompose([[1, 3, 6]]),
            ([[0, 1, 3, 4, 6], [2, 5]],
             [pp.Lattice(5, [[1, 2, 0, 0, 0], [0, 0, 1, 3, 0]]),
              pp.Lattice(2, [[1, 7]])])
        )

    def test_relations_among(self):
        pp = Prime_u8(7)
        self.assertEqual(pp.relations_among([]), pp.Lattice(0))
        self.assertEqual(pp.relations_among([pp.Vector([])]), pp.Lattice(1, [[1]]))
        self.assertEqual(pp.relations_among([pp.Vector([]), pp.Vector([])]), pp.Lattice(2, [[1, 0], [0, 1]]))
        self.assertEqual(pp.relations_among([pp.Vector([1])]), pp.Lattice(1))
        self.assertEqual(pp.relations_among([pp.Vector([0])]), pp.Lattice(1, [[1]]))
        self.assertEqual(pp.relations_among([pp.Vector([0]), pp.Vector([0])]), pp.Lattice(2, [[1, 0], [0, 1]]))
        self.assertEqual(pp.relations_among([pp.Vector([1]), pp.Vector([0])]), pp.Lattice(2, [[0, 1]]))
        self.assertEqual(pp.relations_among([pp.Vector([1]), pp.Vector([1])]), pp.Lattice(2, [[1, 6]]))
        self.assertEqual(pp.relations_among([pp.Vector([1]), pp.Vector([2])]), pp.Lattice(2, [[1, 3]]))
        self.assertEqual(pp.relations_among([pp.Vector([3]), pp.Vector([6])]), pp.Lattice(2, [[1, 3]]))
        self.assertEqual(pp.relations_among([pp.Vector([3]), pp.Vector([6]), pp.Vector([1])]), pp.Lattice(3, [[1, 3, 0], [0, 1, 1]]))
        self.assertEqual(pp.relations_among([pp.Vector([2, 3])]), pp.Lattice(1))
        self.assertEqual(pp.relations_among([pp.Vector([0, 0])]), pp.Lattice(1, [[1]]))
        self.assertEqual(pp.relations_among([pp.Vector([2, 3]), pp.Vector([1, 0])]), pp.Lattice(2))
        self.assertEqual(pp.relations_among([pp.Vector([2, 3]), pp.Vector([4, 6])]), pp.Lattice(2, [[1, 3]]))
        self.assertEqual(pp.relations_among([pp.Vector([3, 5]), pp.Vector([0, 0])]), pp.Lattice(2, [[0, 1]]))
        self.assertEqual(pp.relations_among([pp.Vector([3, 6]), pp.Vector([6, 5]), pp.Vector([1, 2])]), pp.Lattice(3, [[1, 3, 0], [0, 1, 1]]))

    def test_rank_nullity(self):
        random.seed(0)
        for p in (2, 3, 5):
            pp = Prime_u8(p)
            for N in range(5):
                zero = pp.Vector([0] * N)
                for R in range(5):
                    for _ in range(5):
                        vecs = [pp.Vector(random.choices(range(p), k=N)) for _ in range(R)]
                        rank = pp.Lattice(N, vecs).rank
                        relations = pp.relations_among(vecs)
                        nullity = relations.rank
                        self.assertEqual(rank + nullity, R)
                        for rel in relations.get_basis():
                            v = zero.copy()
                            for coeff, vec in zip(rel, vecs):
                                v += coeff * vec
                            self.assertEqual(v, zero)

    def test_coefficients_of(self):
        pp = Prime_u8(5)
        L = pp.Lattice(3, [[1, 1, 2], [0, 1, 2]])
        self.assertEqual(L.coefficients_of(pp.Vector([1, 1, 2])), pp.Vector([1, 0]))
        self.assertEqual(L.coefficients_of(pp.Vector([0, 1, 2])), pp.Vector([0, 1]))
        self.assertEqual(L.coefficients_of(pp.Vector([1, 2, 4])), pp.Vector([1, 1]))
        self.assertEqual(L.coefficients_of(pp.Vector([2, 4, 3])), pp.Vector([2, 2]))
        self.assertEqual(L.coefficients_of(pp.Vector([1, 0, 0])), pp.Vector([1, 4]))

    def test_linear_combination(self):
        pp = Prime_u8(5)
        L = pp.Lattice(3, [[1, 1, 2], [0, 1, 2]])
        self.assertEqual(L.linear_combination(pp.Vector([1, 0])), pp.Vector([1, 1, 2]))
        self.assertEqual(L.linear_combination(pp.Vector([0, 1])), pp.Vector([0, 1, 2]))
        self.assertEqual(L.linear_combination(pp.Vector([1, 1])), pp.Vector([1, 2, 4]))
        self.assertEqual(L.linear_combination(pp.Vector([2, 2])), pp.Vector([2, 4, 3]))
        self.assertEqual(L.linear_combination(pp.Vector([1, 4])), pp.Vector([1, 0, 0]))

    def test_coefficients_of_random(self):
        for p in (2, 3, 5):
            pp = Prime_u8(p)
            for N in range(10):
                for R in range(10):
                    for _ in range(5):
                        L = pp.Lattice(N, maxrank=R)
                        for _ in range(R):
                            L.add_vector(pp.Vector(random.choices(range(p), k=N)))
                        for _ in range(5):
                            w = pp.Vector(random.choices(range(p), k=L.rank))
                            v = L.linear_combination(w)
                            self.assertIn(v, L)
                            self.assertEqual(L.coefficients_of(v), w)
                        for _ in range(5):
                            v = pp.Vector(random.choices(range(p), k=N))
                            if v in L:
                                w = L.coefficients_of(v)
                                self.assertEqual(L.linear_combination(w), v)
                            else:
                                with self.assertRaises(ValueError):
                                    L.coefficients_of(v)

class Test_Mod2(unittest.TestCase):
    def test_construct_Prime2(self):
        Prime2()
        Prime2(2)

    def test_repr(self):
        self.assertEqual(repr(Prime2()), "Prime2()")
        self.assertEqual(repr(Prime2(2)), "Prime2()")
        self.assertEqual(repr(Prime2(2).Vector([])), "Prime2().Vector([])")
        self.assertEqual(repr(Prime2(2).Vector([0])), "Prime2().Vector([0])")
        self.assertEqual(repr(Prime2(2).Vector([1])), "Prime2().Vector([1])")
        self.assertEqual(repr(Prime2(2).Vector([0, 1, 1, 0])), "Prime2().Vector([0, 1, 1, 0])")

    def test_str(self):
        self.assertEqual(str(Prime2()), "Prime2()")
        self.assertEqual(str(Prime2(2)), "Prime2()")
        self.assertEqual(str(Prime2(2).Vector([])), "[]")
        self.assertEqual(str(Prime2(2).Vector([0])), "[0]")
        self.assertEqual(str(Prime2(2).Vector([1])), "[1]")
        self.assertEqual(str(Prime2(2).Vector([0, 1, 1, 0])), "[0110]")

    def test_construct_vectors(self):
        pp = Prime2()
        v = pp.Vector([])
        pp.Vector([0, 1, 1, 0])

    def test_inverse(self):
        pp = Prime2()
        self.assertEqual(pp.inverse(0), None)
        self.assertEqual(pp.inverse(1), 1)
        with self.assertRaises(OverflowError): pp.inverse(2)
        with self.assertRaises(OverflowError): pp.inverse(-1)

    def test_iter(self):
        sentinel = object()
        it = iter(Prime2().Vector([0, 1, 1, 0]))
        self.assertEqual(next(it), 0)
        self.assertEqual(next(it), 1)
        self.assertEqual(next(it), 1)
        self.assertEqual(next(it), 0)
        self.assertIs(next(it, sentinel), sentinel)

    def test_getitem(self):
        v = Prime2().Vector([0, 1, 0, 0])
        self.assertEqual([v[0], v[1], v[2], v[3]], [0, 1, 0, 0])
        self.assertEqual([v[-4], v[-3], v[-2], v[-1]], [0, 1, 0, 0])

    def test_setitem(self):
        v = Prime2().Vector([0, 0, 0, 0])
        self.assertEqual(v.tolist(), [0, 0, 0, 0])
        v[2] = 1
        self.assertEqual(v.tolist(), [0, 0, 1, 0])
        v[2] = 0
        self.assertEqual(v.tolist(), [0, 0, 0, 0])

    def test_vector_tolist(self):
        pp = Prime2()
        self.assertEqual(pp.Vector([]).tolist(), [])
        self.assertEqual(pp.Vector([0, 1, 1, 0]).tolist(), [0, 1, 1, 0])
        self.assertEqual(pp.Vector([1, 1, 0, 1]).tolist(), [1, 1, 0, 1])
        self.assertEqual(pp.Vector([0, 1] * 100).tolist(), [0, 1] * 100)
        self.assertEqual(pp.Vector([1, 0] * 100).tolist(), [1, 0] * 100)

    def test_vector_equality(self):
        pp = Prime2()
        lists = [
            [], [0], [1], [0, 0], [0, 1], [1, 0], [1, 1]
        ] + [
            [a] + [b]*1000 + [c]
            for a in range(2) for b in range(2) for c in range(2)
        ]
        for v in lists:
            self.assertEqual(pp.Vector(v).tolist(), v)
            for w in lists:
                self.assertEqual(pp.Vector(v) == pp.Vector(w), v == w)
                self.assertEqual(pp.Vector(v) != pp.Vector(w), v != w)

    def test_vector_tolist_random(self):
        random.seed(0)
        pp = Prime2()
        for N in range(300):
            data = [random.randrange(2) for _ in range(N)]
            self.assertEqual(pp.Vector(data).tolist(), data)

    def test_arithmetic(self):
        pp = Prime2()
        for N in [0, 1, 2, 20, 300]:
            zero = pp.Vector([0] * (5*N))
            A = pp.Vector([0, 1, 0, 1, 0] * N)
            B = pp.Vector([0, 0, 1, 1, 0] * N)
            C = pp.Vector([0, 1, 1, 0, 0] * N)
            self.assertEqual(A + B, C)
            self.assertEqual(A - B, C)
            self.assertEqual(A*0, zero)
            self.assertEqual(A*1, A)
            self.assertEqual(A*2, zero)
            self.assertEqual(A*3, A)
            A2 = A.copy()
            self.assertIsNot(A2, A)
            self.assertEqual(A2, A)
            A2 += B
            self.assertEqual(A2, C)
            A2 *= 1
            self.assertEqual(A2, C)
            A2 *= 0
            self.assertEqual(A2, zero)

    def test_shuffled_by_action(self):
        pp = Prime2()
        self.assertEqual(pp.Vector([0, 1]).shuffled_by_action(array(WORD, [0, 1])).tolist(), [0, 1])
        self.assertEqual(pp.Vector([0, 1]).shuffled_by_action(array(WORD, [1, 0])).tolist(), [1, 0])
        self.assertEqual(pp.Vector([1, 1]).shuffled_by_action(array(WORD, [1, 1])).tolist(), [0, 0])
        self.assertEqual(pp.Vector([0, 0, 0, 1, 1, 1]).shuffled_by_action(array(WORD, [0, 2, 4, 1, 3, 5])).tolist(), [0, 1, 0, 1, 0, 1])
        self.assertEqual(pp.Vector([0, 0, 0, 1, 1, 1]).shuffled_by_action(array(WORD, [0, 0, 0, 1, 1, 1])).tolist(), [0, 1, 0, 0, 0, 0])
        self.assertEqual(pp.Vector([0, 0, 0, 1, 1, 1]).shuffled_by_action(array(WORD, [0, 0, 0, 1, 1, 1]), 3).tolist(), [0, 1, 0])

class TestLatticeMod2(unittest.TestCase):
    def test_construct_Lattice(self):
        pp = Prime2()
        pp.Lattice(3)
        pp.Lattice(3, [])
        pp.Lattice(3, [[2, 1, 1]])
        pp.Lattice(3, maxrank=1)
        pp.Lattice(3, [[0, 1, 2], [1, 0, 0]], maxrank=-1)

    def test_tolist(self):
        pp = Prime2()
        self.assertEqual(pp.Lattice(3).tolist(), [])
        self.assertEqual(pp.Lattice(3, []).tolist(), [])
        self.assertEqual(pp.Lattice(3, [[0, 1, 1]]).tolist(), [[0, 1, 1]])
        self.assertEqual(pp.Lattice(3, [[1, 0, 0], [0, 1, 1]]).tolist(), [[1, 0, 0], [0, 1, 1]])

    def test_get_basis(self):
        pp = Prime2()
        self.assertEqual(pp.Lattice(3).get_basis(), [])
        self.assertEqual(pp.Lattice(3, []).get_basis(), [])
        self.assertEqual(pp.Lattice(3, [[0, 1, 1]]).get_basis(), [pp.Vector([0, 1, 1])])
        self.assertEqual(pp.Lattice(3, [[1, 0, 0], [0, 1, 1]]).get_basis(), [pp.Vector([1, 0, 0]), pp.Vector([0, 1, 1])])

    def test_full(self):
        pp = Prime2()
        self.assertEqual(pp.Lattice_full(0).tolist(), [])
        self.assertEqual(pp.Lattice_full(1).tolist(), [[1]])
        self.assertEqual(pp.Lattice_full(2).tolist(), [[1, 0], [0, 1]])
        self.assertEqual(pp.Lattice_full(3).tolist(), [[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    def test_maxrank(self):
        pp = Prime2()
        L = pp.Lattice(8, [[1]*8], maxrank=2)
        self.assertEqual(L.maxrank, 2)
        self.assertEqual(L.rank, 1)
        L.add_vector(pp.Vector([0, 0, 0, 0, 1, 1, 1, 1]))
        self.assertEqual(L.rank, 2)
        L.add_vector(pp.Vector([1, 1, 1, 1, 0, 0, 0, 0]))
        self.assertEqual(L.rank, 2)
        with self.assertRaises(IndexError):
            L.add_vector(pp.Vector([0, 1, 0, 0, 0, 0, 0, 0]))

    def test_row_echelon(self):
        pp = Prime2()
        self.assertEqual(pp.Lattice(3, [[1, 1, 1], [1, 1, 1]]).tolist(), [[1, 1, 1]])
        self.assertEqual(pp.Lattice(3, [[1, 1, 1], [1, 1, 0]]).tolist(), [[1, 1, 1], [0, 0, 1]])
        self.assertEqual(pp.Lattice(0, [[], [], []]).tolist(), [])
        self.assertEqual(pp.Lattice(1, []).tolist(), [])
        self.assertEqual(pp.Lattice(1, [[0]]).tolist(), [])
        self.assertEqual(pp.Lattice(1, [[1]]).tolist(), [[1]])
        self.assertEqual(pp.Lattice(1, [[1], [1]]).tolist(), [[1]])
        self.assertEqual(pp.Lattice(2, []).tolist(), [])
        self.assertEqual(pp.Lattice(2, [[0, 0]]).tolist(), [])
        self.assertEqual(pp.Lattice(2, [[1, 0]]).tolist(), [[1, 0]])
        self.assertEqual(pp.Lattice(2, [[1, 1]]).tolist(), [[1, 1]])
        self.assertEqual(pp.Lattice(2, [[0, 1]]).tolist(), [[0, 1]])
        self.assertEqual(pp.Lattice(2, [[0, 0], [1, 0]]).tolist(), [[1, 0]])
        self.assertEqual(pp.Lattice(2, [[1, 0], [0, 1]]).tolist(), [[1, 0], [0, 1]])
        self.assertEqual(pp.Lattice(2, [[1, 1], [1, 0]]).tolist(), [[1, 1], [0, 1]])
        self.assertEqual(pp.Lattice(2, [[0, 1], [1, 1]]).tolist(), [[1, 1], [0, 1]])

    def test_add_vector(self):
        pp = Prime2()
        L = pp.Lattice(4)
        self.assertEqual(L.tolist(), [])
        L.add_vector(pp.Vector([0, 0, 0, 1]))
        self.assertEqual(L.tolist(), [[0, 0, 0, 1]])
        L.add_vector(pp.Vector([0, 1, 1, 1]))
        self.assertEqual(L.tolist(), [[0, 1, 1, 1], [0, 0, 0, 1]])
        L.add_vector(pp.Vector([1, 1, 1, 1]))
        self.assertEqual(L.tolist(), [[1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 0, 1]])
        L.add_vector(pp.Vector([0, 0, 1, 1]))
        self.assertEqual(L.tolist(), [[1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1]])
        L.add_vector(pp.Vector([0, 1, 0, 1]))
        self.assertEqual(L.tolist(), [[1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 1, 1], [0, 0, 0, 1]])

    def test_copy(self):
        pp = Prime2()
        L = pp.Lattice(4, [[0, 1, 1, 1]])
        L2 = L.copy()
        self.assertEqual(L.tolist(), [[0, 1, 1, 1]])
        self.assertEqual(L2.tolist(), [[0, 1, 1, 1]])
        L2.add_vector(pp.Vector([0, 0, 1, 1]))
        self.assertEqual(L.tolist(), [[0, 1, 1, 1]])
        self.assertEqual(L2.tolist(), [[0, 1, 1, 1], [0, 0, 1, 1]])

    def test_reduce(self):
        pp = Prime2()
        L = pp.Lattice(4, [[1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 0, 1]])
        self.assertEqual(L.tolist(), [[1, 1, 1, 1], [0, 1, 1, 1], [0, 0, 0, 1]])
        L.reduce()
        self.assertEqual(L.tolist(), [[1, 0, 0, 0], [0, 1, 1, 0], [0, 0, 0, 1]])

    def test_reduce_large(self):
        pp = Prime2()
        for N in [0, 1, 2, 10, 100, 200]:
            L = pp.Lattice(N)
            for i in range(N):
                L.add_vector(pp.Vector([0] * i + [1] * (N - i)))
            self.assertEqual(L.tolist(), [[0] * i + [1] * (N - i) for i in range(N)])
            L.reduce()
            self.assertEqual(L.tolist(), [[0] * i + [1] + [0] * (N - 1 - i) for i in range(N)])

    def test_comparison_exhaustive(self):
        pp = Prime2()
        # get all row echelon forms
        vectors = list(map(list, itertools.product(range(2), repeat=4)))
        del vectors[0]
        pairs = [[u, v] for u, v in itertools.combinations(vectors, 2) if u.index(1) < v.index(1)]
        triples = [[u, v, w] for u, v in pairs for w in vectors if v.index(1) < w.index(1)]
        spaces = (
            [pp.Lattice(4)]
            + [pp.Lattice(4, [v]) for v in vectors]
            + [pp.Lattice(4, pair) for pair in pairs]
            + [pp.Lattice(4, triple) for triple in triples]
        )
        for a in spaces:
            for b in spaces:
                a_in_b = all(map(b.__contains__, a.get_basis()))
                b_in_a = all(map(a.__contains__, b.get_basis()))
                self.assertEqual(a == b, a_in_b and b_in_a)
                self.assertEqual(a != b, not (a_in_b and b_in_a))
                self.assertEqual(a <= b, a_in_b)
                self.assertEqual(a < b, a_in_b and not b_in_a)
                self.assertEqual(a >= b, b_in_a)
                self.assertEqual(a > b, b_in_a and not a_in_b)

    def test_str(self):
        pp = Prime2()
        self.assertEqual(str(pp.Lattice(0)), "<zero subspace of Z/2Z ^ 0>")
        self.assertEqual(str(pp.Lattice(1)), "<zero subspace of Z/2Z ^ 1>")
        self.assertEqual(str(pp.Lattice(2)), "<zero subspace of Z/2Z ^ 2>")
        self.assertEqual(str(pp.Lattice(1, [[1]])), "[1]\n")
        self.assertEqual(str(pp.Lattice(2, [[1, 0]])), "[10]\n")
        self.assertEqual(str(pp.Lattice(2, [[1, 1], [0, 1]])), "[11]\n[01]\n")

    def test_repr(self):
        pp = Prime2()
        self.assertEqual(repr(pp.Lattice(0, [])), "Prime(2).Lattice(0, [])")
        self.assertEqual(repr(pp.Lattice(1, [])), "Prime(2).Lattice(1, [])")
        self.assertEqual(repr(pp.Lattice(1, [[1]])), "Prime(2).Lattice(1, [[1]])")
        self.assertEqual(repr(pp.Lattice(2, [])), "Prime(2).Lattice(2, [])")
        self.assertEqual(repr(pp.Lattice(2, [[1, 1]])), "Prime(2).Lattice(2, [[1, 1]])")
        self.assertEqual(repr(pp.Lattice(2, [[1, 1], [0, 1]])), "Prime(2).Lattice(2, [[1, 1], [0, 1]])")

    def test_clear(self):
        pp = Prime2()
        for N in [0, 1, 2, 10, 64, 65, 128, 129]:
            L = pp.Lattice_full(N)
            basis = L.get_basis()
            self.assertIsNone(L.clear())
            self.assertEqual(L.get_basis(), [])
            for v in basis:
                L.add_vector(v) # make sure we can still add vectors

    def test_random(self):
        random.seed(123)
        pp = Prime2()
        for N in [0, 1, 2, 3, 4, 5, 10, 20, 64, 65, 128, 129]:
            for _ in range(20):
                mod = pp.Lattice(N)
                vecs = [
                    pp.Vector([random.randrange(2) for _ in range(N)])
                    for _ in range(random.randrange(N + 2))
                ]
                zero = pp.Vector([0] * N)
                vec0 = vecs[0] if vecs else zero
                for vec in vecs:
                    mod.add_vector(vec)
                    self.assertIn(vec, mod)
                    self.assertIn(zero-vec, mod)
                    self.assertIn(3*vec, mod)
                    self.assertIn(vec0 + vec, mod)
                    self.assertIn(3*vec0-4*vec, mod)

                if random.randrange(2):
                    copy = mod.copy()
                    mod.clear()
                    for v in copy.get_basis():
                        mod.add_vector(v)

                copy = mod.copy()
                copy.reduce()
                self.assertEqual(copy, mod)

                lin_combo = pp.Vector([0] * N)
                for vec in vecs:
                    lin_combo += random.randrange(2) * vec
                self.assertIn(lin_combo, mod)

                basis = mod.get_basis()
                for to_remove in basis:
                    submod = pp.Lattice(N)
                    for vec in basis:
                        if vec != to_remove:
                            submod.add_vector(vec)
                    self.assertNotIn(to_remove, submod)
                    also_outside = to_remove
                    for vec in submod.get_basis():
                        also_outside += random.randrange(2) * vec
                    self.assertNotIn(also_outside, submod)

    def test_decompose_simple(self):
        pp = Prime2()
        self.assertEqual(
            pp.Lattice(2, [[1, 0], [0, 1]]).decompose(),
            ([[0], [1]], [pp.Lattice(1, [[1]]), pp.Lattice(1, [[1]])])
        )
        self.assertEqual(
            pp.Lattice(2, [[1, 1], [0, 1]]).decompose(),
            ([[0], [1]], [pp.Lattice(1, [[1]]), pp.Lattice(1, [[1]])])
        )
        self.assertEqual(
            pp.Lattice(2, [[1, 1]]).decompose(),
            ([[0, 1]],
             [pp.Lattice(2, [[1, 1]])])
        )
        self.assertEqual(
            pp.Lattice(4, [[1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 1]]).decompose(),
            ([[0, 2, 3], [1]],
             [pp.Lattice(3, [[1, 0, 1], [0, 1, 1]]), pp.Lattice(1, [[1]])])
        )
        self.assertEqual(pp.Lattice(0).decompose(), ([], []))
        self.assertEqual(pp.Lattice(1).decompose(), ([[0]], [pp.Lattice(1)]))
        self.assertEqual(pp.Lattice(1, [[1]]).decompose(), ([[0]], [pp.Lattice(1, [[1]])]))
        self.assertEqual(
            pp.Lattice(10).decompose(),
            ([[i] for i in range(10)],
             [pp.Lattice(1, maxrank=0) for _ in range(10)])
        )
        self.assertEqual(
            pp.Lattice(5, [[0, 1, 0, 0, 0], [0, 0, 0, 1, 1]]).decompose(),
            ([[0], [1], [2], [3, 4]],
             [pp.Lattice(1), pp.Lattice(1, [[1]]), pp.Lattice(1), pp.Lattice(2, [[1, 1]])],
             )
        )
        self.assertEqual(
            pp.Lattice(5, [[0, 0, 1, 1, 1]]).decompose(),
            ([[0], [1], [2, 3, 4]],
             [pp.Lattice(1), pp.Lattice(1), pp.Lattice(3, [[1, 1, 1]])])
        )
        self.assertEqual(
            pp.Lattice(10, [[0,1,0,0,0,0,0,0,0,1],
                            [0,0,1,0,0,0,0,0,0,0],
                            [0,0,0,1,0,0,0,0,0,1],
                            [0,0,0,0,1,0,1,0,0,0],
                            [0,0,0,0,0,1,0,0,0,0],
                            [0,0,0,0,0,0,0,1,0,0]]).decompose(),
            ([[0], [1, 3, 9], [2], [4, 6], [5], [7], [8]],
             [pp.Lattice(1),
              pp.Lattice(3, [[1, 0, 1], [0, 1, 1]]),
              pp.Lattice(1, [[1]]),
              pp.Lattice(2, [[1, 1]]),
              pp.Lattice(1, [[1]]),
              pp.Lattice(1, [[1]]),
              pp.Lattice(1),])
        )

    def test_decompose_exhaustive(self):
        pp = Prime2()
        vectors = [pp.Vector(list(v)) for v in itertools.product(range(2), repeat=4)]
        lattices = (
            [pp.Lattice(4, maxrank=0), pp.Lattice_full(4)]
            + [pp.Lattice(4, [v], maxrank=1) for v in vectors]
            + [pp.Lattice(4, [v, w], maxrank=2) for v, w in itertools.combinations(vectors, 2)]
            + [pp.Lattice(4, [u, v, w], maxrank=3) for u, v, w in itertools.combinations(vectors, 3)]
        )
        for L in lattices:
            indexes, summands = L.decompose()
            self.assertEqual(len(indexes), len(summands))
            self.assertEqual(sum(map(len, indexes)), 4)
            self.assertEqual({x for arr in indexes for x in arr}, set(range(4)))
            self.assertEqual(sum(summand.rank for summand in summands), L.rank)
            self.assertEqual(sum(summand.ambient_dimension for summand in summands), L.ambient_dimension)
            for summand in summands:
                self.assertEqual(summand.maxrank, summand.rank)
                self.assertEqual(summand.decompose(), ([list(range(summand.ambient_dimension))], [summand]))
            reconstruction = pp.Lattice(4)
            for arr, summand in zip(indexes, summands):
                self.assertEqual(summand.ambient_dimension, len(arr))
                for vec0 in summand.get_basis():
                    vec = pp.Vector([0] * 4)
                    for i, x in enumerate(arr):
                        vec[x] = vec0[i]
                    reconstruction.add_vector(vec)
            self.assertEqual(L, reconstruction)

    def test_decompose_keeping_together(self):
        pp = Prime2()
        L = pp.Lattice(7, [[1,1,0,0,0,0,0],
                           [0,0,1,0,0,1,0],
                           [0,0,0,1,1,0,0]])
        self.assertEqual(L.decompose([]),
            ([[0, 1], [2, 5], [3, 4], [6]],
             [pp.Lattice(2, [[1, 1]]),
              pp.Lattice(2, [[1, 1]]),
              pp.Lattice(2, [[1, 1]]),
              pp.Lattice(1)])
        )
        self.assertEqual(L.decompose([[1, 5, 4, 6]]), ([[0,1,2,3,4,5,6]], [L]))
        self.assertEqual(L.decompose([[1, 5]]),
            ([[0, 1, 2, 5], [3, 4], [6]],
             [pp.Lattice(4, [[1, 1, 0, 0], [0, 0, 1, 1]]),
              pp.Lattice(2, [[1, 1]]),
              pp.Lattice(1)])
        )
        self.assertEqual(L.decompose([[0, 4]]),
            ([[0, 1, 3, 4], [2, 5], [6]],
             [pp.Lattice(4, [[1, 1, 0, 0], [0, 0, 1, 1]]),
              pp.Lattice(2, [[1, 1]]),
              pp.Lattice(1)])
        )
        self.assertEqual(L.decompose([[1, 6]]),
            ([[0, 1, 6], [2, 5], [3, 4]],
             [pp.Lattice(3, [[1, 1, 0]]),
              pp.Lattice(2, [[1, 1]]),
              pp.Lattice(2, [[1, 1]])])
        )
        self.assertEqual(L.decompose([[0, 6], [5, 4]]),
            ([[0, 1, 6], [2, 3, 4, 5]],
             [pp.Lattice(3, [[1, 1, 0]]),
              pp.Lattice(4, [[1, 0, 0, 1], [0, 1, 1, 0]])])
        )
        self.assertEqual(L.decompose([[0, 3], [2, 6]]),
            ([[0, 1, 3, 4], [2, 5, 6]],
             [pp.Lattice(4, [[1, 1, 0, 0], [0, 0, 1, 1]]),
              pp.Lattice(3, [[1, 1, 0]])])
        )
        self.assertEqual(L.decompose([[0, 2], [3, 6]]),
            ([[0, 1, 2, 5], [3, 4, 6]],
             [pp.Lattice(4, [[1, 1, 0, 0], [0, 0, 1, 1]]),
              pp.Lattice(3, [[1, 1, 0]])])
        )
        self.assertEqual(L.decompose([[1, 3, 6]]),
            ([[0, 1, 3, 4, 6], [2, 5]],
             [pp.Lattice(5, [[1, 1, 0, 0, 0], [0, 0, 1, 1, 0]]),
              pp.Lattice(2, [[1, 1]])])
        )

    def test_relation_among(self):
        pp = Prime2()
        self.assertEqual(pp.relations_among([]), pp.Lattice(0))
        self.assertEqual(pp.relations_among([pp.Vector([])]), pp.Lattice(1, [[1]]))
        self.assertEqual(pp.relations_among([pp.Vector([]), pp.Vector([])]), pp.Lattice(2, [[1, 0], [0, 1]]))
        self.assertEqual(pp.relations_among([pp.Vector([1])]), pp.Lattice(1))
        self.assertEqual(pp.relations_among([pp.Vector([0])]), pp.Lattice(1, [[1]]))
        self.assertEqual(pp.relations_among([pp.Vector([0]), pp.Vector([0])]), pp.Lattice(2, [[1, 0], [0, 1]]))
        self.assertEqual(pp.relations_among([pp.Vector([0]), pp.Vector([1])]), pp.Lattice(2, [[1, 0]]))
        self.assertEqual(pp.relations_among([pp.Vector([1]), pp.Vector([0])]), pp.Lattice(2, [[0, 1]]))
        self.assertEqual(pp.relations_among([pp.Vector([1]), pp.Vector([1])]), pp.Lattice(2, [[1, 1]]))
        self.assertEqual(pp.relations_among([pp.Vector([0, 0])]), pp.Lattice(1, [[1]]))
        self.assertEqual(pp.relations_among([pp.Vector([0, 1])]), pp.Lattice(1))
        self.assertEqual(pp.relations_among([pp.Vector([1, 0])]), pp.Lattice(1))
        self.assertEqual(pp.relations_among([pp.Vector([1, 1])]), pp.Lattice(1))
        self.assertEqual(pp.relations_among([pp.Vector([1, 0]), pp.Vector([0, 1])]), pp.Lattice(2))
        self.assertEqual(pp.relations_among([pp.Vector([1, 0]), pp.Vector([0, 1]), pp.Vector([1, 1])]), pp.Lattice(3, [[1, 1, 1]]))

    def test_rank_nullity(self):
        random.seed(0)
        pp = Prime2()
        for N in [0, 1, 2, 3, 4, 63, 64, 65, 100, 127, 128, 129]:
            zero = pp.Vector([0] * N)
            for R in range(5):
                for _ in range(5):
                    vecs = [pp.Vector(random.choices(range(2), k=N)) for _ in range(R)]
                    rank = pp.Lattice(N, vecs).rank
                    relations = pp.relations_among(vecs)
                    nullity = relations.rank
                    self.assertEqual(rank + nullity, R)
                    for rel in relations.get_basis():
                        v = zero.copy()
                        for coeff, vec in zip(rel, vecs):
                            v += coeff * vec
                        self.assertEqual(v, zero)



if __name__ == "__main__":
    unittest.main()

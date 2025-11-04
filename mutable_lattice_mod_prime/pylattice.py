from functools import lru_cache
from bisect import bisect_left
import itertools

@lru_cache(maxsize=5)
def assert_prime(p):
    if p in {2,3,5,7}:
        return
    if p <= 10:
        raise ValueError(f"{p} is not prime")
    if p % 2 == 0:
        return False
    assert p.bit_length() <= 32
    # deterministic miller-rabin
    d = p - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 7, 61):
        x = pow(a, d, p)
        for _ in range(s):
            y = (x * x) % p
            if y == 1 and x != 1 and x != p - 1:
                raise ValueError(f"{p} is not prime")
            x = y
        if y != 1:
            raise ValueError(f"{p} is not prime")

# notes on bounds
# want: a + b * c < 2**64
# suffices if m:=max(a,b,c)
# 	satisfies m + m * m < 2**64
# it works if m <= 4_294_967_295=2**32-1
# because (2**32-1) + (2**32-1)**2
#         = 2**32 - 1 + 2**64 - 2**33 + 1
#         = 2**64 - 2**32

@lru_cache(maxsize=5)
def get_inverse_table(p):
    assert p <= 256
    inverse = [None] * p
    inverse[1] = 1
    for x in range(2, p):
        if inverse[x] is not None:
            continue
        cycle = [1, x]
        while cycle[-1] != 1:
            cycle.append((cycle[-1] * x) % p)
        for x, y in zip(cycle, reversed(cycle)):
            inverse[x] = y
    return tuple(inverse)

assert get_inverse_table(2) == (None, 1)
assert get_inverse_table(3) == (None, 1, 2)
assert get_inverse_table(5) == (None, 1, 3, 2, 4)
assert get_inverse_table(7) == (None, 1, 4, 5, 2, 3, 6)
assert get_inverse_table(11) == (None, 1, 6, 4, 3, 9, 2, 8, 7, 5, 10)

class LatticeModPrime:

    __slots__ = ["p",
                 "N",
                 "basis",
                 "pivot_location_in_column",
                 "pivot_location_in_row",
                 "zero_columns",
                 "inverses"]

    def __init__(self, p, ambient_dimension):
        assert_prime(p)
        self.p = p
        self.N = N = ambient_dimension
        self.pivot_location_in_column = [None] * N
        self.pivot_location_in_row = []
        self.basis = []
        self.zero_columns = list(range(N))
        self.inverses = get_inverse_table(p)

    def copy(self):
        other = object.__new__(__class__)
        other.N = self.N
        other.pivot_location_in_column = self.pivot_location_in_column.copy()
        other.pivot_location_in_row = self.pivot_location_in_row.copy()
        other.basis = list(map(list.copy, self.basis))
        other.zero_columns = self.zero_columns.copy()
        return other

    def __contains__(self, vec):
        zc = self.zero_columns
        for _ in filter(vec.__getitem__, zc):
            return False
        col_piv = self.pivot_location_in_column
        N = self.N
        basis = self.basis
        vec = vec.copy()
        p = self.p
        for j in filter(vec.__getitem__, range(N)):
            # assert vec[j] % p != 0
            i = col_piv[j]
            if i is None:
                return False
            row = basis[i]
            assert row[j] == 1
            b = vec[j]
            for jj in range(j, N):
                # in c, remember not to subtract!
                vec[jj] = (vec[jj] - b * row[jj]) % p
        return True

    def add_vector(self, vec0):
        col_piv = self.pivot_location_in_column
        row_piv = self.pivot_location_in_row
        N = self.N
        basis = self.basis
        assert len(vec0) == N
        p = self.p
        vec = vec0.copy()

        # TODO: benchmark some possibilities:
        # (1) This stuff, not reduced
        #       (a) with the if-statements
        #       (b) no if-statements
        # (2) This stuff, reduced
        #       (a) with the if-statements
        #       (b) no if-statements
        # (3) For each row, reduce, then check for nonzero index, then insert.
        #       (a) with the if-statements
        #       (b) no if-statements

        for j in filter(vec.__getitem__, range(N)):
            # if vec[j] % p == 0:
            #     continue
            i = col_piv[j]
            if i is None:
                # assert set(map(type, row_piv)) <= {int}, row_piv
                # assert type(j) is int
                where = bisect_left(row_piv, j)
                basis.insert(where, vec)
                row_piv.insert(where, j)
                col_piv[j] = where
                for ii in range(where + 1, len(basis)):
                    assert col_piv[row_piv[ii]] == ii - 1
                    col_piv[row_piv[ii]] = ii                
                zc = self.zero_columns
                s = bisect_left(zc, j)
                zc[s:] = itertools.filterfalse(vec0.__getitem__, zc[s:])

                # assert vec[j] != 0
                # assert vec[j] < p
                # inv = pow(vec[j], p - 2, p)
                inv = self.inverses[vec[j]]
                # assert (inv * vec[j]) % p == 1, (inv, vec[j], p)
                for jj in filter(vec.__getitem__, range(j, N)):
                    vec[jj] = (vec[jj] * inv) % p
                # assert vec[j] == 1, vec[j]
                if 0:
                    # *reduced* row echelon form
                    for j in filter(vec.__getitem__, range(N)):
                        if i is not None:
                            row = basis[i]
                            # assert row[j] == 1
                            b = vec[j]
                            for jj in filter(row.__getitem__, range(j, N)):
                                # in c, remember not to subtract!
                                vec[jj] = (vec[jj] - b * row[jj]) % p
                    for ii in range(where):
                        row = basis[ii]
                        b = row[j]
                        if b:
                            for jj in filter(vec.__getitem__, range(j, N)):
                                row[jj] = (row[jj] - b * vec[jj]) % p
                            # assert not row[j]
                return
            row = basis[i]
            # assert row[j] == 1
            b = vec[j]
            for jj in filter(row.__getitem__, range(j, N)):
                # optional: check if row[jj] is nonzero
                # in c, remember not to subtract!
                vec[jj] = (vec[jj] - b * row[jj]) % p


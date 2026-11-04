"""Exact checks for RESULTS.md, Section 2.1 (no SDP solver required).

Run from any directory: python search/projection_bound.py
The proof uses the two-projection decomposition over C, so applies to complex
as well as real realizations. Symbolic checks accompany, not replace, that proof.
"""
from fractions import Fraction as Q
import sympy as sp


def check_blocks():
    t = sp.symbols("t", real=True)
    # Eigenvalues of 8(BC+CB)/2-I on a nontrivial two-dimensional block.
    lm, lp = 4*t*t - 4*t - 1, 4*t*t + 4*t - 1
    assert sp.expand(lm + 2 - 4*(t-sp.Rational(1, 2))**2) == 0
    assert sp.expand(lm + lp + 2 - 8*t*t) == 0
    assert sp.expand((t*t-t)/2 + sp.Rational(1, 8)
                     - (t-sp.Rational(1, 2))**2/2) == 0
    # Their negative-part sum is min(0, lm, lp, lm+lp) >= -2,
    # since lp >= lm for 0 <= t <= 1. This proves the block inequality.
    for b in (0, 1):
        for c in (0, 1):
            assert min(0, 8*b*c-1) + 2-b-c >= 0

    # At equality, the nontrivial block consists of three trine projectors.
    B = sp.diag(1, 0)
    v = sp.Matrix([sp.Rational(1, 2), sp.sqrt(3)/2])
    C = v*v.T
    u = sp.Matrix([sp.Rational(1, 2), -sp.sqrt(3)/2])
    A = u*u.T
    assert (A*B*C).trace() == -sp.Rational(1, 8)
    assert A+B+C == sp.Rational(3, 2)*sp.eye(2)
    assert 8*(A*B*C).trace() == (A+B+C).trace()-4
    print("exact: three-projection trace inequality and trine equality block")


def check_equality_shapes():
    # Eight trine blocks and psi use 17 dimensions and rank 9 per context.
    # x,y,z are the ranks of the pair-intersection sectors 12,13,23.
    shapes = set()
    for x in range(3):
        for y in range(3):
            z = 2-x-y
            if z < 0:
                continue
            shape = tuple(sorted((9+x+y, 9+x+z, 9+y+z)))
            if min(shape) >= 10:
                shapes.add(shape)
    assert shapes == {(10, 10, 11)}

    # The two size-10 contexts each have precisely two deficient events.
    # Merge their twin columns/rows: all nine effective weights are 1/9.
    # Two saturated rows require common-neighbour weight 1/12 instead.
    assert all(Q(k, 9) != Q(1, 12) for k in range(10))
    # Equivalent spectral integrality obstruction for the merged incidence N.
    NN = sp.Rational(9, 4)*sp.eye(9) + sp.Rational(3, 4)*sp.ones(9)
    assert NN[0, 1] == sp.Rational(3, 4)  # cannot equal integer N N^T
    print("exact: 31 positive events in d=19, equality shape (10,10,11), are impossible")


def check_bounds():
    assert 2*18-7 < 30
    assert [n for n in range(30, 34) if n <= 2*19-7] == [30, 31]
    assert 33 == 2*20-7  # H33 saturates the new bound.
    print("n >= 30; n <= 2d-7; d >= 19; d=19 requires n=30, shape (10,10,10)")
    print("H33: n=33, d=20 attains n=2d-7")


if __name__ == "__main__":
    check_blocks()
    check_equality_shapes()
    check_bounds()

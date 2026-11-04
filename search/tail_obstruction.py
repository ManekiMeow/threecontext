"""Exact algebra behind the exclusion of the rotated-tail d=19 GHZ ansatz.

This checks the matrix identities used in the proof, its support-counting
constants, and the independent no-go for nine identical two-dimensional
projection blocks. The combinatorial arguments are analytic, not a numerical
search over graphs. Run from any directory with Python + SymPy.
"""
from fractions import Fraction as Q
import sympy as sp


def tail_geometry():
    # Canonical trines: A=B-C, B=(1,0), C=(1/2,-sqrt(3)/2).
    a = sp.Matrix([sp.Rational(1, 2), sp.sqrt(3)/2])
    b = sp.Matrix([1, 0])
    c = sp.Matrix([sp.Rational(1, 2), -sp.sqrt(3)/2])
    assert a*a.T+b*b.T+c*c.T == sp.Rational(3, 2)*sp.eye(2)
    assert (a.T*b)[0] == sp.Rational(1, 2)
    assert (a.T*c)[0] == -sp.Rational(1, 2)
    u, v = sp.Rational(3, 5), sp.Rational(4, 5)
    ta, tb, tc = sp.Matrix([u, v]), sp.Matrix([1, 0]), sp.Matrix([0, 1])
    E = ta*ta.T-sp.eye(2)/2
    assert E*E == sp.eye(2)/4
    assert (ta.T*E*tb)[0] == u/2
    assert (ta.T*E*tc)[0] == v/2
    assert (tb.T*E*tc)[0] == u*v

    # In the full context cross matrix Y, the nine shared-sector singular
    # values are 1 and eight copies of 1/2; the tail contributes u.
    Y = sp.diag(1, *([sp.Rational(1, 2)]*8), u)
    p, q = sp.eye(10)[:, 0], sp.eye(10)[:, 9]
    assert Y*Y.T == sp.eye(10)/4+3*p*p.T/4+(u*u-sp.Rational(1, 4))*q*q.T
    print("exact: trine/tail edge coefficients and cross-Gram decomposition")


def disjoint_support_obstruction():
    # Subcubic row and column supports leave at least 3 disjoint other rows.
    assert 10-1-3*(3-1) == 3
    # Every bipartite component of minimum degree 3 has >=6 vertices, so
    # a bipartite graph with these parameters on 10 vertices is connected.
    assert 2*(2*3) > 10
    z = sp.symbols("z", real=True)  # z=u^2 in [1/2,1)
    K = 3/(4*z-1)
    assert sp.simplify(K-1-4*(1-z)/(4*z-1)) == 0
    P, R, S = sp.symbols("P R S", real=True)
    Qmass = 1-P
    # Orthogonality imposes PR=QS, and unit norms impose PR²+QS²=1.
    # Therefore RS=1, contradicting the disjoint-support requirement RS=K>1.
    assert sp.expand(P*R**2+Qmass*S**2-R*S
                     -(R-S)*(P*R-Qmass*S)) == 0
    print("exact: disjoint-row degree bound, connectivity, and RS=1 contradiction")


def identical_block_obstruction():
    # A triple of rank-one projections in C^2 has det Gram=0. With Bargmann
    # invariant -1/9, the sum of squared pair overlaps is 7/9.
    pair_sum = 1+2*Q(-1, 9)
    defect_square = pair_sum-Q(3, 4)
    assert defect_square == Q(1, 36)
    c = sp.symbols("c")
    assert sp.solve((1-c)/36-(1-9*c)/4, c) == [sp.Rational(1, 10)]
    t2 = sp.symbols("t2", positive=True)
    kappa = sp.Rational(1, 2)-1/(9*t2)
    assert sp.simplify((3-2*kappa)/(1-2*kappa)) == 1+9*t2
    # A k-regular support incidence matrix M must satisfy
    # (M M^T)_ij = k(k-1)/9 for distinct i,j.
    degrees = [k for k in range(1, 11) if (k*(k-1)) % 9 == 0]
    assert degrees == [1, 9, 10]
    cosines_squared = {k: Q(10-k, 9*k) for k in degrees}
    assert cosines_squared == {1: Q(1), 9: Q(1, 81), 10: Q(0)}
    # Generic nonparallel, nonorthogonal lines leave only k=9 at each pair.
    # Their Bargmann invariant would have modulus (1/9)^3, not 1/9.
    assert Q(1, 9)**3 != Q(1, 9)
    print("exact: identical-block model forces incompatible 9-regular supports")


if __name__ == "__main__":
    tail_geometry()
    disjoint_support_obstruction()
    identical_block_obstruction()

"""Exact checks for the pairwise-isoclinic obstruction in dimension 19.

No (10,10,10) GHZ realization in C^19 can have all three pairs of the
reduced context spaces S_k intersect psi-perp isoclinic. This includes every
construction made from nine identical two-dimensional triples of lines.

Run with standard Python plus SymPy:
    python search/construction30_isoclinic.py

The proof is analytic; the checks below verify its symbolic identities and the
small integer obstruction. They are not a numerical search for realizations.
"""
from fractions import Fraction as Q
import sympy as sp


def odd_rank_reduction():
    # Gauge two isoclinic pairs as A=[I;0], B=[u I;s I], C=[v I;t W],
    # with u,v,s,t positive and W unitary. Write a=u*v, b=s*t.
    # The remaining pair has overlap a I+b W. Isoclinicity gives W+W*=h I.
    # Its eigenvalues are exp(+/- i theta). The imaginary part of tr(ABC) is
    # a*b*(2*m-9)*sin(theta), where m counts the plus eigenvalues.
    coefficients = [2*m-9 for m in range(10)]
    assert 0 not in coefficients
    # Therefore the real trace -1 forces sin(theta)=0, and W=+I or -I.
    # W=+I makes tr(ABC)=9*a*(a+b)>0, so only W=-I can remain.
    print("exact: odd reduced rank 9 forces W = -I in the isoclinic gauge")


def homogeneous_spectrum():
    x, y, z, lam, c = sp.symbols("x y z lam c", real=True)
    line_gram = sp.Matrix([[1, x, y], [x, 1, z], [y, z, 1]])
    assert sp.expand(line_gram.det() - (1-x*x-y*y-z*z+2*x*y*z)) == 0
    # Line triple product xyz=-1/9 and determinant zero imply sum cos^2=7/9.
    sum_squares = Q(1) + 2*Q(-1, 9)
    assert sum_squares == Q(7, 9)
    # The sum of the three line projections has trace 3 and determinant 20/9.
    determinant = Q(3)-sum_squares
    assert determinant == Q(20, 9)
    polynomial = lam**2-3*lam+sp.Rational(20, 9)
    assert sp.expand(polynomial-(lam-sp.Rational(4, 3))*(lam-sp.Rational(5, 3))) == 0
    # Hence E^2=I/36 on psi-perp. Compare with the per-event identity, Lemma 3.
    probability = sp.solve(sp.Eq((1-c)/36, (1-9*c)/4), c)
    assert probability == [sp.Rational(1, 10)]
    print("exact: F|psi-perp has eigenvalues 4/3,5/3; every probability is 1/10")


def pair_support_obstruction():
    s, entry, k = sp.symbols("s entry k", positive=True)
    # If X is one pair's overlap matrix and s its reduced principal cosine^2,
    # the third-context expansion on an edge gives
    # 0 = 1/10 - (entry-1/10)/(9*s).
    nonzero_entry = sp.solve(sp.Eq(0, sp.Rational(1, 10)
                                  -(entry-sp.Rational(1, 10))/(9*s)), entry)[0]
    assert sp.simplify(nonzero_entry-(1+9*s)/10) == 0
    # X 1 = 1, so each row/column has k=10/(1+9s) ones in its support M.
    s_of_k = (10-k)/(9*k)
    assert sp.simplify(nonzero_entry.subs(s, s_of_k)-1/k) == 0
    # XX*=s I+(1-s)J/10 and X=M/k give the common-neighbour count lambda.
    common = sp.simplify(k*k*(1-s_of_k)/10)
    assert sp.simplify(common-k*(k-1)/9) == 0
    allowed = [j for j in range(2, 10) if Q(j*(j-1), 9).denominator == 1]
    assert allowed == [9]
    forced_s = Q(10-9, 9*9)
    assert 3*forced_s == Q(1, 27) != Q(7, 9)
    print("exact: each pair support would be a symmetric 2-(10,k,k(k-1)/9) design")
    print("exact: 2 <= k <= 9 permits only k=9; sum principal cosine^2=1/27 != 7/9")
    # Claude's concrete spread model is excluded even before the design test.
    assert 10/(1+9*Q(1, 3)) == Q(5, 2)
    print("specific cosines (1/sqrt3,1/3,1/sqrt3): two pair degrees would be 5/2")


if __name__ == "__main__":
    odd_rank_reduction()
    homogeneous_spectrum()
    pair_support_obstruction()
    print("No pairwise-isoclinic (10,10,10) GHZ realization exists in complex dimension 19.")

"""Exact obstructions to uniform 30-event constructions.

Run: python search/event30_hoffman.py

The checks accompany the proofs supplied to RESULTS.md; no SDP tolerance or
graph enumeration is involved. In particular they do NOT exclude arbitrary
weighted or complex (10,10,10) realizations.
"""
from fractions import Fraction as Q
import sympy as sp


def degree_pairs():
    # Uniform c=1/10, neighbourhood weights k/10,l/10, residual (10-k-l)/10.
    # Lemma 14 gives 1/k + 1/l + 1/(10-k-l) <= 1.
    allowed = []
    saturated = []
    for k in range(1, 10):
        for l in range(1, 10-k):
            value = Q(1, k) + Q(1, l) + Q(1, 10-k-l)
            if value <= 1:
                allowed.append((k, l))
            if value == 1:
                saturated.append((k, l))
    assert allowed == [(2, 4), (3, 3), (3, 4), (4, 2), (4, 3), (4, 4)]
    assert saturated == [(2, 4), (4, 2), (4, 4)]
    print("uniform c=1/10: allowed cross-degree pairs", allowed)
    print("three-way Cauchy-Schwarz equality pairs", saturated)


def unweighted_degree_six():
    # H has n=30, degree6, three neighbours in each other class, and
    # G=I+A/3 is PSD. P=(A+3I)(2A-3I)^2 is PSD. Its Perron eigenvalue
    # is729 and its diagonal is27. Removing (729/30)J leaves diagonal27/10,
    # hence 108/5 <= P_ij <=27 for i!=j.
    x = sp.symbols("x")
    assert sp.expand((x+3)*(2*x-3)**2) == 4*x**3 - 27*x + 27
    center, radius = Q(729, 30), Q(27)-Q(729, 30)
    lo, hi = center-radius, center+radius
    assert (lo, hi) == (Q(108, 5), Q(27))
    edge_walks = [w for w in range(100) if lo <= 4*w-27 <= hi]
    nonedge_walks = [w for w in range(100) if lo <= 4*w <= hi]
    assert edge_walks == [13] and nonedge_walks == [6]
    # Thus A^3 = 7A - 6I + 6J. Its diagonal entries are zero as required.
    # (A^4)_ii=7*6+6*6=78. Subtract (A^2)_ii^2=36 to obtain
    # sum of squared common-neighbour counts over the 23 nonneighbours.
    count = 30-1-6
    first_moment = 6*(6-1)
    fourth_diag = 7*6 + 6*6
    second_moment = fourth_diag - 6**2
    integer_lower_bound = 3*first_moment - 2*count
    assert (count, first_moment, second_moment, integer_lower_bound) == (23, 30, 42, 44)
    assert second_moment < integer_lower_bound
    # Every integer mu has (mu-1)(mu-2)>=0, so sum mu^2>=3 sum mu-2 count.
    print("degree6 unweighted Gram: A^3=7A-6I+6J would force")
    print("  23 integer common-neighbour counts with sum30, sum-squares42 >=44: impossible")


def degree_four_per_context():
    # This step needs uniform event probabilities and cross-degrees(4,4),
    # but not any initial assumption on the phases/weights of Gram entries.
    c, alpha, beta, rho = map(sp.Rational, ("1/10", "2/5", "2/5", "1/5"))
    assert c*(1/alpha + 1/beta + 1/rho) == 1
    positive = sp.sqrt(c)/alpha
    negative = -sp.sqrt(c)/rho
    assert sp.simplify(positive-sp.sqrt(10)/4) == 0
    assert sp.simplify(negative+2*positive) == 0
    # Equality forces Pg=positive*Ppsi, Qg=positive*Qpsi,
    # Rg=negative*Rpsi. Consequently Fg=(4/3)g+(sqrt10/6)psi.
    frame_scalar = sp.Rational(4, 3)
    state_scalar = sp.simplify((2*positive/3)/sp.sqrt(c))
    assert state_scalar == sp.Rational(5, 3)
    dimension = sp.simplify((30-state_scalar)/frame_scalar)
    assert dimension == sp.Rational(85, 4) and dimension.q != 1
    print("uniform cross-degrees(4,4): F=(4/3)I+(5/3)|psi><psi|")
    print("  tr(F)=30 forces dimension85/4: impossible, over C as well as R")


if __name__ == "__main__":
    degree_pairs()
    unweighted_degree_six()
    degree_four_per_context()
    print("All regular unweighted Hoffman constructions on (10,10,10) are excluded.")

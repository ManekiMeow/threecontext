"""Angle-dependent probabilities for one isoclinic pair (Theorem 58).

These are necessary conditions for thirty positive events in contexts
(10,10,10), with triangle-free allowed support, complex vectors, and ANY
ambient dimension. They neither construct a GHZ realization nor exclude
every nonuniform isoclinic pair. Allowed edges may have zero actual overlap.

Let E_A,E_B be the full rank-ten context projections, and suppose their
reduced rank-nine spaces in psi-perp are isoclinic with cosine u. The endpoints
u=0 and u=1 are impossible: the first makes every A--B overlap sqrt(c_a c_b)
nonzero, leaving no positive C event; the second gives E_A=E_B, whereas
triangle-freeness makes E_A g_c and E_B g_c orthogonal, forcing both zero.
Thus 0<u<1. Write t=u(1-u) and

    T=(E_A E_B+E_B E_A)/2.

Its spectral decomposition has eigenvalue 1 on psi, eigenvalues

    lambda_minus=(u*u-u)/2=-t/2,
    lambda_plus =(u*u+u)/2,

each ninefold, and zero on the orthogonal complement of the two-context span.
Let Q_minus,Q_plus,Q_zero denote these three projections on psi-perp.
Triangle-freeness of the allowed support gives

    <g_c,T g_c>=Re <E_A g_c,E_B g_c>=0

for each C event. This remains true when allowed edges happen to vanish.

INDIVIDUAL PROBABILITIES AND ANGLES
----------------------------------
For c=|<psi,g_c>|^2, the spectral minimum gives

    0 >= c-t*(1-c)/2,
    c <= q(u):=t/(2+t) <= 1/9.

Equality c=q(u) holds exactly when g_c-sqrt(c)*psi lies in Q_minus, after
positive-state rephasing. Summing the ten probability bounds gives t>=2/9,
equivalently 1/3<=u<=2/3. At either endpoint q=1/10, normalization makes
all C probabilities 1/10 and every reduced C vector lies in Q_minus. The
rank-nine reduced C space is therefore Q_minus. Its cosines with both A
and B are sqrt((1-u)/2), so all three pairs are isoclinic and the entire
three-context span is nineteen-dimensional. Theorem 36 excludes this case.
Consequently

    1/3 < u < 2/3,
    1-9*q(u) < c <= q(u),
    1-9*q(u)=2*(2*u-1)^2/(2+t).

The lower bound is strict. If equality held, the other nine C events would
all attain q(u); their reduced vectors are independent because their Gram
matrix is I-ww*, with determinant equal to the omitted positive probability.
They would span Q_minus and the whole reduced C space, forcing every C event
to attain q(u), hence q=1/10 and the already excluded endpoints. At u=1/2
the lower bound is simply positivity; Corollary 32 separately makes the
upper bound 1/9 strict there.

The strict angle interval was previously known in dimension nineteen. Its
extension to arbitrary ambient dimension and the individual probability
bounds here do not rely on uniform probabilities in A or B.

QUANTITATIVE SPECTRAL CONCENTRATION
---------------------------------
Put R_C=E_C-|psi><psi|, a rank-nine projection, and define the nonnegative
traces

    a=tr(Q_minus R_C), b=tr(Q_plus R_C), r=tr(Q_zero R_C).

Then a+b+r=9 and the sum of the ten zero expectations gives

    1+lambda_minus*a+lambda_plus*b=0,
    b+(1-u)*r/2 = s(u):=9*(1-u)/2-1/u.

Therefore

    0 <= b <= s(u) <= 9/2-3*sqrt(2),
    0 <= r <= 9-2/t <= 1,
    a >= 2/t >= 8.

The maximum for s follows from 9*u/2+1/u>=3*sqrt(2), with equality at
u=sqrt(2)/3. In dimension nineteen, Q_zero=0, so r=0 and b=s(u) exactly.
In that case ||R_C-Q_minus||_F^2=2*b<=9-6*sqrt(2). In larger ambient
dimension the extra term r must be retained. A small projection trace is
not a bound on the rank of an intersection or of a residual subspace.

The exact controls below verify the identities, interval implications,
spectral blocks and equality examples. Local spectral equality controls
are not GHZ examples; no realization is inferred from these necessary
scalar conditions. No numerical optimizer or bounded search is used.
"""

if not __debug__:
    raise RuntimeError('Exact isoclinic probability checks require assertions; do not use -O.')

from fractions import Fraction as F
import json


def angle_bounds(cosine):
    """Return exact rational necessary bounds for a rational interior cosine."""
    u = F(cosine)
    if not F(1, 3) < u < F(2, 3):
        raise ValueError('a thirty-event isoclinic pair requires 1/3 < u < 2/3')
    t = u*(1-u)
    upper = t/(2+t)
    return {
        'cosine': str(u),
        'probability_lower_strict': str(1-9*upper),
        'probability_upper': str(upper),
        'upper_is_strict_by_corollary32': u == F(1, 2),
        'positive_trace_upper': str(F(9, 2)*(1-u)-1/u),
        'extra_space_trace_upper': str(9-2/t),
        'negative_trace_lower': str(2/t),
        'dimension19_positive_trace_equals_upper': True,
    }


def multiply(A, B):
    return [[sum((a*b for a, b in zip(row, column)), F(0))
             for column in zip(*B)] for row in A]


def add(A, B):
    return [[a+b for a, b in zip(row, other)] for row, other in zip(A, B)]


def scale(value, matrix):
    return [[value*x for x in row] for row in matrix]


def block_control():
    """Exact rational spectral blocks at u=3/5; no GHZ claim."""
    u, v = F(3, 5), F(4, 5)
    assert u*u+v*v == 1
    A = [[F(1), F(0)], [F(0), F(0)]]
    B = [[u*u, u*v], [u*v, v*v]]
    minus = scale(F(1, 2), [[1-u, -v], [-v, 1+u]])
    plus = scale(F(1, 2), [[1+u, v], [v, 1-u]])
    identity = [[F(1), F(0)], [F(0), F(1)]]
    zero = [[F(0), F(0)], [F(0), F(0)]]
    for projection in (A, B, minus, plus):
        assert multiply(projection, projection) == projection
    assert add(minus, plus) == identity and multiply(minus, plus) == zero
    T = scale(F(1, 2), add(multiply(A, B), multiply(B, A)))
    lower, upper = (u*u-u)/2, (u*u+u)/2
    assert T == add(scale(lower, minus), scale(upper, plus))
    assert multiply(T, minus) == scale(lower, minus)
    assert multiply(T, plus) == scale(upper, plus)
    # The negative line has squared overlap (1-u)/2 with both pair lines.
    assert sum(multiply(A, minus)[i][i] for i in range(2)) == (1-u)/2
    assert sum(multiply(B, minus)[i][i] for i in range(2)) == (1-u)/2
    q = (u*(1-u))/(2+u*(1-u))
    assert q == F(3, 28) and q+lower*(1-q) == 0
    # Adding a zero eigenspace does not change the spectral minimum.
    for exterior_mass in (F(0), F(1, 10), F(1, 2)):
        assert lower*(1-q-exterior_mass) >= lower*(1-q)
    return {'cosine': '3/5', 'local_spectral_cap': '3/28',
            'scope': 'two projections and one local spectral equality, not a GHZ realization'}


def self_test():
    for denominator in (9, 60, 120):
        for numerator in range(1, denominator):
            u = F(numerator, denominator)
            t = u*(1-u)
            q = t/(2+t)
            assert F(1, 9)-q == 2*(2*u-1)**2/(9*(2+t))
            assert 1-9*q == 2*(2*u-1)**2/(2+t)
            assert 9*t-2 == -(3*u-1)*(3*u-2)
            assert (10*q >= 1) == (F(1, 3) <= u <= F(2, 3))
            assert q-t*(1-q)/2 == 0
            s = F(9, 2)*(1-u)-1/u
            # Exact squared AM-GM identity; the left square root is positive.
            assert (F(9, 2)-s)**2-18 == (9*u*u-2)**2/(4*u*u)
            assert 9-2/t <= 1 and 2/t >= 8
            lower, upper = (u*u-u)/2, (u*u+u)/2
            for r in (F(0), F(1, 7), F(1)):
                b = s-(1-u)*r/2
                a = 9-b-r
                assert 1+lower*a+upper*b == 0
                assert a == 2/t+(1+u)*b/(1-u)
            if F(1, 3) < u < F(2, 3):
                assert s > 0
                assert 9-2/t > 0
                assert angle_bounds(u)['probability_upper'] == str(q)
    for u in (F(1, 3), F(2, 3)):
        assert u*(1-u)/(2+u*(1-u)) == F(1, 10)
        assert F(9, 2)*(1-u)-1/u == 0
        # Nine reduced event vectors at endpoint have positive Gram determinant.
        assert 1-9*F(1, 10) == F(1, 10)
    for u in (F(0), F(1, 3), F(2, 3), F(1)):
        try:
            angle_bounds(u)
        except ValueError:
            pass
        else:
            raise AssertionError('out-of-range cosine accepted')
    assert angle_bounds(F(1, 2))['probability_lower_strict'] == '0'
    assert angle_bounds(F(1, 2))['extra_space_trace_upper'] == '1'
    assert F(9, 2)**2 > 18  # 9/2-3*sqrt(2) is positive.
    return {'scope': 'one isoclinic reduced pair; arbitrary ambient dimension and probabilities',
            'probability_interval': '(1-9*t/(2+t), t/(2+t)] with t=u*(1-u)',
            'cosine_interval': '(1/3, 2/3)',
            'positive_trace_universal_upper': '9/2-3*sqrt(2)',
            'exterior_trace_universal_upper': '1',
            'negative_trace_universal_lower': '8',
            'rank_conclusion': False,
            'rational_block_control': block_control()}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cosine', help='optional exact rational interior cosine')
    args = parser.parse_args()
    print(json.dumps(angle_bounds(args.cosine) if args.cosine is not None else self_test(), indent=2))

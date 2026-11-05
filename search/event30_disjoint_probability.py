"""Exact probability bounds from pairwise disjoint cross-neighborhoods.

Let U contain m>=2 events from one context of any three-context GHZ graph,
and suppose their neighborhoods in another context are pairwise disjoint.
Write E,F for projections onto the two other context spans and choose phases
so <psi,g_i>=sqrt(c_i)>0. Set X_i=sqrt(c_i)g_i, P_i=E X_i, Q_i=F X_i.

The Gram matrix of psi and these3m vectors is PSD. We have XX diagonal c_i,
PP diagonal (disjoint neighborhoods), XP=PP, XQ=QQ, PQ diagonal zero (triangle
freeness), and <psi,X_i>=<psi,P_i>=<psi,Q_i>=c_i. Take its real part and
average simultaneous permutations of the m indices. This is a relaxation of
necessary Gram conditions, not a claim that a graph has those symmetries.

Let c be the mean of the c_i, p,q the mean P/Q diagonal entries, and r,t the
mean off-diagonal QQ/PQ entries. The standard permutation representation gives

    S = [[c,p,q-r], [p,p,-t], [q-r,-t,q-r]] >= 0.

The vectors psi, mean(X_i), mean(P_i), mean(Q_i) give

    T = [[1,c,c,c],
         [c,c/m,p/m,(q+(m-1)r)/m],
         [c,p/m,p/m,(m-1)t/m],
         [c,(q+(m-1)r)/m,(m-1)t/m,(q+(m-1)r)/m]] >= 0.

A rational PSD pair Z_S,Z_T with exact identity
u-c = Tr(Z_S S)+Tr(Z_T T) proves sum_{i in U}c_i <= m*u. Neither symmetry,
uniform event probabilities, nor a dimension bound is assumed. In particular
the saved m=5 certificate proves total probability <=15/32. Saved certificates
for m=6..10 fall below (m-1)/9, so those disjoint-neighborhood subsets cannot
occur among thirty positive events.

For an actual five-vertex automorphism orbit, all per-event probabilities can
be made equal by weighted-Gram averaging. Quotient column counts <=1 toward
another context guarantee disjoint neighborhoods and hence c_i<=3/32. The
second certificate format combines this upper bound with the rational sector
inequalities from event30_probability_lp. It verifies the moment certificate,
the quotient geometry, and the exact Farkas identity using only stdlib.
"""
if not __debug__:
    raise RuntimeError('Exact disjoint-probability verification requires Python assertions; do not use -O.')

from fractions import Fraction as F
import itertools
import json

from event30_probability_lp import constraint as sector_constraint, validate_quotient
from overnight30_screen import exact_psd


def moment_coefficients(m):
    assert type(m) is int and m >= 2
    standard, trivial = [], []
    for k in range(5):
        c, p, q, r, t = [F(int(i == k)) for i in range(5)]
        standard.append([[c, p, q-r], [p, p, -t], [q-r, -t, q-r]])
        trivial.append([[F(0), c, c, c],
                        [c, c/m, p/m, (q+(m-1)*r)/m],
                        [c, p/m, p/m, (m-1)*t/m],
                        [c, (q+(m-1)*r)/m, (m-1)*t/m, (q+(m-1)*r)/m]])
    return standard, trivial


def inner(left, right):
    return sum((a*b for row, other in zip(left, right) for a, b in zip(row, other)), F(0))


def verify_bound_certificate(certificate):
    assert certificate['format'] == 'ghz-disjoint-neighborhood-mean-bound-v1'
    m = certificate['subset_size']
    standard, trivial = moment_coefficients(m)
    upper = F(certificate['mean_probability_upper'])
    assert upper > 0
    A = [list(map(F, row)) for row in certificate['standard_dual']]
    B = [list(map(F, row)) for row in certificate['trivial_dual']]
    assert len(A) == 3 and all(len(row) == 3 for row in A)
    assert len(B) == 4 and all(len(row) == 4 for row in B)
    assert exact_psd(A) and exact_psd(B)
    assert B[0][0] == upper
    assert [inner(A, S)+inner(B, T) for S, T in zip(standard, trivial)] == [F(-1), F(0), F(0), F(0), F(0)]
    return upper


def bound_is_strict(certificate):
    """A positive-definite trivial dual makes the upper bound strict: T00=1."""
    verify_bound_certificate(certificate)
    matrix = [list(map(F, row)) for row in certificate['trivial_dual']]
    for i in range(4):
        pivot = matrix[i][i]
        if pivot <= 0:
            return False
        for j in range(i+1, 4):
            for k in range(j, 4):
                value = matrix[j][k]-matrix[i][j]*matrix[i][k]/pivot
                matrix[j][k] = matrix[k][j] = value
    return True


def disjoint_contexts(quotient):
    """Choose one certified disjoint-neighborhood context for each eligible orbit."""
    validate_quotient(quotient, [5]*6, [[0, 1], [2, 3], [4, 5]])
    result = []
    for i in range(6):
        choices = [k for k in range(3) if k != i//2
                   and all(quotient[j][i] <= 1 for j in range(2*k, 2*k+2))]
        if choices:
            result.append([i, choices[0]])
    return result


def verify_quotient_certificate(certificate):
    assert certificate['format'] == 'ghz30-disjoint-orbit-probability-farkas-v1'
    quotient, sizes, parts = (certificate[k] for k in ('quotient', 'sizes', 'parts'))
    validate_quotient(quotient, sizes, parts)
    assert sizes == [5]*6 and parts == [[0, 1], [2, 3], [4, 5]]
    moment = certificate['mean_bound_certificate']
    upper = verify_bound_certificate(moment)
    assert moment['subset_size'] == 5 and upper <= F(1, 9)
    eligible = {}
    for i, context in certificate['disjoint_contexts']:
        assert type(i) is int and 0 <= i < 6 and i not in eligible
        assert type(context) is int and 0 <= context < 3 and context != i//2
        assert all(quotient[j][i] <= 1 for j in parts[context])
        eligible[i] = context
    total, rhs, mass = [F(0)]*6, F(0), F(0)
    for term in certificate['inequalities']:
        multiplier = F(term['multiplier'])
        assert multiplier >= 0
        descriptor = term['constraint']
        row, bound = sector_constraint(quotient, sizes, parts, descriptor)
        if descriptor['kind'] == 'upper' and descriptor['orbit'] in eligible:
            bound = upper
        total = [a+multiplier*b for a, b in zip(total, row)]
        rhs += multiplier*bound
        mass += multiplier
    assert mass == 1
    equalities = list(map(F, certificate['equalities']))
    assert len(equalities) == 3
    for part, multiplier in zip(parts, equalities):
        for i in part:
            total[i] += multiplier*sizes[i]
        rhs += multiplier
    assert all(value == 0 for value in total) and rhs < 0
    assert str(-rhs) == certificate['contradiction_margin']
    return -rhs


def certify_quotient(quotient, mean_bound_certificate, denominator=8, refinements=8):
    """Discover exact multipliers; a missing certificate is inconclusive."""
    import numpy as np
    from scipy.optimize import linprog
    from sympy import Matrix, Rational
    upper = verify_bound_certificate(mean_bound_certificate)
    assert mean_bound_certificate['subset_size'] == 5 and upper <= F(1, 9)
    sizes, parts = [5]*6, [[0, 1], [2, 3], [4, 5]]
    contexts = disjoint_contexts(quotient)
    eligible = dict(contexts)
    descriptors = [{'kind': kind, 'orbit': i} for i in range(6) for kind in ('lower', 'upper')]
    for i in range(6):
        for a in range(denominator+1):
            for b in range(denominator-a+1):
                descriptors.append({'kind': 'sector', 'orbit': i,
                                    'weights': [str(F(a, denominator)), str(F(b, denominator)),
                                                str(F(denominator-a-b, denominator))]})
    E = [[F(5) if i in part else F(0) for i in range(6)] for part in parts]
    for _ in range(refinements+1):
        rows = []
        for item in descriptors:
            row, bound = sector_constraint(quotient, sizes, parts, item)
            if item['kind'] == 'upper' and item['orbit'] in eligible:
                bound = upper
            rows.append((row, bound))
        A, bounds = zip(*rows)
        result = linprog([0.0]*6+[1.0],
                         A_ub=np.array([[float(a) for a in row]+[-1.0] for row in A]),
                         b_ub=np.array(list(map(float, bounds))),
                         A_eq=np.array([[float(a) for a in row]+[0.0] for row in E]),
                         b_eq=np.ones(3), bounds=[(None, None)]*7, method='highs')
        if not result.success or not np.all(np.isfinite(result.x)):
            return None
        if result.fun > 1e-9:
            active = [i for i, value in enumerate(result.ineqlin.marginals) if value < -1e-9]
            columns = [list(A[i])+[F(1)] for i in active]+[row+[F(0)] for row in E]
            matrix = Matrix([[Rational(v.numerator, v.denominator) for v in column] for column in columns]).T
            try:
                solution, parameters = matrix.gauss_jordan_solve(Matrix([0]*6+[1]))
                solution = solution.subs({symbol: 0 for symbol in parameters})
                multipliers = [F(int(v.p), int(v.q)) for v in solution]
            except (ValueError, AttributeError, ZeroDivisionError):
                return None
            y, z = multipliers[:len(active)], multipliers[len(active):]
            margin = -sum((v*bounds[i] for i, v in zip(active, y)), F(0))-sum(z)
            if any(value < 0 for value in y) or margin <= 0:
                return None
            certificate = {'format': 'ghz30-disjoint-orbit-probability-farkas-v1',
                           'quotient': quotient, 'sizes': sizes, 'parts': parts,
                           'mean_bound_certificate': mean_bound_certificate, 'disjoint_contexts': contexts,
                           'inequalities': [{'constraint': descriptors[i], 'multiplier': str(v)}
                                            for i, v in zip(active, y) if v],
                           'equalities': list(map(str, z)), 'contradiction_margin': str(margin)}
            verify_quotient_certificate(certificate)
            return certificate
        added = False
        for i in range(6):
            others = [part for k, part in enumerate(parts) if k != i//2]
            alpha, beta = [sum(quotient[i][j]*result.x[j] for j in part) for part in others]
            rho = 1-alpha-beta
            if min(alpha, beta, rho) <= 0:
                continue
            reciprocal = [1/alpha, 1/beta, 1/rho]
            total = sum(reciprocal)
            if result.x[i] <= 1/total+1e-9:
                continue
            s, t = [F(value/total).limit_denominator(1000) for value in reciprocal[:2]]
            item = {'kind': 'sector', 'orbit': i, 'weights': list(map(str, (s, t, 1-s-t)))}
            if item not in descriptors:
                descriptors.append(item)
                added = True
        if not added:
            return None
    return None


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path)
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text())
    if certificate['format'] == 'ghz-disjoint-neighborhood-mean-bound-v1':
        value = verify_bound_certificate(certificate)
        result = {'subset_size': certificate['subset_size'], 'mean_upper': str(value),
                  'total_upper': str(certificate['subset_size']*value)}
    else:
        result = {'exact_contradiction_margin': str(verify_quotient_certificate(certificate))}
    print(json.dumps(result, indent=2))

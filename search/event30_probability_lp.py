"""Exact conditional exclusions of automorphism-orbit quotients on 30 events.

The input is an equitable quotient Q of an actual context-preserving group,
not merely an equitable partition. Averaging a state-weighted Gram matrix
gives per-event orbit probabilities x_i, with sum m_i*x_i=1 in each context.
For the two other-context neighbor masses alpha,beta and rho=1-alpha-beta,
Lemma 14 implies x_i <= F(alpha,beta,rho). For any rational s+t+u=1,

    F(alpha,beta,rho) <= s*s*alpha+t*t*beta+u*u*rho.

These valid linear inequalities, 0<=x_i<=1/9, and the context equations
yield rational Farkas certificates. Floating-point linear programming only
discovers multipliers; verification uses standard-library exact arithmetic.
An unverified or absent certificate proves nothing. A quotient certificate
does not itself certify that an input graph has the requisite automorphisms.
"""
from fractions import Fraction as F
import itertools
import json


def validate_quotient(quotient, sizes, parts):
    n = len(sizes)
    assert len(parts) == 3 and sorted(sum(parts, [])) == list(range(n))
    assert all(type(m) is int and 1 <= m <= 10 for m in sizes)
    assert all(sum(sizes[i] for i in part) == 10 for part in parts)
    assert len(quotient) == n and all(len(row) == n for row in quotient)
    color = {i: k for k, part in enumerate(parts) for i in part}
    for i, j in itertools.product(range(n), repeat=2):
        q = quotient[i][j]
        assert type(q) is int and 0 <= q <= sizes[j]
        assert color[i] != color[j] or q == 0
        assert sizes[i]*q == sizes[j]*quotient[j][i]
    return color


def constraint(quotient, sizes, parts, descriptor):
    n, row = len(sizes), [F(0)]*len(sizes)
    kind, i = descriptor['kind'], descriptor['orbit']
    assert type(i) is int and 0 <= i < n
    if kind == 'upper':
        row[i] = F(1)
        return row, F(1, 9)
    if kind == 'lower':
        row[i] = F(-1)
        # From the <1/9 bound on the 10-m events in the same context.
        return row, -F(sizes[i]-1, 9*sizes[i])
    assert kind == 'sector'
    s, t, u = map(F, descriptor['weights'])
    assert s+t+u == 1
    own = next(k for k, part in enumerate(parts) if i in part)
    others = [part for k, part in enumerate(parts) if k != own]
    row[i] = F(1)
    for j in others[0]:
        row[j] -= (s*s-u*u)*quotient[i][j]
    for j in others[1]:
        row[j] -= (t*t-u*u)*quotient[i][j]
    return row, u*u


def verify_certificate(certificate):
    assert certificate['format'] == 'ghz30-orbit-probability-farkas-v1'
    quotient, sizes, parts = (certificate[k] for k in ('quotient', 'sizes', 'parts'))
    validate_quotient(quotient, sizes, parts)
    total, rhs, mass = [F(0)]*len(sizes), F(0), F(0)
    for term in certificate['inequalities']:
        multiplier = F(term['multiplier'])
        assert multiplier >= 0
        row, bound = constraint(quotient, sizes, parts, term['constraint'])
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


def certify(quotient, sizes=None, parts=None, denominator=8, refinements=5):
    """Return a verified certificate, or None; uses scipy/sympy for discovery."""
    import numpy as np
    from scipy.optimize import linprog
    from sympy import Matrix, Rational
    sizes = sizes or [5]*6
    parts = parts or [[0, 1], [2, 3], [4, 5]]
    color = validate_quotient(quotient, sizes, parts)
    n = len(sizes)
    descriptors = [{'kind': kind, 'orbit': i} for i in range(n) for kind in ('lower', 'upper')]
    for i in range(n):
        for a in range(denominator+1):
            for b in range(denominator-a+1):
                descriptors.append({'kind': 'sector', 'orbit': i,
                                    'weights': [str(F(a, denominator)), str(F(b, denominator)),
                                                str(F(denominator-a-b, denominator))]})
    E = [[F(sizes[i]) if i in part else F(0) for i in range(n)] for part in parts]
    for iteration in range(refinements+1):
        rows = [constraint(quotient, sizes, parts, item) for item in descriptors]
        A, bounds = zip(*rows)
        result = linprog([0.0]*n+[1.0],
                         A_ub=np.array([[float(a) for a in row]+[-1.0] for row in A]),
                         b_ub=np.array(list(map(float, bounds))),
                         A_eq=np.array([[float(a) for a in row]+[0.0] for row in E]),
                         b_eq=np.ones(3), bounds=[(None, None)]*(n+1), method='highs')
        if not result.success or not np.all(np.isfinite(result.x)):
            return None
        if result.fun > 1e-9:
            active = [i for i, value in enumerate(result.ineqlin.marginals) if value < -1e-9]
            columns = [list(A[i])+[F(1)] for i in active]+[row+[F(0)] for row in E]
            matrix = Matrix([[Rational(v.numerator, v.denominator) for v in column]
                             for column in columns]).T
            try:
                solution, parameters = matrix.gauss_jordan_solve(Matrix([0]*n+[1]))
                solution = solution.subs({symbol: 0 for symbol in parameters})
                multipliers = [F(int(v.p), int(v.q)) for v in solution]
            except (ValueError, AttributeError, ZeroDivisionError):
                return None
            y, z = multipliers[:len(active)], multipliers[len(active):]
            margin = -sum((v*bounds[i] for i, v in zip(active, y)), F(0))-sum(z)
            if any(value < 0 for value in y) or margin <= 0:
                return None
            certificate = {'format': 'ghz30-orbit-probability-farkas-v1',
                           'quotient': quotient, 'sizes': sizes, 'parts': parts,
                           'inequalities': [{'constraint': descriptors[i], 'multiplier': str(v)}
                                            for i, v in zip(active, y) if v],
                           'equalities': list(map(str, z)), 'contradiction_margin': str(margin)}
            verify_certificate(certificate)
            return certificate
        added = False
        for i in range(n):
            others = [part for k, part in enumerate(parts) if k != color[i]]
            alpha, beta = [sum(quotient[i][j]*result.x[j] for j in part) for part in others]
            rho = 1-alpha-beta
            if min(alpha, beta, rho) <= 0:
                continue
            reciprocal = [1/alpha, 1/beta, 1/rho]
            total = sum(reciprocal)
            if result.x[i] <= 1/total+1e-9:
                continue
            s, t = [F(value/total).limit_denominator(1000) for value in reciprocal[:2]]
            descriptor = {'kind': 'sector', 'orbit': i, 'weights': list(map(str, (s, t, 1-s-t)))}
            if descriptor not in descriptors:
                descriptors.append(descriptor)
                added = True
        if not added:
            return None
    return None


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path, nargs='?')
    args = parser.parse_args()
    if args.certificate:
        print('Verified exact contradiction margin:', verify_certificate(json.loads(args.certificate.read_text())))
    else:
        quotient = [[int(i//2 != j//2) for j in range(6)] for i in range(6)]
        certificate = certify(quotient)
        assert certificate is not None
        print(json.dumps(certificate, indent=2))

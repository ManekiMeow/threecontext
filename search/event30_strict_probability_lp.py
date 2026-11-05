"""Exact strict Farkas certificates for six actual five-event symmetry orbits.

The original sector inequalities remain nonstrict. On thirty positive events,
each per-event orbit probability lies strictly between 4/45 and 1/9. Theorem56
also gives strict weighted degree between two and five for each ordered
context pair. Verified disjoint-neighborhood moment bounds may strengthen an
eligible orbit's upper bound; their strictness is checked, not assumed.

For inequalities a_i*x<=b_i (some strict) and context equations E*x=1, a
certificate supplies nonnegative y and unrestricted z with y*A+z*E=0,
sum(y_i over strict rows)=1, and y*b+sum(z)<=0. Combining the strict rows then
gives 0<y*b+sum(z)<=0. Thus exact zero-margin certificates are sound.

Actual automorphism orbits, rather than merely equitable counts, are required
to average probabilities. Averaging can increase rank; these are event-number
obstructions in every dimension, including complex realizations.
"""
if not __debug__:
    raise RuntimeError('Exact strict Farkas verification requires assertions; do not use -O.')

from fractions import Fraction as F
import itertools
import json

from event30_probability_lp import constraint as old_constraint, validate_quotient
from event30_disjoint_probability import verify_bound_certificate, bound_is_strict, disjoint_contexts

SIZES = [5]*6
PARTS = [[0, 1], [2, 3], [4, 5]]


def geometry(certificate):
    quotient, sizes, parts = (certificate[k] for k in ('quotient', 'sizes', 'parts'))
    validate_quotient(quotient, sizes, parts)
    assert sizes == SIZES and parts == PARTS
    moment = certificate.get('mean_bound_certificate')
    if moment is None:
        assert not certificate.get('disjoint_contexts')
        return quotient, {}, None, False
    upper = verify_bound_certificate(moment)
    assert moment['subset_size'] == 5 and upper <= F(1, 9)
    strict = bound_is_strict(moment)
    eligible = {}
    for i, other in certificate['disjoint_contexts']:
        assert type(i) is int and 0 <= i < 6 and i not in eligible
        assert type(other) is int and 0 <= other < 3 and other != i//2
        assert all(quotient[j][i] <= 1 for j in PARTS[other])
        eligible[i] = other
    return quotient, eligible, upper, strict


def constraint(quotient, descriptor, eligible, moment_upper, moment_strict):
    kind = descriptor['kind']
    if kind in ('weighted_degree_lower', 'weighted_degree_upper'):
        a, b = descriptor['event_context'], descriptor['source_context']
        assert type(a) is int and type(b) is int and 0 <= a < 3 and 0 <= b < 3 and a != b
        row = [F(0)]*6
        for j in PARTS[b]:
            row[j] = SIZES[j]*sum(quotient[j][i] for i in PARTS[a])
        if kind == 'weighted_degree_lower':
            return [-value for value in row], F(-2), True
        return row, F(5), True
    row, upper = old_constraint(quotient, SIZES, PARTS, descriptor)
    strict = kind in ('lower', 'upper')
    if kind == 'upper' and descriptor['orbit'] in eligible:
        upper, strict = moment_upper, moment_strict
    return row, upper, strict


def verify_certificate(certificate):
    assert certificate['format'] == 'ghz30-strict-orbit-probability-farkas-v1'
    quotient, eligible, moment_upper, moment_strict = geometry(certificate)
    total, rhs, strict_mass = [F(0)]*6, F(0), F(0)
    for term in certificate['inequalities']:
        multiplier = F(term['multiplier'])
        assert multiplier >= 0
        row, bound, strict = constraint(quotient, term['constraint'], eligible, moment_upper, moment_strict)
        total = [a+multiplier*b for a, b in zip(total, row)]
        rhs += multiplier*bound
        if strict:
            strict_mass += multiplier
    assert strict_mass == 1
    equalities = list(map(F, certificate['equalities']))
    assert len(equalities) == 3
    for part, multiplier in zip(PARTS, equalities):
        for i in part:
            total[i] += multiplier*SIZES[i]
        rhs += multiplier
    assert all(value == 0 for value in total) and rhs <= 0
    assert str(-rhs) == certificate['contradiction_margin']
    return -rhs


def certify(quotient, mean_bound_certificate=None, denominator=8, refinements=5):
    """Numerically discover multipliers, then require an exact certificate."""
    import numpy as np
    from scipy.optimize import linprog
    from sympy import Matrix, Rational
    header = {'format': 'ghz30-strict-orbit-probability-farkas-v1',
              'quotient': quotient, 'sizes': SIZES, 'parts': PARTS}
    if mean_bound_certificate is not None:
        header.update(mean_bound_certificate=mean_bound_certificate,
                      disjoint_contexts=disjoint_contexts(quotient))
    quotient, eligible, moment_upper, moment_strict = geometry(header)
    descriptors = [{'kind': kind, 'orbit': i} for i in range(6) for kind in ('lower', 'upper')]
    descriptors += [{'kind': kind, 'event_context': a, 'source_context': b}
                    for a, b in itertools.permutations(range(3), 2)
                    for kind in ('weighted_degree_lower', 'weighted_degree_upper')]
    for i in range(6):
        for a in range(denominator+1):
            for b in range(denominator-a+1):
                descriptors.append({'kind': 'sector', 'orbit': i,
                                    'weights': list(map(str, (F(a, denominator), F(b, denominator),
                                                             F(denominator-a-b, denominator))))})
    E = [[F(5) if i in part else F(0) for i in range(6)] for part in PARTS]
    for _ in range(refinements+1):
        rows = [constraint(quotient, item, eligible, moment_upper, moment_strict) for item in descriptors]
        A, bounds, strict = zip(*rows)
        result = linprog([0.0]*6+[-1.0],
                         A_ub=np.array([[float(a) for a in row]+[float(flag)] for row, flag in zip(A, strict)]),
                         b_ub=np.array(list(map(float, bounds))),
                         A_eq=np.array([[float(a) for a in row]+[0.0] for row in E]),
                         b_eq=np.ones(3), bounds=[(None, None)]*7, method='highs')
        if not result.success or not np.all(np.isfinite(result.x)):
            return None
        if result.x[-1] <= 1e-9:
            active = [i for i, value in enumerate(result.ineqlin.marginals) if value < -1e-9]
            columns = [list(A[i])+[F(strict[i])] for i in active]+[row+[F(0)] for row in E]
            matrix = Matrix([[Rational(v.numerator, v.denominator) for v in column] for column in columns]).T
            try:
                solution, parameters = matrix.gauss_jordan_solve(Matrix([0]*6+[1]))
                solution = solution.subs({symbol: 0 for symbol in parameters})
                multipliers = [F(int(value.p), int(value.q)) for value in solution]
            except (ValueError, AttributeError, ZeroDivisionError):
                multipliers = None
            if multipliers is not None:
                y, z = multipliers[:len(active)], multipliers[len(active):]
                rhs = sum((v*bounds[i] for i, v in zip(active, y)), F(0))+sum(z)
                if all(value >= 0 for value in y) and rhs <= 0:
                    certificate = dict(header,
                        inequalities=[{'constraint': descriptors[i], 'multiplier': str(v)}
                                      for i, v in zip(active, y) if v],
                        equalities=list(map(str, z)), contradiction_margin=str(-rhs))
                    verify_certificate(certificate)
                    return certificate
        added = False
        for i in range(6):
            others = [part for k, part in enumerate(PARTS) if k != i//2]
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


def self_test(discovery=False):
    """Boundary contradictions and false-certificate controls, without an SDP."""
    from copy import deepcopy
    controls = []
    for degree in (2, 5):
        quotient = [[0 if i//2 == j//2 else (1 if degree == 2 else 2+int(i % 2 != j % 2))
                     for j in range(6)] for i in range(6)]
        lower = degree == 2
        proof = {'format': 'ghz30-strict-orbit-probability-farkas-v1',
                 'quotient': quotient, 'sizes': SIZES, 'parts': PARTS,
                 'inequalities': [{'constraint': {
                     'kind': 'weighted_degree_lower' if lower else 'weighted_degree_upper',
                     'event_context': 0, 'source_context': 1}, 'multiplier': '1'}],
                 'equalities': ['0', '2' if lower else '-5', '0'],
                 'contradiction_margin': '0'}
        assert verify_certificate(proof) == 0
        if discovery:
            found = certify(quotient)
            assert found is not None and verify_certificate(found) >= 0
        controls.append(proof)
    invalid = []
    bad = deepcopy(controls[0]); bad['inequalities'][0]['multiplier'] = '-1'; invalid.append(bad)
    bad = deepcopy(controls[0]); bad['inequalities'][0]['multiplier'] = '2'; invalid.append(bad)
    bad = deepcopy(controls[0]); bad['equalities'][1] = '1'; invalid.append(bad)
    bad = deepcopy(controls[0]); bad['contradiction_margin'] = '1/1000'; invalid.append(bad)
    bad = deepcopy(controls[0]); bad.update(inequalities=[], equalities=['0']*3); invalid.append(bad)
    # This degree-three identity has a positive right side and proves nothing.
    bad = deepcopy(controls[0])
    bad['quotient'] = [[0 if i//2 == j//2 else 1+int(i % 2 != j % 2)
                        for j in range(6)] for i in range(6)]
    bad['equalities'][1] = '3'; bad['contradiction_margin'] = '-1'; invalid.append(bad)
    for bad in invalid:
        try:
            verify_certificate(bad)
        except AssertionError:
            pass
        else:
            raise AssertionError('Invalid strict contradiction was accepted')
    return {'exact_zero_margin_controls': len(controls), 'false_certificates_rejected': len(invalid),
            'numerical_discovery_controls': discovery}


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path, nargs='?')
    parser.add_argument('--discover-controls', action='store_true')
    args = parser.parse_args()
    if args.certificate is not None:
        print('Verified exact strict contradiction margin:', verify_certificate(json.loads(args.certificate.read_text())))
    else:
        print(json.dumps(self_test(args.discover_controls), indent=2))

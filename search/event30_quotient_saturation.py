"""Exact affine and Cauchy-equality propagation for six five-event orbits.

The quotient must come from an actual context-preserving automorphism group;
an equitable partition alone is insufficient. Averaging its state-weighted
Gram matrix yields orbit probabilities x_i, with x_0+x_1=x_2+x_3=x_4+x_5=1/5.
Rank may increase, so this is an event-number obstruction in every dimension.

The rules below use exact rational arithmetic:

* Constant neighbor mass alpha gives x_i <= h(alpha), where
  h(alpha)=alpha(1-alpha)/(1+3alpha). This bound is nonstrict; only the global
  x_i<1/9 bound is initially strict.
* If the two orbit upper bounds in a context sum to1/5, both are attained.
  A smaller sum, or equality involving a strict upper bound, is impossible.
* When x_i=h(alpha) with both quantities fixed by affine equations,
  Cauchy equality forces the other neighbor mass to (1-alpha)/2.
* Fixed x_i,alpha,beta saturating x_i=F(alpha,beta,1-alpha-beta) force the
  overlap on every incident edge to sqrt(x_i*x_j)/alpha, using the neighbor
  mass toward that edge's other context. Two saturated endpoints of an edge
  must therefore have equal such masses. A mismatch is a contradiction.

Every affine consequence is checked by rational RREF. Certificates contain a
deterministic derivation transcript; verification replays the proof without
any LP, SDP, numerical tolerance or external dependency. Failure to derive a
contradiction is inconclusive. This checker does not itself verify a graph's
automorphism group or membership in the conditional quotient class.
"""
from fractions import Fraction as F
import itertools
import json

from event30_probability_lp import validate_quotient

SIZES = [5]*6
PARTS = [[0, 1], [2, 3], [4, 5]]


class AffineSystem:
    """Small exact RREF system A*x=b, with deterministic pivot ordering."""
    def __init__(self):
        self.rows = []

    def add(self, coefficients, rhs):
        row = list(map(F, coefficients))+[F(rhs)]
        assert len(row) == 7
        for pivot, basis in self.rows:
            scalar = row[pivot]
            row = [a-scalar*b for a, b in zip(row, basis)]
        nonzero = [i for i, value in enumerate(row[:6]) if value]
        if not nonzero:
            return 'inconsistent' if row[-1] else 'redundant'
        pivot = nonzero[0]
        scalar = row[pivot]
        row = [value/scalar for value in row]
        self.rows = [(p, [a-basis[pivot]*b for a, b in zip(basis, row)])
                     for p, basis in self.rows]
        self.rows.append((pivot, row))
        self.rows.sort()
        return 'added'

    def constant(self, coefficients):
        row, value = list(map(F, coefficients)), F(0)
        assert len(row) == 6
        for pivot, basis in self.rows:
            scalar = row[pivot]
            row = [a-scalar*b for a, b in zip(row, basis[:6])]
            value += scalar*basis[-1]
        return None if any(row) else value


def h(alpha):
    alpha = F(alpha)
    assert 0 < alpha < 1
    return alpha*(1-alpha)/(1+3*alpha)


def harmonic(alpha, beta):
    alpha, beta = F(alpha), F(beta)
    rho = 1-alpha-beta
    assert min(alpha, beta, rho) > 0
    return 1/(1/alpha+1/beta+1/rho)


def derive(quotient):
    """Return a deterministic exact exclusion transcript, or an inconclusion."""
    color = validate_quotient(quotient, SIZES, PARTS)
    assert len(quotient) == 6
    system = AffineSystem()
    units = [[F(int(i == j)) for j in range(6)] for i in range(6)]
    masses = {(i, k): [F(quotient[i][j]) if j in part else F(0) for j in range(6)]
              for i in range(6) for k, part in enumerate(PARTS) if k != color[i]}
    for part in PARTS:
        assert system.add([int(i in part) for i in range(6)], F(1, 5)) == 'added'
    upper, strict, trace = [F(1, 9)]*6, [True]*6, []

    def finish(terminal):
        return {'status': 'excluded', 'certificate': {
            'format': 'ghz30-quotient-saturation-v1', 'quotient': quotient,
            'sizes': SIZES, 'parts': PARTS, 'trace': trace, 'terminal': terminal}}

    def equation(coefficients, rhs, reason):
        status = system.add(coefficients, rhs)
        if status == 'added':
            trace.append({'rule': reason, 'coefficients': list(map(str, coefficients)), 'rhs': str(rhs)})
        return status

    for _ in range(12):
        changed = False
        for i in range(6):
            probability = system.constant(units[i])
            if probability is not None:
                if probability <= 0 or probability >= F(1, 9):
                    return finish({'kind': 'constant_probability_outside_range', 'orbit': i,
                                   'probability': str(probability)})
                if probability > upper[i] or (probability == upper[i] and strict[i]):
                    return finish({'kind': 'constant_probability_exceeds_upper', 'orbit': i,
                                   'probability': str(probability), 'upper': str(upper[i]),
                                   'strict_upper': strict[i]})
                if probability < upper[i]:
                    upper[i], strict[i], changed = probability, False, True
                    trace.append({'rule': 'constant_probability_bound', 'orbit': i,
                                  'upper': str(probability), 'strict': False})
            for k in range(3):
                if k == color[i]:
                    continue
                alpha = system.constant(masses[i, k])
                if alpha is None:
                    continue
                if not 0 < alpha < 1:
                    return finish({'kind': 'constant_neighbor_mass_outside_range', 'orbit': i,
                                   'neighbor_context': k, 'mass': str(alpha)})
                bound = h(alpha)
                if bound < upper[i]:
                    upper[i], strict[i], changed = bound, False, True
                    trace.append({'rule': 'constant_neighbor_bound', 'orbit': i,
                                  'neighbor_context': k, 'mass': str(alpha),
                                  'upper': str(bound), 'strict': False})
                if probability is not None and probability > bound:
                    return finish({'kind': 'constant_probability_exceeds_neighbor_bound', 'orbit': i,
                                   'probability': str(probability), 'neighbor_context': k,
                                   'mass': str(alpha), 'upper': str(bound)})
                if probability is not None and probability == bound:
                    other = 3-color[i]-k
                    rhs = (1-alpha)/2
                    reason = {'kind': 'one_sector_cauchy_equality', 'orbit': i,
                              'neighbor_context': k, 'mass': str(alpha), 'probability': str(probability)}
                    status = equation(masses[i, other], rhs, reason)
                    if status == 'inconsistent':
                        return finish({'kind': 'inconsistent_cauchy_equation', 'reason': reason,
                                       'coefficients': list(map(str, masses[i, other])), 'rhs': str(rhs)})
                    changed |= status == 'added'
        for k, part in enumerate(PARTS):
            total = sum((upper[i] for i in part), F(0))
            if total < F(1, 5) or (total == F(1, 5) and any(strict[i] for i in part)):
                return finish({'kind': 'context_upper_budget', 'context': k,
                               'orbit_upper_sum': str(total), 'strict_bounds': [strict[i] for i in part]})
            if total == F(1, 5):
                for i in part:
                    reason = {'kind': 'context_budget_equality', 'context': k,
                              'orbit_upper_bounds': [str(upper[j]) for j in part],
                              'strict_bounds': [strict[j] for j in part]}
                    status = equation(units[i], upper[i], reason)
                    if status == 'inconsistent':
                        return finish({'kind': 'inconsistent_context_equality', 'reason': reason,
                                       'orbit': i, 'probability': str(upper[i])})
                    changed |= status == 'added'
        saturated = {}
        for i in range(6):
            other_contexts = [k for k in range(3) if k != color[i]]
            alpha, beta = [system.constant(masses[i, k]) for k in other_contexts]
            probability = system.constant(units[i])
            if alpha is None or beta is None:
                continue
            if min(alpha, beta, 1-alpha-beta) <= 0:
                return finish({'kind': 'constant_sector_outside_simplex', 'orbit': i,
                               'masses': [str(alpha), str(beta), str(1-alpha-beta)]})
            if probability is None:
                continue
            bound = harmonic(alpha, beta)
            if probability > bound:
                return finish({'kind': 'constant_probability_exceeds_three_sector_bound', 'orbit': i,
                               'probability': str(probability), 'masses': [str(alpha), str(beta)],
                               'upper': str(bound)})
            if probability == bound:
                saturated[i] = {'probability': probability, 'masses': dict(zip(other_contexts, (alpha, beta)))}
        for i, j in itertools.combinations(sorted(saturated), 2):
            if quotient[i][j] == 0:
                continue
            a = saturated[i]['masses'][color[j]]
            b = saturated[j]['masses'][color[i]]
            if a != b:
                return finish({'kind': 'saturated_edge_mismatch', 'orbits': [i, j],
                               'probabilities': [str(saturated[v]['probability']) for v in (i, j)],
                               'neighbor_masses': [str(a), str(b)], 'quotient_edge_count': quotient[i][j]})
        if not changed:
            return {'status': 'not_excluded', 'reason': 'no further exact affine consequence',
                    'trace': trace, 'upper_bounds': list(map(str, upper)),
                    'strict_bounds': strict, 'saturated_orbits': sorted(saturated)}
    return {'status': 'not_excluded', 'reason': 'bounded derivation limit', 'trace': trace}


def verify_certificate(certificate):
    assert certificate['format'] == 'ghz30-quotient-saturation-v1'
    assert certificate['sizes'] == SIZES and certificate['parts'] == PARTS
    actual = derive(certificate['quotient'])
    assert actual['status'] == 'excluded'
    assert actual['certificate'] == certificate
    return certificate['terminal']


def self_test():
    from event30_c5_quotient import Q
    result = derive(Q)
    assert result['status'] == 'excluded'
    assert result['certificate']['terminal']['kind'] == 'saturated_edge_mismatch'
    verify_certificate(result['certificate'])
    regular = [[0 if i//2 == j//2 else 1+(i % 2 != j % 2) for j in range(6)] for i in range(6)]
    assert derive(regular)['status'] == 'not_excluded'
    # Uniform probabilities1/10 strictly satisfy the sector inequalities for
    # regular degree triples333 and334. These are relaxation controls only.
    for degrees in ((3, 3, 3), (3, 3, 4)):
        for entries in itertools.product(*(range(degree+1) for degree in degrees)):
            quotient = [[0]*6 for _ in range(6)]
            for (a, b), degree, diagonal in zip(((0, 1), (0, 2), (1, 2)), degrees, entries):
                block = [[diagonal, degree-diagonal], [degree-diagonal, diagonal]]
                for i, j in itertools.product(range(2), repeat=2):
                    quotient[2*a+i][2*b+j] = quotient[2*b+j][2*a+i] = block[i][j]
            assert derive(quotient)['status'] == 'not_excluded'
    for contexts in itertools.permutations(range(3)):
        for flips in itertools.product(range(2), repeat=3):
            order = [2*contexts[k]+(i ^ flips[k]) for k in range(3) for i in range(2)]
            quotient = [[Q[order[i]][order[j]] for j in range(6)] for i in range(6)]
            assert derive(quotient)['status'] == 'excluded'
    system = AffineSystem()
    assert system.add([1, 1, 0, 0, 0, 0], F(1, 5)) == 'added'
    assert system.constant([2, 2, 0, 0, 0, 0]) == F(2, 5)
    assert system.constant([1, 0, 0, 0, 0, 0]) is None
    assert system.add([1, 0, 0, 0, 0, 0], F(1, 10)) == 'added'
    assert system.constant([0, 1, 0, 0, 0, 0]) == F(1, 10)
    assert system.add([0, 1, 0, 0, 0, 0], F(1, 9)) == 'inconsistent'
    import copy
    for field in ('trace', 'terminal'):
        bad = copy.deepcopy(result['certificate'])
        if field == 'trace':
            bad[field].pop()
        else:
            bad[field]['neighbor_masses'][0] = '1/3'
        try:
            verify_certificate(bad)
        except AssertionError:
            pass
        else:
            raise AssertionError('tampered saturation certificate accepted')
    return {'controls': 'exact RREF, 48 Proposition46 relabelings, 144 feasible relaxations and tamper rejection passed'}


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('self-test')
    verify = sub.add_parser('verify')
    verify.add_argument('certificate', type=Path)
    run = sub.add_parser('derive')
    run.add_argument('quotient', type=Path)
    args = parser.parse_args()
    if args.command == 'self-test':
        output = self_test()
    elif args.command == 'verify':
        output = verify_certificate(json.loads(args.certificate.read_text()))
    else:
        record = json.loads(args.quotient.read_text())
        output = derive(record.get('quotient', record) if isinstance(record, dict) else record)
    print(json.dumps(output, indent=2))

"""Exact per-vertex probability obstructions for unrestricted GHZ30 supports.

The trivial automorphism group has thirty singleton orbits, so the sector
inequalities of event30_probability_lp apply without symmetry or uniformity.
This separate certificate format also allows Theorem56's strict weighted
degree bounds and Theorem53's verified disjoint-neighborhood subset bounds.
Corollary57b also bounds probabilities in any closed degree-two subset by1/10,
even when that subset contains no whole context; this upper bound is nonstrict.
Singleton probabilities satisfy 0<c_v<1/9; the five-orbit lower bound4/45 is
NOT used. Context normalization is exact.

Each row is a*x<=b, possibly strict. Nonnegative multipliers y and unrestricted
context multipliers z certify y*A+z*E=0 and y*b+sum(z)<0, or equality with a
positive multiplier on a strict row. The multipliers y sum to one. All
verification uses exact standard-library arithmetic and checks the saved
graph. Numerical LP failure or absence of a certificate proves nothing.
"""
if not __debug__:
    raise RuntimeError('Exact vertex-probability verification requires assertions; do not use -O.')

from fractions import Fraction as F
import itertools
import json

from overnight30_screen import normalize_graph, base_validation
from event30_probability_lp import constraint as orbit_constraint
from event30_disjoint_probability import verify_bound_certificate, bound_is_strict


def geometry(graph, moment_certificates=(), closure=None):
    assert normalize_graph(graph) == graph
    adj, colors = base_validation(graph)
    quotient = [[int(j in adj[i]) for j in range(30)] for i in range(30)]
    bounds = {}
    for proof in moment_certificates:
        m = proof['subset_size']
        assert type(m) is int and 2 <= m <= 10 and m not in bounds
        bounds[m] = (m*verify_bound_certificate(proof), bound_is_strict(proof))
    closed = set()
    if closure is not None:
        vertices, choices = closure['vertices'], closure['choices']
        assert isinstance(vertices, list) and vertices and vertices == sorted(set(vertices))
        assert all(type(v) is int and 0 <= v < 30 for v in vertices)
        assert isinstance(choices, list) and len(choices) == len(vertices)
        assert all(isinstance(choice, list) and len(choice) == 2 for choice in choices)
        assert [choice[0] for choice in choices] == vertices
        closed = set(vertices)
        for v, k in choices:
            assert type(v) is int and type(k) is int and 0 <= k < 3 and k != colors[v]
            neighbors = adj[v] & set(graph['parts'][k])
            assert len(neighbors) <= 2 and neighbors <= closed
    return adj, colors, quotient, bounds, closed


def constraint(graph, data, descriptor):
    adj, colors, quotient, bounds, closed = data
    kind, row = descriptor['kind'], [F(0)]*30
    if kind == 'closed_upper':
        v = descriptor['vertex']
        assert type(v) is int and v in closed
        row[v] = F(1)
        return row, F(1, 10), False
    if kind in ('weighted_degree_lower', 'weighted_degree_upper'):
        a, b = descriptor['event_context'], descriptor['source_context']
        assert type(a) is int and type(b) is int and 0 <= a < 3 and 0 <= b < 3 and a != b
        for v in graph['parts'][b]:
            row[v] = F(len(adj[v] & set(graph['parts'][a])))
        if kind == 'weighted_degree_lower':
            return [-x for x in row], F(-2), True
        return row, F(5), True
    if kind == 'disjoint_subset':
        vertices, other = descriptor['vertices'], descriptor['other_context']
        assert isinstance(vertices, list) and vertices == sorted(set(vertices))
        assert 2 <= len(vertices) <= 10 and all(type(v) is int and 0 <= v < 30 for v in vertices)
        own = colors[vertices[0]]
        assert all(colors[v] == own for v in vertices)
        assert type(other) is int and 0 <= other < 3 and other != own
        seen, target = set(), set(graph['parts'][other])
        for v in vertices:
            neighbors = adj[v] & target
            assert not seen.intersection(neighbors)
            seen.update(neighbors)
            row[v] = F(1)
        upper, strict = bounds[len(vertices)]
        return row, upper, strict
    row, upper = orbit_constraint(quotient, [1]*30, graph['parts'], descriptor)
    return row, upper, kind in ('lower', 'upper')


def verify_certificate(proof):
    assert proof['format'] == 'ghz30-vertex-probability-farkas-v1'
    graph = proof['graph']
    data = geometry(graph, proof['mean_bound_certificates'], proof.get('degree_two_closed_subset'))
    total, rhs, mass, strict_mass, used = [F(0)]*30, F(0), F(0), F(0), set()
    for term in proof['inequalities']:
        multiplier = F(term['multiplier'])
        assert multiplier > 0
        descriptor = term['constraint']
        row, bound, strict = constraint(graph, data, descriptor)
        total = [a+multiplier*b for a, b in zip(total, row)]
        rhs += multiplier*bound
        mass += multiplier
        strict_mass += multiplier*strict
        if descriptor['kind'] == 'disjoint_subset':
            used.add(len(descriptor['vertices']))
    assert mass == 1 and used == set(data[3])
    equalities = list(map(F, proof['equalities']))
    assert len(equalities) == 3
    for part, multiplier in zip(graph['parts'], equalities):
        for v in part:
            total[v] += multiplier
        rhs += multiplier
    assert all(value == 0 for value in total)
    assert rhs < 0 or (rhs == 0 and strict_mass > 0)
    assert proof['contradiction_margin'] == str(-rhs)
    assert proof['strict_multiplier_mass'] == str(strict_mass)
    return -rhs


def descriptors_for(graph, data, denominator):
    adj, colors, _, bounds, closed = data
    descriptors = [{'kind': kind, 'orbit': v} for v in range(30) for kind in ('lower', 'upper')]
    descriptors += [{'kind': 'closed_upper', 'vertex': v} for v in sorted(closed)]
    descriptors += [{'kind': kind, 'event_context': a, 'source_context': b}
                    for a, b in itertools.permutations(range(3), 2)
                    for kind in ('weighted_degree_lower', 'weighted_degree_upper')]
    for v in range(30):
        for a in range(denominator+1):
            for b in range(denominator-a+1):
                descriptors.append({'kind': 'sector', 'orbit': v,
                                    'weights': list(map(str, (F(a, denominator), F(b, denominator),
                                                             F(denominator-a-b, denominator))))})
    for own, part in enumerate(graph['parts']):
        seen_subsets = set()
        for other in range(3):
            if own == other:
                continue
            target = set(graph['parts'][other])
            unions, eligible = [set() for _ in range(1024)], [False]*1024
            eligible[0] = True
            for mask in range(1, 1024):
                low = mask & -mask
                i, rest = low.bit_length()-1, mask ^ low
                neighbors = adj[part[i]] & target
                unions[mask] = unions[rest] | neighbors
                eligible[mask] = eligible[rest] and not unions[rest].intersection(neighbors)
                if eligible[mask] and mask.bit_count() in bounds and mask not in seen_subsets:
                    seen_subsets.add(mask)
                    descriptors.append({'kind': 'disjoint_subset',
                                        'vertices': [part[i] for i in range(10) if mask >> i & 1],
                                        'other_context': other})
    return descriptors


def certify(record, denominator=4, refinements=5, time_limit=3.0):
    """Discover a candidate LP identity and return only an exact verified proof."""
    import time
    import numpy as np
    from scipy.optimize import linprog
    from sympy import Matrix, Rational
    from event30_disjoint_partition import catalog
    from event30_degree_two_cover import maximal_closed_subset
    started = time.monotonic()
    graph = normalize_graph(record)
    _, moment_proofs = catalog()
    vertices, choices = maximal_closed_subset(graph)
    closure = {'vertices': vertices, 'choices': choices} if vertices else None
    data = geometry(graph, moment_proofs.values(), closure)
    descriptors = descriptors_for(graph, data, denominator)
    E = [[F(int(v in part)) for v in range(30)] for part in graph['parts']]
    for iteration in range(refinements+1):
        if time.monotonic()-started >= time_limit:
            return None
        rows = [constraint(graph, data, item) for item in descriptors]
        A, bounds, strict = zip(*rows)
        # First seek an ordinary negative-margin Farkas identity. If it is
        # inconclusive, a second LP also detects strict boundary contradictions.
        for strict_mode in (False, True):
            last_column = [float(flag) if strict_mode else -1.0 for flag in strict]
            result = linprog([0.0]*30+[-1.0 if strict_mode else 1.0],
                             A_ub=np.array([[float(x) for x in row]+[last] for row, last in zip(A, last_column)]),
                             b_ub=np.array(list(map(float, bounds))),
                             A_eq=np.array([[float(x) for x in row]+[0.0] for row in E]),
                             b_eq=np.ones(3), bounds=[(None, None)]*31,
                             options={'time_limit': max(.01, time_limit-(time.monotonic()-started))}, method='highs')
            if not result.success or not np.all(np.isfinite(result.x)):
                continue
            can_certify = result.x[-1] <= 1e-9 if strict_mode else result.x[-1] > 1e-9
            if can_certify:
                active = [i for i, value in enumerate(result.ineqlin.marginals) if value < -1e-9]
                normalizers = [F(strict[i]) if strict_mode else F(1) for i in active]
                columns = [list(A[i])+[normalizer] for i, normalizer in zip(active, normalizers)]+[row+[F(0)] for row in E]
                matrix = Matrix([[Rational(x.numerator, x.denominator) for x in column] for column in columns]).T
                try:
                    solution, parameters = matrix.gauss_jordan_solve(Matrix([0]*30+[1]))
                    solution = solution.subs({symbol: 0 for symbol in parameters})
                    multipliers = [F(int(x.p), int(x.q)) for x in solution]
                except (ValueError, AttributeError, ZeroDivisionError):
                    multipliers = None
                if multipliers is not None:
                    y, z = multipliers[:len(active)], multipliers[len(active):]
                    mass = sum(y)
                    rhs = sum((v*bounds[i] for i, v in zip(active, y)), F(0))+sum(z)
                    strict_mass = sum(v*strict[i] for i, v in zip(active, y))
                    if mass > 0 and all(v >= 0 for v in y) and (rhs < 0 or (rhs == 0 and strict_mass > 0)):
                        used = sorted({len(descriptors[i]['vertices']) for i, v in zip(active, y)
                                       if v and descriptors[i]['kind'] == 'disjoint_subset'})
                        proof = {'format': 'ghz30-vertex-probability-farkas-v1', 'graph': graph,
                                 'mean_bound_certificates': [moment_proofs[m] for m in used],
                                 'inequalities': [{'constraint': descriptors[i], 'multiplier': str(v/mass)}
                                                  for i, v in zip(active, y) if v],
                                 'equalities': [str(v/mass) for v in z],
                                 'contradiction_margin': str(-rhs/mass),
                                 'strict_multiplier_mass': str(strict_mass/mass)}
                        if any(v and descriptors[i]['kind'] == 'closed_upper' for i, v in zip(active, y)):
                            proof['degree_two_closed_subset'] = closure
                        verify_certificate(proof)
                        return proof
            if time.monotonic()-started >= time_limit:
                return None
        if not result.success or not np.all(np.isfinite(result.x)):
            return None
        added = False
        for v in range(30):
            others = [part for k, part in enumerate(graph['parts']) if k != data[1][v]]
            alpha, beta = [sum(result.x[j] for j in part if j in data[0][v]) for part in others]
            rho = 1-alpha-beta
            if min(alpha, beta, rho) <= 0:
                continue
            reciprocal = [1/alpha, 1/beta, 1/rho]
            total = sum(reciprocal)
            if result.x[v] <= 1/total+1e-9:
                continue
            s, t = [F(value/total).limit_denominator(1000) for value in reciprocal[:2]]
            item = {'kind': 'sector', 'orbit': v, 'weights': list(map(str, (s, t, 1-s-t)))}
            if item not in descriptors:
                descriptors.append(item); added = True
        if not added:
            return None
    return None


def self_test():
    from copy import deepcopy
    from pathlib import Path
    parts = [list(range(10*k, 10*k+10)) for k in range(3)]
    edges = [(i, 10+(i+s) % 10) for i in range(10) for s in (0, 1)]
    edges += [(i, 20+(i+s) % 10) for i in range(10) for s in (0, 1)]
    edges += [(10+i, 20+(i+s) % 10) for i in range(10) for s in (3, 4)]
    graph = normalize_graph({'parts': parts, 'edges': edges})
    proof = {'format': 'ghz30-vertex-probability-farkas-v1', 'graph': graph,
             'mean_bound_certificates': [],
             'inequalities': [{'constraint': {'kind': 'weighted_degree_lower',
                                              'event_context': 0, 'source_context': 1}, 'multiplier': '1'}],
             'equalities': ['0', '2', '0'], 'contradiction_margin': '0', 'strict_multiplier_mass': '1'}
    assert verify_certificate(proof) == 0
    data = geometry(graph)
    assert constraint(graph, data, {'kind': 'lower', 'orbit': 0})[1] == 0
    bad = []
    altered = deepcopy(proof); altered['inequalities'][0]['multiplier'] = '-1'; bad.append(altered)
    altered = deepcopy(proof); altered['equalities'][1] = '1'; bad.append(altered)
    altered = deepcopy(proof); altered['strict_multiplier_mass'] = '0'; bad.append(altered)
    altered = deepcopy(proof); altered['contradiction_margin'] = '1/10'; bad.append(altered)
    altered = deepcopy(proof); altered['inequalities'] = []; bad.append(altered)
    for altered in bad:
        try:
            verify_certificate(altered)
        except AssertionError:
            pass
        else:
            raise AssertionError('false vertex-probability certificate accepted')
    historical = json.loads((Path(__file__).with_name('overnight30_reviews')/'2026-11-05T0100'/
                             '0711be2a3c6aff655ef8_certificate.json').read_text())
    moment = {'format': 'ghz30-vertex-probability-farkas-v1', 'graph': historical['graph'],
              'mean_bound_certificates': historical['witness']['mean_bound_certificates'],
              'inequalities': [{'constraint': dict(kind='disjoint_subset', **block), 'multiplier': '1/2'}
                               for block in historical['witness']['blocks']],
              'equalities': ['0', '-1/2', '0'], 'contradiction_margin': '1/32',
              'strict_multiplier_mass': '1'}
    assert verify_certificate(moment) == F(1, 32)
    wrong = deepcopy(moment)
    wrong['inequalities'][0]['constraint']['other_context'] = 2
    try:
        verify_certificate(wrong)
    except AssertionError:
        pass
    else:
        raise AssertionError('non-disjoint subset accepted')
    from event30_degree_two_cover import balanced_fixture, maximal_closed_subset
    fixture = balanced_fixture()
    vertices, choices = maximal_closed_subset(fixture)
    closure = {'vertices': vertices, 'choices': choices}
    closed_data = geometry(fixture, closure=closure)
    assert constraint(fixture, closed_data, {'kind': 'closed_upper', 'vertex': 0})[1:] == (F(1, 10), False)
    missing = deepcopy(closure)
    missing['vertices'].remove(29)
    missing['choices'] = [choice for choice in missing['choices'] if choice[0] != 29]
    changed = normalize_graph({'parts': fixture['parts'], 'edges': fixture['edges']+[[0, 12]]})
    for graph, closed in ((fixture, missing), (changed, closure)):
        try:
            geometry(graph, closure=closed)
        except AssertionError:
            pass
        else:
            raise AssertionError('induced-only or nonclosed neighborhood accepted')
    return {'strict_zero_margin_control': True, 'singleton_lower_bound': '0',
            'moment_negative_margin_control': '1/32', 'closed_subset_nonstrict_cap': '1/10',
            'false_certificates_rejected': len(bad)+3}


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', nargs='?', type=Path)
    args = parser.parse_args()
    print(json.dumps({'verified_margin': str(verify_certificate(json.loads(args.certificate.read_text())))}
                     if args.certificate else self_test(), indent=2))

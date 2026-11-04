"""Exact interval propagation for thirty-event GHZ probabilities.

No symmetry or equal probabilities are assumed. Start with the strict bounds
c_i < U_i = 1/9 (Corollary 32). For a proper nonempty neighborhood N in another
context, its probability mass lies strictly between

    max(0, 1-sum_{j not in N} U_j) and min(1, sum_{j in N} U_j).

Lemma 14 bounds c_i by F(alpha,beta,1-alpha-beta), the reciprocal sum of three
reciprocals. Its exact maximum over the resulting rectangle in the probability
simplex is rational: test equal thirds, and on each rectangle edge optimize the
other coordinate by clamping half the remaining mass. F is zero on the simplex
boundary and has no other positive interior stationary point.

If this maximum is below 1/9, its maximizer activates a strict rectangle bound,
so the new event upper bound is still strict. If the maximum is 1/9, Corollary
32 supplies strictness. Upward rounding and taking the minimum with the old
bound preserve strictness. Updates are synchronous. A context budget <=1 or
an empty strict sector domain is therefore an exact exclusion.

An empty neighborhood would make the event orthogonal to the state. A full
neighborhood forces the other neighborhood empty by triangle-freeness; these
degenerate cases are certified directly before interval propagation.

Verification uses only Python's standard library. Search failure, a fixed point,
or reaching the round limit means only that this necessary test did not exclude
the graph. These bounds apply to allowed supports even when allowed edges have
zero overlap. Do not run an assertion-based certificate verifier with python -O.
"""
import argparse
from fractions import Fraction as F
import itertools
import json
from pathlib import Path

from overnight30_screen import adjacency, base_validation, normalize_graph


def harmonic(alpha, beta):
    """Continuous extension of F to the closed two-dimensional simplex."""
    alpha, beta = F(alpha), F(beta)
    rho = 1-alpha-beta
    assert alpha >= 0 and beta >= 0 and rho >= 0
    if min(alpha, beta, rho) == 0:
        return F(0)
    return 1/(1/alpha+1/beta+1/rho)


def rectangle_maximum(a, A, b, B):
    """Exact supremum, or None when the strict sector domain is empty."""
    a, A, b, B = map(F, (a, A, b, B))
    assert 0 <= a <= 1 and 0 <= A <= 1 and 0 <= b <= 1 and 0 <= B <= 1
    if a >= A or b >= B or a+b >= 1:
        return None
    candidates = [F(0)]
    if a <= F(1, 3) <= A and b <= F(1, 3) <= B:
        candidates.append(F(1, 9))
    for alpha in (a, A):
        upper = min(B, 1-alpha)
        if b <= upper:
            beta = min(upper, max(b, (1-alpha)/2))
            candidates.append(harmonic(alpha, beta))
    for beta in (b, B):
        upper = min(A, 1-beta)
        if a <= upper:
            alpha = min(upper, max(a, (1-beta)/2))
            candidates.append(harmonic(alpha, beta))
    result = max(candidates)
    assert 0 < result <= F(1, 9)
    return result


def round_up(value, denominator):
    assert type(denominator) is int and denominator > 0
    value = F(value)
    scaled = value*denominator
    return F(-(-scaled.numerator//scaled.denominator), denominator)


def prepare(record):
    graph = normalize_graph(record)
    adj, colors = base_validation(graph)
    neighborhoods = []
    for vertex in range(30):
        entry = []
        for k, part in enumerate(graph['parts']):
            if k == colors[vertex]:
                continue
            neighbors = sorted(adj[vertex].intersection(part))
            complement = sorted(set(part)-set(neighbors))
            if not neighbors or not complement:
                return graph, None, {'kind': 'degenerate_neighborhood', 'vertex': vertex,
                                     'context': k, 'neighbors': neighbors}
            entry.append((neighbors, complement))
        neighborhoods.append(entry)
    return graph, neighborhoods, None


def context_failure(graph, upper):
    for k, part in enumerate(graph['parts']):
        bound = sum((upper[i] for i in part), F(0))
        if bound <= 1:
            return {'kind': 'context_budget', 'context': k, 'bound': str(bound)}
    return None


def next_round(neighborhoods, upper, denominator):
    updated = []
    for vertex, entries in enumerate(neighborhoods):
        intervals = [(max(F(0), 1-sum((upper[i] for i in complement), F(0))),
                      min(F(1), sum((upper[i] for i in neighbors), F(0))))
                     for neighbors, complement in entries]
        (a, A), (b, B) = intervals
        maximum = rectangle_maximum(a, A, b, B)
        if maximum is None:
            return None, {'kind': 'empty_sector', 'vertex': vertex,
                          'intervals': [[str(lo), str(hi)] for lo, hi in intervals]}
        updated.append(min(upper[vertex], round_up(maximum, denominator)))
    return updated, None


def screen(record, max_rounds=3, denominator=1000000):
    """Return a verified exclusion or inconclusive propagated bounds."""
    assert type(max_rounds) is int and max_rounds >= 0
    assert type(denominator) is int and denominator > 0
    graph, neighborhoods, witness = prepare(record)
    rounds, upper = [], [F(1, 9)]*30
    if witness is None:
        for _ in range(max_rounds):
            updated, witness = next_round(neighborhoods, upper, denominator)
            if witness is not None:
                break
            rounds.append(list(map(str, updated)))
            witness = context_failure(graph, updated)
            unchanged = updated == upper
            upper = updated
            if witness is not None or unchanged:
                break
    result = {'status': 'not_excluded', 'rounds': len(rounds),
              'upper_bounds': list(map(str, upper)),
              'context_budgets': [str(sum((upper[i] for i in part), F(0)))
                                  for part in graph['parts']]}
    if witness is not None:
        certificate = {'format': 'ghz30-probability-propagation-v1', 'graph': graph,
                       'rounding_denominator': denominator, 'rounds': rounds,
                       'witness': witness}
        verify_certificate(certificate)
        result.update(status='excluded', certificate=certificate)
    return result


def verify_certificate(certificate):
    assert certificate['format'] == 'ghz30-probability-propagation-v1'
    graph, neighborhoods, witness = prepare(certificate['graph'])
    assert graph == certificate['graph']
    denominator = certificate['rounding_denominator']
    assert type(denominator) is int and denominator > 0
    upper = [F(1, 9)]*30
    if witness is not None:
        assert not certificate['rounds'] and certificate['witness'] == witness
        return witness
    for index, recorded in enumerate(certificate['rounds']):
        assert witness is None
        updated, failure = next_round(neighborhoods, upper, denominator)
        assert failure is None
        assert recorded == list(map(str, updated))
        upper = updated
        witness = context_failure(graph, upper)
        if witness is not None:
            assert index == len(certificate['rounds'])-1
    if witness is None:
        _, witness = next_round(neighborhoods, upper, denominator)
    assert witness is not None and witness == certificate['witness']
    return witness


def self_test():
    assert rectangle_maximum(0, 1, 0, 1) == F(1, 9)
    assert rectangle_maximum(F(4, 9), F(5, 9), F(4, 9), F(5, 9)) == F(2, 27)
    assert rectangle_maximum(0, F(1, 9), 0, 1) == F(2, 27)
    assert rectangle_maximum(F(1, 2), F(2, 3), F(1, 2), F(2, 3)) is None
    assert rectangle_maximum(F(1, 3), F(1, 3), 0, 1) is None
    assert round_up(F(1, 9), 1000000) == F(13889, 125000)
    # A balanced C6 blow-up is excluded after one step.
    parts = [list(range(10*k, 10*k+10)) for k in range(3)]
    clusters = [list(range(offset, offset+5)) for offset in (0, 10, 20, 5, 15, 25)]
    edges = [(u, v) for i in range(6) for u in clusters[i] for v in clusters[(i+1) % 6]]
    excluded = screen({'parts': parts, 'edges': edges})
    assert excluded['status'] == 'excluded' and excluded['rounds'] == 1
    verify_certificate(excluded['certificate'])
    # A triangle-free (3,3)-degree support passes this necessary relaxation.
    edges = [(i, 10+(i+s) % 10) for i in range(10) for s in (0, 1, 2)]
    edges += [(i, 20+(i+s) % 10) for i in range(10) for s in (0, 1, 2)]
    edges += [(10+i, 20+(i+s) % 10) for i in range(10) for s in (3, 4, 5)]
    admitted = screen({'parts': parts, 'edges': edges})
    assert admitted['status'] == 'not_excluded'
    assert admitted['upper_bounds'] == ['1/9']*30
    # Empty and full allowed neighborhoods cannot occur with positive events.
    empty = screen({'parts': parts, 'edges': []})
    assert empty['certificate']['witness']['kind'] == 'degenerate_neighborhood'
    verify_certificate(empty['certificate'])
    import copy
    bad = copy.deepcopy(excluded['certificate'])
    bad['rounds'][0][0] = '1/9'
    try:
        verify_certificate(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('tampered propagation transcript accepted')
    # Exact dense-grid lower controls for every admissible small rational box.
    grid = [F(k, 6) for k in range(7)]
    for a, A in itertools.combinations(grid, 2):
        for b, B in itertools.combinations(grid, 2):
            bound = rectangle_maximum(a, A, b, B)
            if bound is None:
                assert a+b >= 1
                continue
            for alpha in grid:
                for beta in grid:
                    if a <= alpha <= A and b <= beta <= B and alpha+beta <= 1:
                        assert harmonic(alpha, beta) <= bound
    return {'controls': 'exact rectangle, rounding, exclusion, feasible relaxation, tampering, dense grid passed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('self-test')
    verify = sub.add_parser('verify')
    verify.add_argument('certificate', type=Path)
    run = sub.add_parser('screen')
    run.add_argument('graph', type=Path)
    run.add_argument('--rounds', type=int, default=3)
    run.add_argument('--denominator', type=int, default=1000000)
    args = parser.parse_args()
    if args.command == 'self-test':
        result = self_test()
    elif args.command == 'verify':
        result = verify_certificate(json.loads(args.certificate.read_text()))
    else:
        result = screen(json.loads(args.graph.read_text()), args.rounds, args.denominator)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

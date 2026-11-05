"""Exact one-sided regularity obstruction for thirty-event GHZ supports.

If every event in context B has q neighbors in context A, then q is3 or4.
No regularity on the A side or on either other pair support is assumed.

Let alpha_i be A-event i's B-neighbor probability mass. Then sum_A alpha_i=q.
The one-sector bound c_i<=h(alpha_i), h(a)=a(1-a)/(1+3a), and strict concavity
h''(a)=-8/(1+3a)^3 give1<=10h(q/10). Integers leaveq=2,3,4,5. The endpoint
casesq=0,10 have degenerate neighborhoods and are also impossible.

Atq=2 or5, equality forces every A probability to1/10, alpha=q/10, and the
other two sector masses beta=rho=(1-alpha)/2. For E=Pi_B,F=Pi_C, put
s=(1+alpha)/(1-alpha) and L=sE+2F. Equality gives L=I on the reduced A-space
orthogonal to psi. The projection identity

    s(LE+EL)=L^2-2L+s(s+2)E

therefore gives compression kappa=1/s^2 there. Since Epsi=psi, distinct
A-events have <g_i,E g_j>=(1-kappa)/10. Cauchy equality forces every allowed
A--B overlap to sqrt(c_i*c_b)/alpha, so every pair of A-events has common-B
probability mass (1-kappa)q^2/100:1/45 forq2 and2/9 forq5.

Each B probability is strictly below1/9. Thus every A-pair needs at least one
or three common neighbors, respectively. The B-side supplies only10*C(q,2)
common-neighbor incidences:10<45 or100<135. These contradictions allow complex
vectors and allowed zero-overlap edges; equality would force those allowed
edge overlaps nonzero. The proof uses no dimension or automorphism assumption.

More generally, every oriented context pair satisfies the strict weighted
degree bound2<D=sum_{b in B} c_b*d_A(b)<5. Indeed sum_A alpha_i=D, so the
same Jensen argument gives D in[2,5]. At either endpoint A is uniform and
the equality overlaps are sqrt(10)*sqrt(c_b)/D. Expanding psi in context A
then gives sqrt(c_b)=(d_A(b)/D)*sqrt(c_b), forcing every d_A(b)=D. This
reduces either endpoint to the already excluded one-sided regular case.

In particular, each pair support has at least20 edges. All B-degrees are
positive. If they are nonconstant, put delta=min_b d_A(b)>=1 and let M be
their sum, the pair edge count. Since c_b<1/9 and at least one d_A(b)>delta,

    D=delta+sum_b c_b*(d_A(b)-delta)<(M-delta)/9<=(M-1)/9.

M<=19 would contradict D>2. If the degrees are constant, they are3 or4,
so M is30 or40. These corollaries need neither regularity nor symmetry;
the graph witness API below currently implements the one-sided regular test.
"""
if not __debug__:
    raise RuntimeError('Exact one-sided regularity verification requires assertions; do not use -O.')

import argparse
from fractions import Fraction as F
import itertools
import json
from math import comb
from pathlib import Path

from overnight30_screen import base_validation, normalize_graph


def h(alpha):
    alpha = F(alpha)
    assert 0 <= alpha <= 1
    return alpha*(1-alpha)/(1+3*alpha)


def arithmetic():
    budgets = {q: 10*h(F(q, 10)) for q in range(11)}
    assert [q for q, bound in budgets.items() if bound >= 1] == [2, 3, 4, 5]
    for numerator in range(101):
        D = F(numerator, 10)
        assert 10*h(D/10)-1 == -(D-2)*(D-5)/(10+3*D)
    for edge_count in range(10, 20):
        for minimum in range(1, edge_count//10+1):
            assert F(edge_count-minimum, 9) <= F(edge_count-1, 9) <= 2
    from event30_regular_blocks import add, multiply, scale
    cases = []
    for q in (2, 5):
        assert budgets[q] == 1
        alpha = F(q, 10)
        beta = rho = (1-alpha)/2
        assert 1/(1/alpha+1/beta+1/rho) == F(1, 10)
        s, t = (1+alpha)/(1-alpha), F(2)
        E, G = {'E': F(1)}, {'F': F(1)}
        L = add(scale(E, s), scale(G, t))
        assert scale(add(multiply(L, E), multiply(E, L)), s) == add(multiply(L, L), scale(L, -t), scale(E, s*(s+t)))
        kappa = (1-t)/(s*(2-s-t))
        assert kappa == 1/(s*s)
        assert kappa+(1-kappa)/10 == F(1, q)
        common_mass = (1-kappa)*F(q*q, 100)
        needed = int(9*common_mass)+1
        available = 10*comb(q, 2)
        required = comb(10, 2)*needed
        assert available < required
        cases.append({'regular_degree': q, 'sector_masses': list(map(str, (alpha, beta, rho))),
                      'projection_coefficient': str(s), 'compression': str(kappa),
                      'common_neighbor_mass': str(common_mass), 'required_common_neighbors': needed,
                      'available_pair_incidences': available, 'required_pair_incidences': required})
    assert [case['common_neighbor_mass'] for case in cases] == ['1/45', '2/9']
    return {'remaining_one_sided_degrees': [3, 4], 'strict_weighted_degree_interval': '(2,5)',
            'minimum_pair_edges': 20, 'jensen_budgets': {q: str(v) for q, v in budgets.items()},
            'equality_case_exclusions': cases}


def find_obstruction(record):
    graph = normalize_graph(record)
    adj, _ = base_validation(graph)
    parts = list(map(set, graph['parts']))
    for event_context, regular_context in itertools.permutations(range(3), 2):
        degrees = {len(adj[v] & parts[event_context]) for v in parts[regular_context]}
        if len(degrees) == 1:
            degree = next(iter(degrees))
            if degree not in (3, 4):
                return {'format': 'ghz30-one-sided-regular-v1', 'graph': graph,
                        'witness': {'event_context': event_context, 'regular_context': regular_context,
                                    'regular_degree': degree}}
    return None


def verify_certificate(certificate):
    assert certificate['format'] == 'ghz30-one-sided-regular-v1'
    graph = normalize_graph(certificate['graph'])
    assert graph == certificate['graph']
    adj, _ = base_validation(graph)
    witness = certificate['witness']
    event_context, regular_context = witness['event_context'], witness['regular_context']
    degree = witness['regular_degree']
    assert type(event_context) is int and 0 <= event_context < 3
    assert type(regular_context) is int and 0 <= regular_context < 3 and event_context != regular_context
    assert type(degree) is int and 0 <= degree <= 10 and degree not in (3, 4)
    events = set(graph['parts'][event_context])
    assert all(len(adj[v] & events) == degree for v in graph['parts'][regular_context])
    return True


def self_test():
    report = arithmetic()
    parts = [list(range(10*k, 10*k+10)) for k in range(3)]
    # Every B-event has two A-neighbors, but the A-side degrees are1,2,3.
    edges = [[i, 10+j] for j in range(10) for i in (j, (j+1) % 10)]
    edges.remove([1, 10])
    edges.append([2, 10])
    edges += [[i, 20+i] for i in range(10)]
    edges += [[10+(i+4) % 10, 20+i] for i in range(10)]
    graph = normalize_graph({'parts': parts, 'edges': edges})
    certificate = find_obstruction(graph)
    assert certificate['witness'] == {'event_context': 0, 'regular_context': 1, 'regular_degree': 2}
    verify_certificate(certificate)
    # The worker reconstructs this standard global format for compact cuts.
    from overnight30_screen import verify_certificate as verify_global
    from overnight30_symmetry import dominance_clause, replay_structural_dominance
    from event30_hall import graph6
    global_certificate = {'format': 'overnight30-certificate-v1', 'graph': graph,
                          'method': 'one_sided_regular', 'witness': certificate['witness']}
    assert verify_global(global_certificate)
    orbits = [[(u, v)] for u in range(30) for v in range(u+1, 30) if u//10 != v//10]
    cut = {'clause': dominance_clause(graph, orbits), 'reason': 'verified_structural_dominance',
           'witness': {'graph6': graph6(graph), 'method': 'one_sided_regular',
                       'witness': certificate['witness']}}
    assert replay_structural_dominance(cut, orbits) == cut['clause']
    adj, _ = base_validation(graph)
    assert {len(adj[i] & set(parts[1])) for i in parts[0]} == {1, 2, 3}
    # All three pairs are3-regular and triangle-free: this necessary test passes.
    edges = [[i, 10+(i+s) % 10] for i in range(10) for s in (0, 1, 2)]
    edges += [[i, 20+(i+s) % 10] for i in range(10) for s in (0, 1, 2)]
    edges += [[10+i, 20+(i+s) % 10] for i in range(10) for s in (3, 4, 5)]
    assert find_obstruction({'parts': parts, 'edges': edges}) is None
    import copy
    for changes in ({'regular_degree': 3}, {'event_context': 1, 'regular_context': 0}, {'regular_degree': 5}):
        bad = copy.deepcopy(certificate)
        bad['witness'].update(changes)
        try:
            verify_certificate(bad)
        except AssertionError:
            pass
        else:
            raise AssertionError('tampered oriented one-sided certificate accepted')
    report['controls'] = 'irregular A-side, allowed3-regular relaxation, compact global dominance replay, orientation and degree tampering passed'
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('self-test')
    verify = sub.add_parser('verify')
    verify.add_argument('certificate', type=Path)
    check = sub.add_parser('check')
    check.add_argument('graph', type=Path)
    args = parser.parse_args()
    if args.command == 'self-test':
        result = self_test()
    elif args.command == 'verify':
        result = {'verified': verify_certificate(json.loads(args.certificate.read_text()))}
    else:
        result = find_obstruction(json.loads(args.graph.read_text()))
        if result is not None:
            verify_certificate(result)
    print(json.dumps(result, indent=2))

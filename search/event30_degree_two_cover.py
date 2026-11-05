"""Exact closed degree-two neighborhood obstruction (Corollary 57b).

Scope: thirty positive-probability events in three contexts of ten, with
triangle-free allowed support. No dimension, symmetry, uniformity, or real
overlap assumption is made. Extra allowed edges may have zero overlap.

Suppose S contains a whole context and, for every vertex v in S, at least one
other-context allowed neighborhood N(v) has size at most two and lies wholly
in S. Degrees here are those of the ORIGINAL support, never its subgraph on S.
Let M=max_{v in S} c_v. Every event has c_v<1/9. For a maximizing vertex,
the selected neighbor mass alpha is at most 2M<2/9. The one-sector bound

    c_v <= h(alpha),  h(a)=a(1-a)/(1+3a),

and h'(a)=(1+a)(1-3a)/(1+3a)^2>0 on [0,2/9] imply

    M <= h(2M),   h(2M)-M = M(1-10M)/(1+6M).

Thus every probability in S is at most 1/10. The full context A contained
in S must be uniform. For each of its events, equality now forces its chosen
neighbor mass alpha=1/5, exactly two neighbors of probability 1/10, and the
other sector masses beta=rho=2/5. The other neighborhood need not have four
vertices: its individual probabilities have not been assumed uniform.

With event phases chosen to give positive state overlaps, saturation gives

    g_v=sqrt(1/10)*(15/2 psi_small + 5 psi_large - 5/2 psi),

where psi_small and psi_large are state projections on the two respective
neighbor sets. Let I_B,I_C be the probability masses of the intersections
of two A-events' B- and C-neighborhoods. Orthogonality of two events whose
small sides both lie in B gives 3 I_B+I_C=1/5. Since I_B is an integer
multiple of 1/10, their two-element B-neighborhoods are disjoint. For events
whose small sides lie in opposite contexts, expanding each event in turn
gives

    3 I_B+4 I_C=2/5,   4 I_B+3 I_C=2/5.

Hence I_B=I_C=2/35, impossible: each intersection lies in a small neighbor
set whose members all have probability 1/10. All ten A-events consequently
have the same small-side orientation, with pairwise disjoint two-element
sets inside a ten-event context. This is impossible.

Taking S to be all thirty vertices excludes a global degree-two cover. In
particular, every realizable support has a vertex with at least three
neighbors in each other context. Triangle-freeness makes any three of those
neighbors on either side a 3-by-3 zero rectangle in that context pair.

The discovery routine starts with all vertices and deterministically peels
those having no original low-degree neighborhood still contained in S.
Any closed subset remains throughout this iteration, so the fixed point is
the greatest such subset. A graph is excluded exactly when it contains a
whole context. The saved witness includes the graph, S, a contained context,
and each vertex's chosen other context. Verification independently rechecks
the original degree and complete-neighborhood containment; it does not trust
the peeling history. All controls and certificate replay use standard Python.
"""
if not __debug__:
    raise RuntimeError('Exact degree-two closure verification requires assertions; do not use python -O.')

import argparse
from copy import deepcopy
from fractions import Fraction as F
import itertools
import json
from pathlib import Path

from overnight30_screen import base_validation, certificate, normalize_graph


def maximal_closed_subset(graph, adj=None, colors=None):
    """Return the greatest closed subset and its deterministic local choices."""
    if adj is None or colors is None:
        adj, colors = base_validation(graph)
    parts = list(map(set, graph['parts']))
    # These neighborhoods never change when vertices are peeled.
    options = {
        v: [(k, adj[v] & part) for k, part in enumerate(parts)
            if k != colors[v] and len(adj[v] & part) <= 2]
        for v in range(graph['n'])}
    remaining = set(range(graph['n']))
    while True:
        removed = {v for v in remaining
                   if not any(neighbors <= remaining for _, neighbors in options[v])}
        if not removed:
            break
        remaining.difference_update(removed)
    choices = [[v, next(k for k, neighbors in options[v] if neighbors <= remaining)]
               for v in sorted(remaining)]
    return sorted(remaining), choices


def find_obstruction(record):
    graph = normalize_graph(record)
    if graph['n'] != 30 or list(map(len, graph['parts'])) != [10, 10, 10]:
        return None
    adj, colors = base_validation(graph)
    vertices, choices = maximal_closed_subset(graph, adj, colors)
    remaining = set(vertices)
    for context, part in enumerate(graph['parts']):
        if set(part) <= remaining:
            return certificate(graph, 'degree_two_closure', vertices=vertices,
                               context=context, choices=choices)
    return None


def verify_certificate(proof):
    assert proof['format'] == 'overnight30-certificate-v1'
    assert proof['method'] == 'degree_two_closure'
    graph = normalize_graph(proof['graph'])
    assert graph == proof['graph']
    adj, colors = base_validation(graph)
    witness = proof['witness']
    vertices, context, choices = witness['vertices'], witness['context'], witness['choices']
    assert isinstance(vertices, list) and vertices
    assert all(type(v) is int and 0 <= v < 30 for v in vertices)
    assert vertices == sorted(set(vertices))
    remaining = set(vertices)
    assert type(context) is int and 0 <= context < 3
    assert set(graph['parts'][context]) <= remaining
    assert isinstance(choices, list) and len(choices) == len(vertices)
    assert all(isinstance(choice, list) and len(choice) == 2 for choice in choices)
    assert [choice[0] for choice in choices] == vertices
    parts = list(map(set, graph['parts']))
    for v, k in choices:
        assert type(v) is int and type(k) is int and 0 <= k < 3 and k != colors[v]
        neighbors = adj[v] & parts[k]
        assert len(neighbors) <= 2 and neighbors <= remaining
    return True


def arithmetic():
    def h(a):
        return a*(1-a)/(1+3*a)
    assert h(F(1, 5)) == F(1, 10)
    assert (1+F(2, 9))*(1-3*F(2, 9)) > 0
    for numerator in range(1, 101):
        mass = F(numerator, 900)
        assert h(2*mass)-mass == mass*(1-10*mass)/(1+6*mass)
    small = [F(i, 10) for i in range(3)]
    assert [mass for mass in small if F(1, 5)-3*mass >= 0] == [F(0)]
    opposite = [(a, b) for a, b in itertools.product(small, repeat=2)
                if 3*a+4*b == F(2, 5) and 4*a+3*b == F(2, 5)]
    assert not opposite and F(2, 35) not in small
    assert 10*2 > 10
    return {'closed_subset_probability_cap': '1/10',
            'opposite_orientation_intersection_mass': '2/35',
            'scope': 'all dimensions; complex overlaps; thirty positive events'}


def balanced_fixture():
    """Alternating degree-two/four links on six five-vertex context orbits."""
    parts = [list(range(10*k, 10*k+10)) for k in range(3)]
    cycle = [0, 2, 4, 1, 3, 5]
    edges = []
    for position, a in enumerate(cycle):
        b = cycle[(position+1) % len(cycle)]
        degree = 2 if position % 2 == 0 else 4
        edges.extend((5*a+i, 5*b+(i+s) % 5)
                     for i in range(5) for s in range(degree))
    return normalize_graph({'parts': parts, 'edges': edges})


def positive_controls():
    """Replay rational H33/H36 primals and check they stay outside this scope."""
    from overnight30_screen import verify_certificate as verify_global
    encoded = Path(__file__).with_name('h33.g6').read_text().strip()
    if encoded.startswith('>>graph6<<'):
        encoded = encoded[10:]
    assert ord(encoded[0])-63 == 33
    bits = [(ord(char)-63) >> shift & 1 for char in encoded[1:]
            for shift in range(5, -1, -1)]
    pairs = [(i, j) for j in range(33) for i in range(j)]
    assert len(bits) >= len(pairs)
    edges = [edge for edge, bit in zip(pairs, bits) if bit]
    h33 = normalize_graph({'parts': [list(range(9)), list(range(9, 21)), list(range(21, 33))],
                           'edges': edges})
    proof33 = certificate(h33, 'theta_primal',
                          weights=[str(F(1, 9) if v < 9 else F(1, 12)) for v in range(33)],
                          edge_values=[[u, v, str(F(1, 36) if u < 9 or v < 9 else F(1, 24))]
                                       for u, v in h33['edges']])
    shifts = (((0, 0), (0, 3), (0, 4), (1, 4)),
              ((0, 0), (0, 3), (0, 4), (1, 4)),
              ((0, 1), (1, 0), (1, 1), (1, 3)))
    edges = [(12*k+6*a+b, 12*((k+1) % 3)+6*((a+s) % 2)+(b+t) % 6)
             for k in range(3) for a in range(2) for b in range(6)
             for s, t in shifts[k]]
    h36 = normalize_graph({'parts': [list(range(12*k, 12*k+12)) for k in range(3)],
                           'edges': edges})
    proof36 = certificate(h36, 'theta_primal', weights=['1/12']*36,
                          edge_values=[[u, v, '1/48'] for u, v in h36['edges']])
    for graph, proof in ((h33, proof33), (h36, proof36)):
        assert verify_global(proof)
        assert find_obstruction(graph) is None
        invalid = certificate(graph, 'degree_two_closure', vertices=list(range(graph['n'])),
                              context=0, choices=[[v, (v//12+1) % 3] for v in range(graph['n'])])
        try:
            verify_certificate(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError('out-of-scope graph accepted by closure verifier')
    return ['H33 exact rational primal', 'H36 exact rational primal']


def self_test():
    report = arithmetic()
    graph = balanced_fixture()
    adj, colors = base_validation(graph)
    parts = list(map(set, graph['parts']))
    for source, target in itertools.permutations(range(3), 2):
        degrees = [len(adj[v] & parts[target]) for v in parts[source]]
        assert sorted(degrees) == [2]*5+[4]*5 and sum(degrees) == 30
    proof = find_obstruction(graph)
    assert proof['witness']['vertices'] == list(range(30))
    assert verify_certificate(proof)
    from overnight30_screen import verify_certificate as verify_global
    assert verify_global(json.loads(json.dumps(proof)))
    # One added edge destroys coverage in the A0--B0 block. Its original
    # degree-two dependency cycle peels all ten vertices, leaving a proper
    # closed subset that still contains the whole third context.
    changed = normalize_graph({'parts': graph['parts'], 'edges': graph['edges']+[[0, 12]]})
    partial = find_obstruction(changed)
    assert partial['witness']['vertices'] == list(range(5, 10))+list(range(15, 30))
    assert partial['witness']['context'] == 2 and verify_certificate(partial)
    assert verify_global(json.loads(json.dumps(partial)))
    # Changing each low-degree orbit link similarly leaves no closed subset.
    nonclosed = normalize_graph({'parts': graph['parts'],
                                'edges': graph['edges']+[[0, 12], [7, 20], [15, 27]]})
    base_validation(nonclosed)
    assert maximal_closed_subset(nonclosed) == ([], [])
    assert find_obstruction(nonclosed) is None
    bad = []
    modified = deepcopy(proof); modified['graph'] = changed; bad.append(modified)
    modified = deepcopy(partial); modified['witness']['context'] = 0; bad.append(modified)
    modified = deepcopy(proof); modified['witness']['choices'][0][1] = 0; bad.append(modified)
    modified = deepcopy(proof); modified['witness']['choices'][0][1] = 2; bad.append(modified)
    modified = deepcopy(proof); modified['witness']['vertices'].remove(29)
    modified['witness']['choices'] = [item for item in modified['witness']['choices'] if item[0] != 29]
    bad.append(modified)
    modified = deepcopy(proof); modified['witness']['choices'].pop(); bad.append(modified)
    modified = deepcopy(proof); modified['witness']['vertices'].append(29); bad.append(modified)
    for modified in bad:
        try:
            verify_certificate(modified)
        except (AssertionError, ValueError):
            pass
        else:
            raise AssertionError('tampered degree-two closure certificate accepted')
    report.update(global_cover_vertices=30, proper_closed_subset_vertices=20,
                  nonclosed_fixed_point_vertices=0, tampered_certificates_rejected=len(bad),
                  positive_controls=positive_controls())
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('self-test')
    verify = sub.add_parser('verify'); verify.add_argument('certificate', type=Path)
    check = sub.add_parser('check'); check.add_argument('graph', type=Path)
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

"""Exact, phase-independent exclusion of one six-orbit C5 quotient class.

Assume a context-preserving group has six transitive orbits of size5, labeled
A0,A1,B0,B1,C0,C1, with the quotient degrees Q below. A simultaneous5-cycle
on these six orbits suffices. An equitable partition alone does not suffice.

If a quantum realization existed, average its weighted Gram matrix
K_ij=<sqrt(c_i)g_i,sqrt(c_j)g_j> over this group. Positivity, support zeros,
within-context orthogonality, normalization and the equal context sums are
linear constraints, hence survive averaging. Factoring the averaged matrix
therefore gives a realization with probabilities constant on these orbits.
The original probabilities need not be uniform; the averaged realization
can have higher dimension. A contradiction in every dimension excludes the
original realization too.

Write orbit probabilities x,1/5-x,y,1/5-y,z,1/5-z. Every B0 event has total
A-neighbor probability1/5, and every B1 event has C-neighbor probability1/5.
Lemma14 with one sector1/5 gives c<=h(1/5)=1/10, with equality only when the
remaining sectors are2/5,2/5. Thus y<=1/10 and1/5-y<=1/10 imply y=1/10.
Equality at B0 gives3/5-2z=2/5; at B1 it gives1/5+2x=2/5. Hence x=z=1/10.
All averaged probabilities are uniform.

The A1 events have degree pair(2,4), and B1 events have degree pair(4,2).
Both saturate Lemma14. On every A1--B1 edge, equality requires normalized
overlap1/2 from the A1 endpoint and1/4 from the B1 endpoint, after positive
state rephasing. Q[A1,B1]=1 ensures such an edge exists: contradiction.
Only quotient degrees enter, so every support with this group action and
quotient is excluded, regardless of the chosen cyclic edge phases.

Run directly for exact arithmetic checks and, when available, a check of the
saved C5-double search leader. An optional JSON path checks another support
with contiguous five-vertex orbits and the simultaneous cyclic generator.
No numerical theta or Gram output is used in this proof checker.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json


Q = [
    [0, 0, 1, 3, 1, 2],
    [0, 0, 1, 1, 3, 1],
    [1, 1, 0, 0, 1, 3],
    [3, 1, 0, 0, 1, 1],
    [1, 3, 1, 1, 0, 0],
    [2, 1, 3, 1, 0, 0],
]


def check_arithmetic():
    assert all(Q[i][j] == Q[j][i] for i in range(6) for j in range(6))
    assert all(Q[i][j] == 0 for i in range(6) for j in range(6) if i // 2 == j // 2)
    assert Q[2][:2] == [1, 1] and Q[3][4:] == [1, 1]
    alpha = F(1, 5)
    assert alpha * (1 - alpha) / (1 + 3 * alpha) == F(1, 10)
    # For alpha1/5, the numerator of 1/10-F(alpha,t,4/5-t)
    # is5(t-2/5)^2 after multiplication by its positive denominator.
    assert [F(5), F(-4), F(4, 5)] == [F(5), -10 * F(2, 5), 5 * F(2, 5) ** 2]
    z = (F(3, 5) - F(2, 5)) / 2
    x = (F(2, 5) - F(1, 5)) / 2
    y = F(1, 10)
    assert x == y == z == F(1, 10)
    a1_degrees = (sum(Q[1][2:4]), sum(Q[1][4:6]))
    b1_degrees = (sum(Q[3][:2]), sum(Q[3][4:6]))
    assert a1_degrees == (2, 4) and b1_degrees == (4, 2)
    for q, r in (a1_degrees, b1_degrees):
        assert F(1, q) + F(1, r) + F(1, 10 - q - r) == 1
    assert Q[1][3] > 0 and F(1, a1_degrees[0]) != F(1, b1_degrees[0])
    return {'forced_averaged_orbit_probabilities': [str(x)] * 6,
            'contradictory_edge_orbits': ['A1', 'B1'],
            'forced_overlap_values': ['1/2', '1/4']}


def check_graph(path):
    record = json.loads(path.read_text())
    graph = record.get('graph', record)
    edges = {tuple(sorted(edge)) for edge in graph['edges']}
    assert all(0 <= u < v < 30 for u, v in edges)
    parts = [list(range(10 * k, 10 * k + 10)) for k in range(3)]
    assert graph.get('parts', parts) == parts
    adjacency = [set() for _ in range(30)]
    for u, v in edges:
        assert u // 10 != v // 10
        adjacency[u].add(v)
        adjacency[v].add(u)
    assert all(not adjacency[u].intersection(adjacency[v]) for u, v in edges)
    generator = lambda u: 5 * (u // 5) + (u + 1) % 5
    assert {tuple(sorted((generator(u), generator(v)))) for u, v in edges} == edges
    orbits = [set(range(5 * i, 5 * i + 5)) for i in range(6)]
    for i, orbit in enumerate(orbits):
        assert all([len(adjacency[u] & other) for other in orbits] == Q[i] for u in orbit)
    return {'path': str(path), 'vertices': 30, 'edges': len(edges),
            'cyclic_symmetry_and_quotient_verified': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graph', type=Path, nargs='?')
    args = parser.parse_args()
    result = check_arithmetic()
    default = Path(__file__).with_name('overnight30_reviews') / '2026-11-04T2200' / 'snapshot_00023_c5_double_best.json'
    path = args.graph if args.graph is not None else default
    if args.graph is not None or path.exists():
        result['checked_graph'] = check_graph(path)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

"""Seven faithful C3 x C3 T/D action types, and an excluded DDD control.

Scope and classification proof
------------------------------
The setting is a triangle-free support graph with three independent
ten-event contexts, positive probabilities<1/9 summing to1 per context,
and the established nonempty, distinct cross-neighborhood conditions.
This starts after regular9-orbit contexts have been excluded. Each context
has type T=3+3+3+1 or D=3+3+1+1+1+1, with a faithful action of G=C3 x C3.
A moving orbit is G/K for a kernel line K. Distinct kernel lines K,L give
one diagonal-G orbit on (G/K)x(G/L), so their support block is all or none.
The same is true for a moving/fixed block. Consequently, the three vertices
of G/K can have distinct neighborhoods in another context only if that
context also contains a G/K orbit. Neighborhood injectivity in both
directions therefore forces all three contexts to have the SAME KERNEL SET.

Faithfulness leaves two or three kernels for T, and exactly two for D.
If a D occurs, every context has those same two kernels X,Y; a T repeats
one of them. Up to kernel exchange and context permutation, TTT has repeat
patterns XXX or XXY, TTD has equal or opposed repeated kernels, and TDD
and DDD each have one type. The remaining TTT type uses three distinct
kernels in every context. GL(2,3) is transitive on unordered pairs and
triples of the four kernel lines. These are the seven ACTION types below;
this is not an enumeration or exclusion of their possible support graphs.

For a context pair A,B, the number of possible support edge orbits is
o(A)o(B)+2*sum_K m_A(K)m_B(K), where o(T)=4,o(D)=6. Every pair of vertex
orbits contributes one, except two moving orbits with the same kernel
contribute three. We check the resulting counts by explicit permutations.

Further DDD necessary conditions
--------------------------------
Each same-kernel block has one or two of its three circulant edge orbits:
zero or three repeat neighborhoods. At most one of the three context
pairs per kernel can use two: two size2 subsets of Z3 already sum to Z3,
so together with the nonempty third block they create a triangle.
Full cross-kernel blocks form a matching on the six moving orbit nodes:
X_A--Y_B and X_A--Y_C would form triangles with the nonempty Y_B--Y_C
block, and likewise at the other endpoint. A fixed vertex touches each
kernel in at most one other context, hence at most two moving orbits.
Every context pair has a fixed--fixed edge: otherwise its four fixed
vertices have just three possible nonempty neighborhoods among two
moving orbits. At most two fixed vertices per context have no fixed
neighbors anywhere, since only the two split moving patterns remain.
The four fixed probabilities have total>1/3, forcing at least two>2/27
and at least three>1/18, using the individual strict upper bound1/9.

The saved DDD prism is a CONTROL, not a feasible quantum example. It
passes the structural fast screen, all same-context length3 path checks,
and local probability inequalities, but its saved rational theta dual
has bound6383/2500<3. This excludes only this graph, not the DDD action.
The helper regenerates the graph and replays that exact certificate.
No numerical solver is imported or run; default replay uses stdlib only.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path


FOLDER = Path(__file__).with_name('event30_td_controls')
KERNEL_COORDINATES = {'X': (1, 0), 'Y': (0, 1), 'Z': (1, 1)}
ACTION_TYPES = (
    ('TTT_three_kernels', ('XYZ', 'XYZ', 'XYZ'), 66),
    ('TTT_same_repeat', ('XXY', 'XXY', 'XXY'), 78),
    ('TTT_mixed_repeat', ('XXY', 'XXY', 'XYY'), 74),
    ('TTD_same_repeat', ('XXY', 'XXY', 'XY'), 86),
    ('TTD_opposed_repeat', ('XXY', 'XYY', 'XY'), 84),
    ('TDD', ('XXY', 'XY', 'XY'), 100),
    ('DDD', ('XY', 'XY', 'XY'), 120),
)


def action_generators(kernels):
    """Each letter labels a 3-orbit; the remaining context vertices are fixed."""
    generators = [list(range(30)), list(range(30))]
    for context, word in enumerate(kernels):
        assert len(word) in (2, 3) and len(set(word)) >= 2
        for orbit, kernel in enumerate(word):
            for generator, shift in zip(generators, KERNEL_COORDINATES[kernel]):
                start = 10*context+3*orbit
                for k in range(3):
                    generator[start+k] = start+(k+shift) % 3
    return generators


def support_orbit_count(generators):
    remaining = {(10*a+i, 10*b+j) for a, b in combinations(range(3), 2)
                 for i, j in product(range(10), repeat=2)}
    count = 0
    while remaining:
        start = min(remaining)
        orbit, todo = {start}, [start]
        while todo:
            u, v = todo.pop()
            for permutation in generators:
                edge = (permutation[u], permutation[v])
                if edge not in orbit:
                    orbit.add(edge)
                    todo.append(edge)
        remaining.difference_update(orbit)
        count += 1
    return count


def check_action_types():
    result = []
    for name, kernels, expected in ACTION_TYPES:
        assert all(set(word) == set(kernels[0]) for word in kernels)
        multiplicities = list(map(Counter, kernels))
        # len(word) moving orbits and 10-3*len(word) singleton orbits.
        orbit_counts = [10-2*len(word) for word in kernels]
        formula = sum(orbit_counts[a]*orbit_counts[b]
                      + 2*sum(multiplicities[a][k]*multiplicities[b][k]
                              for k in set(kernels[0]))
                      for a, b in combinations(range(3), 2))
        generators = action_generators(kernels)
        assert formula == support_orbit_count(generators) == expected
        for part in [range(10*i, 10*i+10) for i in range(3)]:
            restrictions = set()
            for a, b in product(range(3), repeat=2):
                images = list(part)
                for generator, exponent in zip(generators, (a, b)):
                    for _ in range(exponent):
                        images = [generator[v] for v in images]
                restrictions.add(tuple(images))
            assert len(restrictions) == 9  # Faithful in every context.
        result.append({'name': name, 'kernels': list(kernels), 'edge_orbits': expected})
    return result


def prism_graph():
    parts = [list(range(10*i, 10*i+10)) for i in range(3)]
    X = [[10*i+k for k in range(3)] for i in range(3)]
    Y = [[10*i+3+k for k in range(3)] for i in range(3)]
    R = [[10*i+6+k for k in range(2)] for i in range(3)]
    L = [[10*i+8+k for k in range(2)] for i in range(3)]
    r = [R[j % 3][j//3] for j in range(6)]
    ell = [L[(j+1) % 3][j//3] for j in range(6)]
    edges = set()

    def add(u, v):
        edges.add(tuple(sorted((u, v))))

    def full(left, right):
        for u, v in product(left, right):
            add(u, v)

    # The same-kernel matching phases sum to1 in Z3, preventing triangles.
    for orbit in (X, Y):
        for k in range(3):
            add(orbit[0][k], orbit[1][k])
            add(orbit[1][k], orbit[2][k])
            add(orbit[2][k], orbit[0][(k+1) % 3])
    for i in range(3):
        full(X[i], Y[(i+1) % 3])
        full(L[i], X[(i+1) % 3])
        full(R[i], Y[(i-1) % 3])
    for j in range(6):
        add(r[j], r[(j+1) % 6])
        add(ell[j], ell[(j+1) % 6])
        add(r[j], ell[j])
    labels = {}
    for name, blocks in [('X', X), ('Y', Y), ('R', R), ('L', L)]:
        for i, block in enumerate(blocks):
            for k, vertex in enumerate(block):
                labels[str(vertex)] = f'{name}_{i}_{k}'
    return {'n': 30, 'parts': parts, 'edges': [list(edge) for edge in sorted(edges)],
            'generators': action_generators(('XY', 'XY', 'XY')), 'labels': labels,
            'probabilities': [str(F(21, 200) if v % 10 < 6 else F(37, 400))
                              for v in range(30)],
            'family': 'DDD moving core plus hexagonal-prism fixed core'}


def check_prism():
    from overnight30_screen import fast_screen, normalize_graph, verify_certificate
    graph = prism_graph()
    assert graph == json.loads((FOLDER/'ddd_prism.json').read_text())
    parts, edges = graph['parts'], {tuple(edge) for edge in graph['edges']}
    adj = [set() for _ in range(30)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    assert len(edges) == 99
    assert all(not (adj[u] & adj[v]) for u, v in edges)
    path_checks = 0
    for part in parts:
        assert all(not (adj[u] & set(part)) for u in part)
        for other in parts:
            if other == part:
                continue
            neighborhoods = [frozenset(adj[v] & set(other)) for v in part]
            assert all(neighborhoods) and len(set(neighborhoods)) == 10
        for u, v in combinations(part, 2):
            assert any(adj[a] & adj[v] for a in adj[u]), (u, v)
            path_checks += 1
    for permutation in graph['generators']:
        assert {tuple(sorted((permutation[u], permutation[v]))) for u, v in edges} == edges
        assert all({permutation[u] for u in part} == set(part) for part in parts)
    probabilities = list(map(F, graph['probabilities']))
    assert all(sum(probabilities[v] for v in part) == 1 for part in parts)
    assert all(0 < p < F(1, 9) for p in probabilities)
    local_bounds = set()
    for part in parts:
        others = [set(other) for other in parts if other != part]
        for v in part:
            alpha, beta = [sum(probabilities[w] for w in adj[v] & other) for other in others]
            rho = 1-alpha-beta
            assert rho > 0
            bound = 1/(1/alpha+1/beta+1/rho)
            assert probabilities[v] < bound
            local_bounds.add((str(probabilities[v]), str(bound)))
    screen = fast_screen(graph)
    assert screen['status'] == 'not_excluded'
    assert all(row['components'] == 1 and row['perfect_matching']
               and row['bridge_sizes'] == [10] for row in screen['pair_diagnostics'])
    certificate = json.loads((FOLDER/'ddd_prism_certificate.json').read_text())
    assert certificate['graph'] == normalize_graph(graph)
    assert certificate['method'] == 'theta_dual'
    assert F(certificate['witness']['bound']) == F(6383, 2500) < 3
    assert verify_certificate(certificate)
    return {'graph_id': screen['graph_id'], 'edges': len(edges),
            'same_context_length3_path_checks': path_checks,
            'local_probability_bounds': sorted(local_bounds),
            'fast_screen_status': screen['status'], 'exact_dual_bound': '6383/2500',
            'scope': 'This saved graph is excluded; no family-wide DDD exclusion is asserted.'}


def main():
    if not __debug__:
        raise RuntimeError('Exact checks require assertions; do not run with python -O')
    print(json.dumps({'action_types': check_action_types(), 'prism_control': check_prism()}, indent=2))


if __name__ == '__main__':
    main()

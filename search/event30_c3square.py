r"""Two 9+1 context orbits are impossible when the third context has a fixed point.

The statement holds for any context-preserving automorphism group, not just
C3 x C3. Suppose it acts transitively on the nine nonfixed vertices in each
of A and B, fixing a and b, and has at least one fixed vertex in C. Other
C-orbit sizes are unrestricted. All support edges mean *allowed* nonzero
overlaps; the proof still works if some allowed overlaps actually vanish.

Use the necessary support properties from Proposition39: every cross-context
neighborhood is nonempty and proper, and vertices in one context have distinct
neighborhoods in another. Thus a's B-neighborhood is either {b} or B\{b},
and similarly for b and for each C-fixed vertex. The A--B moving block is
nonempty, since otherwise its nine rows would be identical.

If ab is an edge, a C-fixed vertex c must see one of {a},{b} and the opposite
nine-orbit. Seeing both singletons makes the fixed triangle; seeing both
nine-orbits makes a triangle through the nonempty moving block. Say c sees
{a} and B\{b}. Its residual sector mass is c_b-c_a>0. But b has only a as
an A-neighbor, and Lemma14 gives c_b<c_a, contradiction.

If ab is absent, a sees all of B\{b} and b sees all of A\{a}. Every C-fixed
vertex must see {a} and {b}: either mixed choice makes a triangle through a
or b and a moving vertex, and both moving choices also make a triangle.
Neither a nor b can touch a nonfixed C-orbit. For example an a--orbit edge
would, by triangle-freeness, leave every vertex in that orbit with no moving
B-neighbor. Its B-neighborhood would therefore be the same singleton {b},
contradicting neighborhood injectivity. Hence N_C(a)=N_C(b)=the set of
C-fixed vertices, of total probability s>0. At a the residual sector mass
is c_b-s>0, while Lemma14 at b gives c_b<s. Again a contradiction.

The strict comparisons use F(alpha,beta,rho)<alpha for three positive sector
masses, immediately from its reciprocal sum; equivalently h(t)<t for0<t<1.
No symmetry averaging, probability uniformity, or dimension bound is used.

For C3 x C3 on10, all orbit sizes are1,3,9, so C necessarily has a fixed
point. This excludes RRT as well as RRR, where R=9+1 and T=3+3+3+1.

RTT WITH THREE DISTINCT KERNEL LINES: COMPLETE REDUCTION

Let A be regular on9 points plus a fixed point a, and let B,C each have
three transitive3-orbits plus fixed points b,c. Suppose the three3-orbits in
each context have distinct kernel lines. Their kernel-line sets must agree:
between unequal kernel lines the diagonal action on two3-orbits is transitive,
so such blocks are all-or-none and cannot distinguish points. In fact those
BC blocks must be empty: nonempty A-neighborhoods are unions of cosets, and
cosets of two different lines always intersect, making every BC edge a triangle.

For each matching kernel K, the A--B_K, A--C_K and B_K--C_K cyclic supports
have sizes u,v,s in{1,2}; zero or full supports violate neighborhood
injectivity. Triangle-freeness in F3 gives u+v+s<=4. Also a transitive3-orbit
cannot have cross-degree1: each event would have c<h(1/9)=2/27, so its mass
would be<2/9, whereas the other seven events have total mass<7/9.

First b and c cannot both see A's nine-orbit. If b sees it, b cannot touch
any C3-orbit (each vertex has an A-moving neighbor), so N_C(b)={c}, and c
sees only a in A. Each C_K then has B-neighbors only in B_K, so s=2. If a
touches B_K, it cannot touch C_K, and c cannot touch B_K. Thus that C6 is an
entire BC component with a touching only its B-side. This holds for every
B-neighbor of a, leaving a unbridged, contrary to the pair-component lemma.
Therefore b and c both see only a: ab=ac=1 and bc=0.

The fixed point a must touch moving orbits in both B and C, or the singleton
neighbor bounds would give c_a<c_b<c_a (or c_c). If it touches p B-orbits
and q C-orbits its degrees are1+3p,1+3q. The strict probability bound and
positive residual sector imply their sum<=10. Hence p=q=1. The selected
kernels differ, since the corresponding BC block is nonempty. Label them
K1 for a--B_K1 and K2 for a--C_K2, and the remaining kernel K0.

For K1,K2, the absence of the opposite fixed-point edge and prohibition of
degree1 force s=2 and u=v=1. The BC component containing b is group-invariant,
so its B,C counts are multiples of3, apart from its fixed vertices. Balanced
pair components therefore force it also to contain c. With no cross-kernel
moving edges, the only possible connection is through K0. Thus b--C_K0 and
c--B_K0 are present. Also b--C_K1 must be present: otherwise K1 is an isolated
C6, whose common component vector is orthogonal to g_a. In a's expansion
against B, all B_K1 contributions cancel, leaving only b. This gives
c_a=c_b*|<g_a,g_b>|^2<=c_b, contradicting c_b<c_a. Similarly c--B_K2 is present.

All fixed edges are now forced. At K0 the only size triples are(1,1,2),
(1,1,1),(2,1,1),(1,2,1). GL(2,3) is transitive on ordered triples of distinct
kernel lines; choose their quotient coordinates as K1:x,K2:y,K0:x+y.
Independent translations of the six fiber labels normalize each singleton
A-support to{0}, and the missing element of each size2 A-support to0.
The K1,K2 BC supports are then uniquely{1,2}. For K0 with two singleton
A-supports, the BC support is{1,2} or one of{1},{2}. The latter two are
exchanged by simultaneous negation of A and all six fiber coordinates;
this preserves the other two{1,2} supports. If one A-support has size2,
the K0 BC support is uniquely{0}. Exchanging B,C together with x,y swaps
the wide-B and wide-C cases. Thus exactly three canonical support types
remain. The code below checks all48 GL matrices, all45 allowed local cyclic
triples, and the explicit global isomorphisms. The saved exact duals exclude
the canonical graphs, and hence this entire RTT three-kernel family, in
every dimension. The two-distinct-kernel RTT family is not classified here.

Running directly checks the finite fixed-point neighborhood case split,
the exact strict-bound identity, and the RTT phase classification. Optional
--verify-certificates checks the three saved rational duals using stdlib only.
An optional JSON input accepts an edge graph,
three parts, and a list of generators (permutations of range30); it verifies
all structural assumptions and returns the corresponding exact obstruction.
"""
from fractions import Fraction as Q
from itertools import product, permutations
import argparse
import json
from pathlib import Path


def fixed_cases():
    # aa denotes adjacency of c to a; !aa then means adjacency to all A-moving.
    allowed = {}
    for ab in (0, 1):
        valid = []
        for ca, cb in product((0, 1), repeat=2):
            fixed_triangle = ab and ca and cb
            moving_triangle = not ca and not cb
            a_triangle = ca and not cb and not ab
            b_triangle = cb and not ca and not ab
            if not (fixed_triangle or moving_triangle or a_triangle or b_triangle):
                valid.append((ca, cb))
        allowed[ab] = valid
    assert allowed == {0: [(1, 1)], 1: [(0, 1), (1, 0)]}
    # t-h(t)=4t^2/(1+3t)>0; numerator coefficients are exact.
    assert Q(1) + Q(3) == Q(4)
    for t in (Q(1, 9), Q(1, 3), Q(2, 3)):
        h = t * (1 - t) / (1 + 3 * t)
        assert t - h == 4 * t * t / (1 + 3 * t) > 0
    return allowed


def group_orbits(part, generators):
    unseen, result = set(part), []
    while unseen:
        orbit, todo = {min(unseen)}, [min(unseen)]
        while todo:
            u = todo.pop()
            new = {permutation[u] for permutation in generators} - orbit
            orbit.update(new)
            todo.extend(new)
        unseen.difference_update(orbit)
        result.append(sorted(orbit))
    return result


def check_graph(record):
    parts = record.get('parts', [list(range(10*k, 10*k+10)) for k in range(3)])
    assert list(map(len, parts)) == [10, 10, 10]
    assert sorted(v for part in parts for v in part) == list(range(30))
    edges = {tuple(sorted(edge)) for edge in record['edges']}
    generators = record['generators']
    assert generators and all(sorted(p) == list(range(30)) for p in generators)
    for p in generators:
        assert all({p[v] for v in part} == set(part) for part in parts)
        assert {tuple(sorted((p[u], p[v]))) for u, v in edges} == edges
    adj = [set() for _ in range(30)]
    for u, v in edges:
        assert 0 <= u < v < 30
        assert not any(u in part and v in part for part in parts)
        adj[u].add(v)
        adj[v].add(u)
    assert all(not adj[u] & adj[v] for u, v in edges)
    for left, right in product(parts, repeat=2):
        if left is right:
            continue
        neighborhoods = [tuple(sorted(adj[v] & set(right))) for v in left]
        assert all(0 < len(row) < 10 for row in neighborhoods)
        assert len(set(neighborhoods)) == 10
    orbits = [group_orbits(part, generators) for part in parts]
    regular = [k for k in range(3) if sorted(map(len, orbits[k])) == [1, 9]]
    assert len(regular) >= 2
    A, B = regular[:2]
    C = 3 - A - B
    a = next(orbit[0] for orbit in orbits[A] if len(orbit) == 1)
    b = next(orbit[0] for orbit in orbits[B] if len(orbit) == 1)
    fixed_C = {orbit[0] for orbit in orbits[C] if len(orbit) == 1}
    assert fixed_C
    moving_A, moving_B = set(parts[A]) - {a}, set(parts[B]) - {b}
    assert any(adj[v] & moving_B for v in moving_A)
    if b in adj[a]:
        c = min(fixed_C)
        assert (bool(a in adj[c]), bool(b in adj[c])) in ((False, True), (True, False))
        return {'case': 'adjacent_fixed_points', 'pair': [A, B], 'a': a, 'b': b, 'c': c}
    assert adj[a] & set(parts[C]) == fixed_C == adj[b] & set(parts[C])
    return {'case': 'nonadjacent_fixed_points', 'pair': [A, B], 'a': a, 'b': b,
            'fixed_third_vertices': sorted(fixed_C)}


def canonical_rtt(mode, phase=1):
    """Build the canonical RTT graph with coordinates K0=x+y,K1=x,K2=y."""
    assert mode in ('all_cycles', 'one_matching', 'one_wide_B', 'one_wide_C')
    assert phase in (1, 2)
    edges = set()
    def add(u, v):
        edges.add(tuple(sorted((u, v))))
    for x, y in product(range(3), repeat=2):
        point = 3*x+y
        for k, label in enumerate(((x+y) % 3, x, y)):
            for j in range(3):
                if (j != label if mode == 'one_wide_B' and k == 0 else j == label):
                    add(point, 10+3*k+j)
                if (j != label if mode == 'one_wide_C' and k == 0 else j == label):
                    add(point, 20+3*k+j)
    for k, i, j in product(range(3), repeat=3):
        allowed = i != j
        if k == 0 and mode == 'one_matching':
            allowed = (j-i) % 3 == phase
        if k == 0 and mode in ('one_wide_B', 'one_wide_C'):
            allowed = i == j
        if allowed:
            add(10+3*k+i, 20+3*k+j)
    add(9, 19)
    add(9, 29)
    for j in range(3):
        for u, v in ((9, 13+j), (9, 26+j), (19, 20+j), (19, 23+j),
                     (29, 10+j), (29, 16+j)):
            add(u, v)
    return {'n': 30, 'parts': [list(range(10*k, 10*k+10)) for k in range(3)],
            'edges': [list(edge) for edge in sorted(edges)]}


def edge_image(graph, permutation):
    return {tuple(sorted((permutation[u], permutation[v]))) for u, v in graph['edges']}


def check_rtt_classification():
    def projective(vector):
        first = next(x for x in vector if x)
        return tuple((x*first) % 3 for x in vector)  # 1 and2 are self-inverse.
    lines = [(1, 0), (0, 1), (1, 1), (1, 2)]
    images, matrices = set(), []
    for a, b, c, d in product(range(3), repeat=4):
        if (a*d-b*c) % 3:
            matrices.append((a, b, c, d))
            images.add(tuple(projective(((x*a+y*c) % 3, (x*b+y*d) % 3))
                             for x, y in ((1, 0), (0, 1), (1, 1))))
    assert len(matrices) == 48 and images == set(permutations(lines, 3))
    subsets = [frozenset(i for i in range(3) if mask & (1 << i)) for mask in range(1, 7)]
    normalized, admissible = set(), 0
    for S, T, U in product(subsets, repeat=3):
        if U & {(s-t) % 3 for s in S for t in T}:
            continue
        admissible += 1
        assert len(S)+len(T)+len(U) <= 4
        offset_B = next(iter(S if len(S) == 1 else set(range(3))-S))
        offset_C = next(iter(T if len(T) == 1 else set(range(3))-T))
        Sn = tuple(sorted((s-offset_B) % 3 for s in S))
        Tn = tuple(sorted((t-offset_C) % 3 for t in T))
        Un = tuple(sorted((u+offset_C-offset_B) % 3 for u in U))
        normalized.add((Sn, Tn, Un))
    assert admissible == 45
    assert normalized == {((0,), (0,), (1,)), ((0,), (0,), (2,)),
                          ((0,), (0,), (1, 2)), ((1, 2), (0,), (0,)),
                          ((0,), (1, 2), (0,))}
    negation = list(range(30))
    swap = list(range(30))
    for x, y in product(range(3), repeat=2):
        negation[3*x+y] = 3*((-x) % 3)+(-y) % 3
        swap[3*x+y] = 3*y+x
    for offset in (10, 20):
        for k, j in product(range(3), repeat=2):
            negation[offset+3*k+j] = offset+3*k+(-j) % 3
            swap[offset+3*k+j] = (30-offset)+3*({0: 0, 1: 2, 2: 1}[k])+j
    swap[19], swap[29] = 29, 19
    assert sorted(negation) == sorted(swap) == list(range(30))
    assert edge_image(canonical_rtt('one_matching', 1), negation) == {
        tuple(edge) for edge in canonical_rtt('one_matching', 2)['edges']}
    assert edge_image(canonical_rtt('one_wide_B'), swap) == {
        tuple(edge) for edge in canonical_rtt('one_wide_C')['edges']}
    for mode in ('all_cycles', 'one_matching', 'one_wide_B'):
        graph = canonical_rtt(mode)
        adj = [set() for _ in range(30)]
        for u, v in graph['edges']:
            assert u // 10 != v // 10
            adj[u].add(v)
            adj[v].add(u)
        assert all(not adj[u] & adj[v] for u, v in graph['edges'])
        for dx, dy in ((1, 0), (0, 1)):
            translation = list(range(30))
            for x, y in product(range(3), repeat=2):
                translation[3*x+y] = 3*((x+dx) % 3)+(y+dy) % 3
            for offset in (10, 20):
                for k, j in product(range(3), repeat=2):
                    shift = ((dx+dy) % 3, dx, dy)[k]
                    translation[offset+3*k+j] = offset+3*k+(j+shift) % 3
            assert edge_image(graph, translation) == {tuple(edge) for edge in graph['edges']}
    return {'GL2_matrices': 48, 'ordered_kernel_triples': len(images),
            'admissible_local_support_triples': admissible,
            'normalized_local_phase_types': len(normalized), 'canonical_global_graphs': 3}


def verify_rtt_certificates():
    from overnight30_screen import verify_certificate
    folder = Path(__file__).with_name('event30_c3square_certificates')
    results = []
    for mode in ('all_cycles', 'one_matching', 'one_wide_B'):
        graph = json.loads((folder/(mode+'_graph.json')).read_text())
        cert = json.loads((folder/(mode+'_certificate.json')).read_text())
        assert graph == canonical_rtt(mode) == cert['graph']
        assert cert['method'] == 'theta_dual' and Q(cert['witness']['bound']) < 3
        verify_certificate(cert)
        results.append({'graph': mode, 'edges': len(graph['edges']),
                        'exact_theta_upper_bound': cert['witness']['bound']})
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('graph', type=Path, nargs='?')
    parser.add_argument('--verify-certificates', action='store_true')
    args = parser.parse_args()
    result = {'allowed_fixed_point_cases': fixed_cases(),
              'scope': 'two transitive9+fixed contexts and a fixed vertex in the third'}
    result['RTT_three_distinct_kernel_classification'] = check_rtt_classification()
    if args.verify_certificates:
        result['RTT_exact_certificates'] = verify_rtt_certificates()
    if args.graph is not None:
        result['verified_obstruction'] = check_graph(json.loads(args.graph.read_text()))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()

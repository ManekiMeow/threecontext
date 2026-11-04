"""Exact, cheap orbit restrictions for thirty-event GHZ support graphs.

Run with standard Python: python search/event30_orbits.py

These are necessary combinatorial conditions, not quantum realizations or a
nonexistence theorem for whole C2/C3 families. No Gram symmetrization or rank
argument is used, so the restrictions apply in every dimension over C.

Every neighbourhood in another context is nonempty and proper. Nonemptiness
follows from positive event probability. If N_A(b)=A and b~c, triangle-freeness
forces N_A(c)=empty, impossible. Distinct rows then embed context B
equivariantly into the nonempty proper subsets of context A. In particular,
fixed_B(K) <= 2**orbits_A(K)-2 for every subgroup K.

For an action with exactly one p-cycle per context and all remaining vertices
fixed, the three moving-to-moving blocks are cyclic supports S0,S1,S2. Each
is nonempty/proper and has no nonzero translation stabilizer, since otherwise
two moved vertices have equal neighbourhoods in that context. A length-three
path between two vertices of one moving orbit cannot use a fixed internal
vertex: invariance would then close an original triangle. Therefore, with
T=S0+S1+S2 modulo p, the contraction condition gives
    0 not in T, and T union(-T) = all nonzero residues.
For p=2 this forces C6. For p=3 it gives C9 or K3,3 with a perfect matching
subdivided. For p=5 it excludes250 of the1725 triangle-free moving patterns.
"""
from collections import Counter
from itertools import combinations, product
from math import lcm


def partitions(n, largest=None):
    if n == 0:
        yield ()
        return
    for first in range(min(n, largest or n), 0, -1):
        for rest in partitions(n-first, first):
            yield (first,)+rest


def permutation(cycles):
    answer, offset = [], 0
    for length in cycles:
        answer.extend(offset+(i+1) % length for i in range(length))
        offset += length
    return answer


def proper_subset_orbits(cycles):
    perm = permutation(cycles)
    unseen = set(range(1, (1 << len(perm))-1))
    counts = Counter()
    while unseen:
        start = min(unseen)
        orbit, mask = [], start
        while mask not in orbit:
            orbit.append(mask)
            mask = sum(1 << perm[i] for i in range(len(perm)) if mask >> i & 1)
        assert mask == start
        unseen.difference_update(orbit)
        counts[len(orbit)] += 1
    return counts


def cyclic_type_pairs():
    types = list(partitions(10))
    assert len(types) == 42
    capacities = {kind: proper_subset_orbits(kind) for kind in types}
    excluded = []
    for A, B in product(types, repeat=2):
        if lcm(*A) != lcm(*B):
            continue  # already excluded by faithfulness on each context
        wanted = Counter(B)
        if any(wanted[d] > capacities[A][d] for d in wanted):
            excluded.append((A, B))
    assert excluded == [
        ((10,), (5, 2, 2, 1)),
        ((10,), (5, 2, 1, 1, 1)),
        ((6, 4), (4, 3, 1, 1, 1)),
        ((5, 5), (5, 1, 1, 1, 1, 1)),
    ]
    print("faithful cyclic types on10 vertices: equivariant-neighbourhood exclusions")
    for A, B in excluded:
        print(" ", A, "->", B)


def moved_graph(p, supports):
    adj = [set() for _ in range(3*p)]
    for k, support in enumerate(supports):
        for x, shift in product(range(p), support):
            u, v = k*p+x, ((k+1) % 3)*p+(x+shift) % p
            adj[u].add(v)
            adj[v].add(u)
    return adj


def connected(adj):
    seen, pending = {0}, [0]
    while pending:
        u = pending.pop()
        for v in adj[u]-seen:
            seen.add(v)
            pending.append(v)
    return len(seen) == len(adj)


def single_cycle_patterns(p):
    supports = [frozenset(i for i in range(p) if bits >> i & 1)
                for bits in range(1, (1 << p)-1)]
    supports = [S for S in supports if not any(
        {(x+t) % p for x in S} == S for t in range(1, p))]
    counts, triangle_free, failing_path = Counter(), 0, 0
    for S in product(supports, repeat=3):
        total = {(a+b+c) % p for a,b,c in product(*S)}
        if 0 in total:
            continue
        triangle_free += 1
        if total | {(-x) % p for x in total} != set(range(1,p)):
            failing_path += 1
            continue
        shape = tuple(sorted(map(len, S)))
        counts[shape] += 1
        if p in (2,3):
            adj = moved_graph(p, S)
            assert connected(adj)
            assert all(not (adj[u] & adj[v]) for u in range(3*p) for v in adj[u])
            if shape == (1,1,1):
                assert all(len(neighbours) == 2 for neighbours in adj)  # C_(3p)
            else:
                assert p == 3 and shape == (1,1,2)
                low = {u for u in range(9) if len(adj[u]) == 2}
                high = set(range(9))-low
                assert len(low) == 3 and all(len(adj[u]) == 3 for u in high)
                suppressed = {tuple(sorted((u,v))) for u in high for v in adj[u] & high}
                suppressed.update(tuple(sorted(adj[u])) for u in low)
                assert len(suppressed) == 9
                assert all(sum(u in edge for edge in suppressed) == 3 for u in high)
                assert all(u//3 != v//3 for u,v in suppressed)
                assert len({u//3 for u in high}) == 2  # exactly K3,3
    expected = {
        2: (4, 0, {(1,1,1): 4}),
        3: (45, 0, {(1,1,1): 18, (1,1,2): 27}),
        5: (1725, 250, {(1,1,2): 300, (1,1,3): 300, (1,1,4): 75,
                       (1,2,2): 450, (1,2,3): 300, (2,2,2): 50}),
    }
    assert (triangle_free, failing_path, dict(counts)) == expected[p]
    print(f"single C{p} orbit per context: triangle-free {triangle_free}, "
          f"contraction failures {failing_path}, surviving size types {dict(sorted(counts.items()))}")


if __name__ == "__main__":
    cyclic_type_pairs()
    for prime in (2,3,5):
        single_cycle_patterns(prime)

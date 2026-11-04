"""Exact, solver-independent rejection of a four-context KS skeleton.

Generate the historical d=6 skeleton with the SAT dependency, then verify its
combinatorics and a bipartite-neighbourhood obstruction using only Python's
standard library. Uniform blow-ups preserve the obstruction and supply examples
above the current dimension/degree lower bounds.

  python ks/skeleton_certificate.py generate ks/skeleton_d6_certificate.json
  python ks/skeleton_certificate.py verify ks/skeleton_d6_certificate.json --blowup 4
  python ks/skeleton_certificate.py shrikhande ks/shrikhande_blowup_certificate.json
  python ks/skeleton_certificate.py verify ks/shrikhande_blowup_certificate.json --blowup 5

Edges mean potentially NON-orthogonal. The rejection remains valid when some
allowed overlaps vanish; no nonzero-overlap assumptions are used by the proof.
"""
import argparse
import contextlib
from functools import lru_cache
import io
import itertools
import json
from pathlib import Path
import runpy
import sys


def adjacency(d, edges):
    adj = [set() for _ in range(4*d)]
    for u, v in edges:
        if not (0 <= u < v < 4*d) or u//d == v//d:
            raise ValueError("invalid cross-basis edge")
        adj[u].add(v)
        adj[v].add(u)
    return adj


def bipartition(adj, vertices):
    vertices = set(vertices)
    color = {}
    for root in sorted(vertices):
        if root in color:
            continue
        color[root] = 0
        stack = [root]
        while stack:
            u = stack.pop()
            for v in adj[u] & vertices:
                if v not in color:
                    color[v] = 1-color[u]
                    stack.append(v)
                elif color[v] == color[u]:
                    return None
    return [sorted(v for v in color if color[v] == k) for k in range(2)]


def independence_number(adj):
    """Exact bitset recursion on the small base graph, not its blow-up."""
    closed = [(1 << v) | sum(1 << u for u in adj[v]) for v in range(len(adj))]

    @lru_cache(None)
    def alpha(mask):
        if not mask:
            return 0
        v = max((i for i in range(len(adj)) if mask >> i & 1),
                key=lambda i: (mask & closed[i]).bit_count())
        return max(alpha(mask & ~(1 << v)), 1+alpha(mask & ~closed[v]))

    return alpha((1 << len(adj))-1)


def verify(data, blowup=1):
    if blowup < 1:
        raise ValueError("blowup must be positive")
    d = data["dimension"]
    adj = adjacency(d, data["edges"])
    # Independently check the full historical K4/C2/C4 skeleton requirements.
    bases = [range(k*d, (k+1)*d) for k in range(4)]
    for T in itertools.product(*bases):
        if all(v in adj[u] for u, v in itertools.combinations(T, 2)):
            raise ValueError("K4 found")
    for u, v in itertools.combinations(range(4*d), 2):
        for k, basis in enumerate(bases):
            if k in (u//d, v//d):
                continue
            common = blowup*len(adj[u] & adj[v] & set(basis))
            if v in adj[u] and common == 0:
                raise ValueError("C2 violation")
            if v not in adj[u] and common == 1:
                raise ValueError("C4 violation")
    for ks in itertools.combinations(range(4), 3):
        if not any(all(v in adj[u] for u, v in itertools.combinations(T, 2))
                   for T in itertools.product(*(bases[k] for k in ks))):
            raise ValueError("a triple of bases has no triangle")
    w = data["witness_vertex"]
    left, right = map(set, data["neighbourhood_bipartition"])
    if left & right or left | right != adj[w]:
        raise ValueError("certificate does not partition the neighbourhood")
    for part in (left, right):
        if any(adj[u] & part for u in part):
            raise ValueError("certificate is not a bipartition")
    minimum = min(len(adj[u] & set(bases[k])) for u in range(4*d)
                  for k in range(4) if k != u//d)
    total = min(map(len, adj))
    # Cloning vertices multiplies the above common-neighbour counts. Two clones
    # of the same vertex have >= minimum*blowup common neighbours per other basis.
    if minimum*blowup < 2:
        raise ValueError("clone C4 condition not established")
    print(f"exact: d={d*blowup}, {4*d*blowup} vertices; K4-free, C2/C4/TRI hold")
    print(f"min neighbours: {minimum*blowup} per other basis, {total*blowup} total")
    alpha = independence_number(adj)*blowup
    print(f"maximum orthogonal clique: {alpha} (target dimension {d*blowup})")
    print(f"witness neighbourhood: {len(adj[w])*blowup} vertices, bipartite")
    print("quantum realization IMPOSSIBLE: the three contexts require total Born")
    print("weight 3, but the two orthogonal colour classes allow at most 2 (Bessel).")


def generate(path):
    old_args = sys.argv[:]
    try:
        sys.argv = ["skeleton_sat.py", "6", "12", "4"]
        with contextlib.redirect_stdout(io.StringIO()):
            state = runpy.run_path(str(Path(__file__).with_name("skeleton_sat.py")))
    finally:
        sys.argv = old_args
    if not state["res"]:
        raise RuntimeError("historical skeleton unexpectedly UNSAT")
    model = {v for v in state["s"].get_model() if v > 0}
    edges = []
    for k, l in itertools.combinations(range(4), 2):
        for i, j in itertools.product(range(6), repeat=2):
            if state["pool"].obj2id[("e", k, i, l, j)] in model:
                edges.append([6*k+i, 6*l+j])
    state["s"].delete()
    adj = adjacency(6, edges)
    for w in range(24):
        partition = bipartition(adj, adj[w])
        if partition is not None:
            break
    else:
        raise RuntimeError("this solver's skeleton has no bipartite neighbourhood")
    data = {"dimension": 6, "source": "skeleton_sat.py 6 12 4 (historical relaxed bounds)",
            "edges": edges, "witness_vertex": w, "neighbourhood_bipartition": partition}
    verify(data)
    path.write_text(json.dumps(data, indent=2) + "\n")


def shrikhande(path):
    # Reorder Z4^2 by parity so the four covering bases are contiguous.
    points = [(x, y) for a, b in itertools.product(range(2), repeat=2)
              for x, y in itertools.product(range(4), repeat=2) if (x % 2, y % 2) == (a, b)]
    steps = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    edges = [[i, j] for i, p in enumerate(points) for j, q in enumerate(points) if i < j
             and ((q[0]-p[0]) % 4, (q[1]-p[1]) % 4) in steps]
    adj = adjacency(4, edges)
    data = {"dimension": 4, "source": "Shrikhande Cayley graph on Z4^2, parity classes",
            "recommended_blowup": 5, "edges": edges, "witness_vertex": 0,
            "neighbourhood_bipartition": bipartition(adj, adj[0])}
    verify(data, blowup=5)
    path.write_text(json.dumps(data, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("generate", "shrikhande", "verify"))
    parser.add_argument("path", type=Path)
    parser.add_argument("--blowup", type=int, default=1)
    args = parser.parse_args()
    if args.mode == "generate":
        generate(args.path)
    elif args.mode == "shrikhande":
        shrikhande(args.path)
    else:
        verify(json.loads(args.path.read_text()), args.blowup)


if __name__ == "__main__":
    main()

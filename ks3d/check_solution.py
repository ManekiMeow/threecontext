"""Check a SAT+nauty 'Solution' bitstring: decode the graph, test 010-colourability, and try to
embed it (seed rays fixed; other rays forced by cross products of two known orthogonal rays).

usage: python3 check_solution.py SEED_VECTORS_FILE BITSTRING
"""
import itertools
import sys
from math import gcd

from pysat.solvers import Cadical153


def canon(v):
    g = 0
    for x in v:
        g = gcd(g, abs(x))
    if g == 0:
        return (0, 0, 0)
    v = tuple(x // g for x in v)
    return v if next(x for x in v if x) > 0 else tuple(-x for x in v)


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def main():
    seed = [tuple(map(int, l.split())) for l in open(sys.argv[1]) if l.strip()]
    bits = sys.argv[2].strip()
    n = next(k for k in range(2, 80) if k * (k - 1) // 2 == len(bits))
    E, pos = set(), 0
    for j in range(n):          # column by column: (0,1), (0,2), (1,2), (0,3), ...
        for i in range(j):
            if bits[pos] == "1":
                E.add((i, j))
            pos += 1
    adj = lambda a, b: (min(a, b), max(a, b)) in E
    T = [t for t in itertools.combinations(range(n), 3) if adj(t[0], t[1]) and adj(t[0], t[2]) and adj(t[1], t[2])]
    deg = [sum(adj(i, j) for j in range(n) if j != i) for i in range(n)]
    print(f"n={n} edges={len(E)} triangles={len(T)} min degree={min(deg)}")
    cl = [[-(i + 1), -(j + 1)] for i, j in E] + [[a + 1, b + 1, c + 1] for a, b, c in T]
    with Cadical153(bootstrap_with=cl) as s:
        print("010-colourable:", s.solve())
    vec = {i: canon(v) for i, v in enumerate(seed)}
    changed = True
    while changed:
        changed = False
        for v in range(n):
            if v in vec:
                continue
            known = [u for u in vec if adj(u, v)]
            for a, b in itertools.combinations(known, 2):
                c = canon(cross(vec[a], vec[b]))
                if c != (0, 0, 0):
                    vec[v] = c
                    changed = True
                    break
    free = [v for v in range(n) if v not in vec]
    print("determined by cross products:", len(vec), "free:", free)
    bad = []
    for a, b in itertools.combinations(sorted(vec), 2):
        d = sum(x * y for x, y in zip(vec[a], vec[b]))
        if adj(a, b) and d != 0:
            bad.append(("edge but not orthogonal", a, b))
        if not adj(a, b) and d == 0:
            bad.append(("orthogonal but no edge", a, b))
        if vec[a] == vec[b]:
            bad.append(("same ray", a, b))
    print("conflicts among determined rays:", bad[:10], len(bad))
    for v in free:
        print(f"  free vertex {v}: neighbours {[u for u in range(n) if adj(u, v)]}")
    for v in range(len(seed), n):
        if v in vec:
            print(f"  vertex {v}: {vec[v]}")


if __name__ == "__main__":
    main()

"""Ray universes in R^3 and their orthogonality structure (exact integer arithmetic)."""
from itertools import product, permutations
from math import gcd


def canon(v):
    """Primitive representative of the ray through integer vector v (first nonzero entry > 0)."""
    g = 0
    for x in v:
        g = gcd(g, abs(x))
    v = tuple(x // g for x in v)
    for x in v:
        if x:
            return v if x > 0 else tuple(-y for y in v)
    raise ValueError("zero vector")


def integer_rays(k):
    """All rays spanned by integer vectors with entries in [-k, k]."""
    return sorted({canon(v) for v in product(range(-k, k + 1), repeat=3) if any(v)})


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def structure(rays):
    """Return (edges, triangles) of the orthogonality graph, as index pairs / triples."""
    idx = {r: i for i, r in enumerate(rays)}
    n = len(rays)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if dot(rays[i], rays[j]) == 0]
    tris = set()
    for i, j in edges:
        w = canon(cross(rays[i], rays[j]))
        k = idx.get(w)
        if k is not None:
            tris.add(tuple(sorted((i, j, k))))
    return edges, sorted(tris)


def prune_to_triangles(rays):
    """Drop rays lying in no orthogonal basis of the set (they can always be coloured 0)."""
    _, tris = structure(rays)
    keep = sorted({i for t in tris for i in t})
    return [rays[i] for i in keep]


def signed_perms():
    """The 48 signed permutation matrices, as functions on integer vectors."""
    out = []
    for p in permutations(range(3)):
        for s in product((1, -1), repeat=3):
            out.append(lambda v, p=p, s=s: canon(tuple(s[i] * v[p[i]] for i in range(3))))
    return out

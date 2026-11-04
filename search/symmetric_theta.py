"""RESULTS.md Section 2.5: symmetric 30-event supports have theta(Hbar) < 3.

A GHZ realization with support graph H (positive events, edges = non-orthogonal pairs) gives a feasible point of the
Lovasz SDP for theta(Hbar) with value 3 (B = (1/3) Gram of sqrt(c_x) g_x).  So theta(Hbar) < 3 rules out every
realization, in every dimension, whose support is contained in H.  theta is monotone in the edge set, so only maximal
triangle-free supports need to be checked.

For each group below, acting in the same way on the three 10-event contexts, the script lists all maximal
triangle-free invariant supports (up to symmetry where stated), computes theta(Hbar) and, with --certify, an exact
rational upper bound < 3 (dualcert.exact_upper: a rational dual matrix and an exact PSD check).

  Z10  regular (circulant supports)          D5   regular (Cayley-type supports)
  A5, S5 on the 10 two-subsets of {1..5}     F20 = AGL(1,5) on its 10 cosets of an involution
  Z2^4 : Z5 (even sign changes, 5 blocks of size 2)

usage: python3 search/symmetric_theta.py [--certify]        (about 5 minutes with --certify)
"""
import itertools
import sys
import time

import numpy as np
import cvxpy as cp


# ---------------------------------------------------------------- theta(Hbar)
class ThetaGHZ:
    """max sum(B) s.t. B psd, tr B = 1, B_ij = 0 unless ij is an edge of H (or i = j)."""

    def __init__(self, n):
        self.n = n
        self.B = cp.Variable((n, n), symmetric=True)
        self.mask = cp.Parameter((n, n))
        self.prob = cp.Problem(cp.Maximize(cp.sum(self.B)),
                               [self.B >> 0, cp.trace(self.B) == 1, cp.multiply(self.mask, self.B) == 0])

    def __call__(self, edges):
        mask = np.ones((self.n, self.n))
        for u, v in edges:
            mask[u, v] = mask[v, u] = 0
        np.fill_diagonal(mask, 0)
        self.mask.value = mask
        self.prob.solve(solver="CLARABEL")
        return self.prob.value


# ---------------------------------------------------------------- regular actions (Cayley-type supports)
def cyclic_table(m):
    return np.array([[(a + b) % m for b in range(m)] for a in range(m)])


def dihedral_table(m):
    els = [(i, j) for j in range(2) for i in range(m)]
    idx = {e: k for k, e in enumerate(els)}
    return np.array([[idx[((a[0] + (-1) ** a[1] * b[0]) % m, (a[1] + b[1]) % 2)] for b in els] for a in els])


def cayley_classes(T):
    """Maximal triangle-free triples (S_AB, S_BC, S_AC) of subsets of the group (A_x ~ B_xs for s in S_AB, etc.),
    triangle-free iff S_AB S_BC and S_AC are disjoint.  Classes up to relabelling B and C by right multiplication,
    group automorphisms and permutations of the contexts."""
    g = len(T)
    full = (1 << g) - 1
    inv = [int(np.where(T[a] == 0)[0][0]) for a in range(g)]

    def img(f, S):
        out = 0
        for x in range(g):
            if S >> x & 1:
                out |= 1 << int(f(x))
        return out

    L = np.array([[img(lambda x, a=a: T[a, x], S) for S in range(1 << g)] for a in range(g)], dtype=np.int64)
    R = np.array([[img(lambda x, t=t: T[x, t], S) for S in range(1 << g)] for t in range(g)], dtype=np.int64)
    masks = np.arange(1 << g)
    raw = []
    for SAB in range(1, 1 << g):
        if not SAB & 1:          # relabel B so that the identity lies in S_AB (S_AB is never empty when maximal)
            continue
        elems = [a for a in range(g) if SAB >> a & 1]
        prod = np.bitwise_or.reduce(L[elems], axis=0)
        SAC = full & ~prod
        ok = (SAC != 0) & (masks != 0)
        for s in range(g):
            if not SAB >> s & 1:
                ok &= (L[s] & SAC) != 0
            Rs = int(np.bitwise_or.reduce([1 << int(T[a, s]) for a in elems]))
            ok &= ((masks >> s) & 1).astype(bool) | ((Rs & SAC) != 0)
        raw += [(SAB, int(b), int(SAC[b])) for b in np.nonzero(ok)[0]]
    # automorphisms: brute force over images of a generating pair is enough for g = 10
    auts = []
    gens = [1, g // 2] if g % 2 == 0 else [1]
    for perm in itertools.permutations(range(g), len(gens)):
        # extend homomorphism from generators by BFS; keep if well defined and bijective
        f = {0: 0}
        frontier = [0]
        good = True
        while frontier and good:
            x = frontier.pop()
            for gi, im in zip(gens, perm):
                y, fy = int(T[x, gi]), int(T[f[x], im])
                if y in f:
                    if f[y] != fy:
                        good = False
                        break
                else:
                    f[y] = fy
                    frontier.append(y)
        if good and len(f) == g and len(set(f.values())) == g and \
                all(f[int(T[x, y])] == int(T[f[x], f[y]]) for x in range(g) for y in range(g)):
            auts.append([f[x] for x in range(g)])
    inv_set = lambda S: img(lambda x: inv[x], S)
    seen, classes = set(), []
    for t3 in raw:
        if t3 in seen:
            continue
        orb, fr = set(), [t3]
        while fr:
            x = fr.pop()
            if x in orb:
                continue
            orb.add(x)
            a, b, c = x
            fr += [(inv_set(a), c, b), (c, inv_set(b), a)]
        for a, b, c in orb:
            for f in auts:
                a2, b2, c2 = (img(lambda x: f[x], S) for S in (a, b, c))
                for t in range(g):
                    for u in range(g):
                        seen.add((int(R[t][a2]), int(R[u][L[inv[t]][b2]]), int(R[u][c2])))
        classes.append(t3)

    def edges(t3):
        SAB, SBC, SAC = t3
        E = []
        for x in range(g):
            for s in range(g):
                if SAB >> s & 1: E.append((x, g + int(T[x, s])))
                if SBC >> s & 1: E.append((g + x, 2 * g + int(T[x, s])))
                if SAC >> s & 1: E.append((x, 2 * g + int(T[x, s])))
        return E
    return [edges(t3) for t3 in classes], len(raw)


# ---------------------------------------------------------------- non-regular transitive actions (orbital unions)
def closure(gens, n):
    G, fr = {tuple(range(n))}, [tuple(range(n))]
    while fr:
        a = fr.pop()
        for h in gens:
            b = tuple(h[a[i]] for i in range(n))
            if b not in G:
                G.add(b)
                fr.append(b)
    return G


def orbital_supports(gens, n=10):
    G = closure(gens, n)
    assert len({p[0] for p in G}) == n
    orbs, seen = [], set()
    for x in range(n):
        for y in range(n):
            if (x, y) not in seen:
                o = {(p[x], p[y]) for p in G}
                seen |= o
                orbs.append(o)
    r = len(orbs)

    def adj_of(S3):
        adj = np.zeros((3 * n, 3 * n), bool)
        for S, k, l in zip(S3, (0, 1, 0), (1, 2, 2)):
            for i in S:
                for x, y in orbs[i]:
                    adj[k * n + x, l * n + y] = adj[l * n + y, k * n + x] = True
        return adj

    def tri_free(adj):
        A, B, C = (adj[:n, n:2 * n].astype(int), adj[n:2 * n, 2 * n:].astype(int), adj[:n, 2 * n:].astype(int))
        return not np.any((A @ B) * C)
    subsets = [frozenset(s) for k in range(r + 1) for s in itertools.combinations(range(r), k)]
    good = {S3 for S3 in itertools.product(subsets, repeat=3) if tri_free(adj_of(S3))}
    maximal = [S3 for S3 in good if not any(
        i not in S3[k] and tuple(S3[j] | {i} if j == k else S3[j] for j in range(3)) in good
        for k in range(3) for i in range(r))]
    out = []
    for S3 in maximal:
        adj = adj_of(S3)
        out.append([(u, v) for u in range(3 * n) for v in range(u + 1, 3 * n) if adj[u, v]])
    return out, len(G), r


def two_subsets_gens(alternating):
    P = list(itertools.combinations(range(5), 2))
    idx = {p: i for i, p in enumerate(P)}
    perm = lambda s: tuple(idx[tuple(sorted((s[a], s[b])))] for a, b in P)
    if alternating:
        return [perm((1, 2, 0, 3, 4)), perm((0, 1, 3, 4, 2)), perm((1, 0, 3, 2, 4))]
    return [perm((1, 0, 2, 3, 4)), perm((1, 2, 3, 4, 0))]


def f20_gens():
    els = [(a, b) for a in range(1, 5) for b in range(5)]
    mul = lambda g, h: ((g[0] * h[0]) % 5, (g[0] * h[1] + g[1]) % 5)
    cosets = []
    for x in els:
        c = frozenset(mul(x, h) for h in [(1, 0), (4, 0)])
        if c not in cosets:
            cosets.append(c)
    return [tuple(cosets.index(frozenset(mul(g, x) for x in c)) for c in cosets) for g in [(1, 1), (2, 0)]]


def z2z5_gens():
    pts = [(i, e) for i in range(5) for e in range(2)]
    idx = {p: k for k, p in enumerate(pts)}
    return [tuple(idx[((i + 1) % 5, e)] for i, e in pts), tuple(idx[(i, (e + (i < 2)) % 2)] for i, e in pts)]


def main():
    certify = "--certify" in sys.argv
    if certify:
        import networkx as nx
        from dualcert import exact_upper
    th = ThetaGHZ(30)
    families = []
    t0 = time.time()
    for name, T in (("Z10 regular", cyclic_table(10)), ("D5 regular", dihedral_table(5))):
        sup, nraw = cayley_classes(T)
        families.append((name, sup, f"{nraw} maximal triples with 1 in S_AB, {len(sup)} classes"))
    for name, gens in (("A5 on 2-subsets", two_subsets_gens(True)), ("S5 on 2-subsets", two_subsets_gens(False)),
                       ("F20 on 10 cosets", f20_gens()), ("Z2^4:Z5", z2z5_gens())):
        sup, order, rank = orbital_supports(gens)
        families.append((name, sup, f"order {order}, rank {rank}, {len(sup)} maximal orbital supports"))
    allok = True
    for name, sup, info in families:
        vals = [th(E) for E in sup]
        line = f"{name:18s} {info}: max theta {max(vals):.6f}"
        if certify:
            worst = 0
            for E in sup:
                H = nx.Graph()
                H.add_nodes_from(range(30))
                H.add_edges_from(E)
                bound, _ = exact_upper(H)
                if bound is None or bound >= 3:
                    allok = False
                    worst = None
                    break
                worst = max(worst, bound)
            line += f"; exact upper bounds < 3 for all: {worst is not None}, largest {worst}"
        print(line, f"({time.time() - t0:.0f}s)", flush=True)
    if certify:
        print("ALL CERTIFIED" if allok else "SOME SUPPORT NOT CERTIFIED")


if __name__ == "__main__":
    main()

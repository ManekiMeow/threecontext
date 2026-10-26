# H33: a three-context GHZ graph on 33 vertices with a real realization in dimension 20 (RESULTS.md 1.0, 2).
# Contexts: A = the 9 points of the affine plane AG(2,3), B and C = two copies of its 12 lines.
# Edges of H (non-orthogonal pairs): point a ~ line b (b in B or C) iff a lies on b;  b in B ~ c in C iff the lines
# are parallel and distinct.  Gram matrix: 1/sqrt(12) on incidences, 1/2 on parallel pairs, 0 on non-edges.
# Exact checks in rational arithmetic, using K = D^-1 G D with D = diag(1 on A, sqrt3 on B, C), which is rational and
# similar to G:  K (K - 3/2)(K - 2)(K - 3) = 0, so G is PSD with spectrum in {0, 3/2, 2, 3};  rank 20;  the state
# psi = (1/3) sum_A g_a = (1/sqrt12) sum_B g_b = (1/sqrt12) sum_C g_c, with probabilities 1/9 (A) and 1/12 (B, C).
import itertools, numpy as np, networkx as nx, sympy as sp
from pysat.solvers import Cadical153

pts = list(itertools.product(range(3), repeat=2))
lines = []
for v in [(1, 0), (0, 1), (1, 1), (1, 2)]:
    seen = set()
    for p in pts:
        L = frozenset(((p[0] + t * v[0]) % 3, (p[1] + t * v[1]) % 3) for t in range(3))
        if L not in seen: seen.add(L); lines.append((v, L))
assert len(lines) == 12
A = [("A", p) for p in pts]; Bv = [("B", i) for i in range(12)]; Cv = [("C", i) for i in range(12)]
V = A + Bv + Cv; idx = {v: k for k, v in enumerate(V)}; n = len(V)
H = nx.Graph(); H.add_nodes_from(V)
for i, (v, L) in enumerate(lines):
    for p in L: H.add_edge(("A", p), ("B", i)); H.add_edge(("A", p), ("C", i))
    for j, (w, M) in enumerate(lines):
        if w == v and j != i: H.add_edge(("B", i), ("C", j))
assert sum(1 for c in nx.enumerate_all_cliques(H) if len(c) == 3) == 0
print(f"H33: {n} vertices, {H.number_of_edges()} edges, triangle-free, degrees {sorted(set(d for _, d in H.degree()))}")

half, r12 = sp.Rational(1, 2), sp.sqrt(3) / 6
cls = lambda v: v[0]
K = sp.zeros(n, n); G = np.zeros((n, n))
for u, w in H.edges():
    for s, t in [(u, w), (w, u)]:
        i, j = idx[s], idx[t]
        if "A" in (cls(s), cls(t)):
            G[i, j] = 1 / np.sqrt(12)
            K[i, j] = half if cls(s) == "A" else sp.Rational(1, 6)    # (sqrt3/6) * d_j / d_i
        else:
            G[i, j] = 0.5; K[i, j] = half
G += np.eye(n); K += sp.eye(n)
I = sp.eye(n)
assert K * (K - sp.Rational(3, 2) * I) * (K - 2 * I) * (K - 3 * I) == sp.zeros(n, n)
tr = [sum((K ** k)[i, i] for i in range(n)) for k in (1, 2, 3)]
m = sp.symbols("m0 m1 m2 m3")
mult = sp.solve([sum(m) - n] + [m[1] * sp.Rational(3, 2) ** k + m[2] * 2 ** k + m[3] * 3 ** k - tr[k - 1]
                                for k in (1, 2, 3)], m)
print("exact: K(K-3/2)(K-2)(K-3) = 0; multiplicities of 0, 3/2, 2, 3:", [mult[x] for x in m], "rank", K.rank())
yA = sp.Matrix([sp.Rational(1, 3) if cls(v) == "A" else 0 for v in V])
for X in "BC":
    z = yA - sp.Matrix([sp.Rational(1, 6) if cls(v) == X else 0 for v in V])  # D^-1 (y_A - y_X)
    assert K * z == sp.zeros(n, 1)
print("exact: (1/3) sum_A g_a = (1/sqrt12) sum_B g_b = (1/sqrt12) sum_C g_c =: psi, |psi| = 1")
prob = (G @ np.array([1 / 3 if cls(v) == "A" else 0 for v in V])) ** 2
for X in "ABC":
    pr = prob[[idx[v] for v in V if cls(v) == X]]
    print(f"  context {X}: {len(pr)} events, probabilities {sorted({round(float(p), 12) for p in pr})}, sum {pr.sum():.12f}")
ev, U = np.linalg.eigh(G); keep = ev > 1e-9
vec = U[:, keep] * np.sqrt(ev[keep])           # rows: real vectors in R^20 with Gram matrix G
psi = vec[[idx[v] for v in A]].sum(0) / 3
print(f"realization in R^{vec.shape[1]}: max |Gram - G| = {np.abs(vec @ vec.T - G).max():.1e}; "
      f"F spectrum = {sorted({round(float(e), 9) for e in ev[keep]})}")
for X in "ABC":
    S = vec[[idx[v] for v in V if cls(v) == X]]
    assert np.allclose(S @ S.T, np.eye(len(S))) and np.isclose(((S @ psi) ** 2).sum(), 1)

# every proper 3-colouring of H33: each class has probability 1 (theta argument), so its sizes (k points, l lines)
# satisfy k/9 + l/12 = 1, and every vertex sees both other colours (rigidity)
s = Cadical153(); x = lambda v, c: 3 * idx[v] + c + 1
for v in V:
    s.add_clause([x(v, c) for c in range(3)])
    for c, d in itertools.combinations(range(3), 2): s.add_clause([-x(v, c), -x(v, d)])
for u, w in H.edges():
    for c in range(3): s.add_clause([-x(u, c), -x(w, c)])
s.add_clause([x(V[0], 0)]); count = 0; shapes = set()
while s.solve():
    mdl = set(l for l in s.get_model() if l > 0); col = {v: next(c for c in range(3) if x(v, c) in mdl) for v in V}
    count += 1
    for c in range(3):
        k = sum(1 for v in A if col[v] == c); l = sum(1 for v in V if col[v] == c) - k
        assert sp.Rational(k, 9) + sp.Rational(l, 12) == 1; shapes.add((k, l))
    assert all({col[u] for u in H[v]} == set(range(3)) - {col[v]} for v in V)
    s.add_clause([-x(v, col[v]) for v in V])
print(f"proper 3-colourings (vertex {V[0]} fixed to colour 0): {count}; class shapes (points, lines): {sorted(shapes)};"
      " every vertex sees both other colours in each")
nx.write_graph6(nx.relabel_nodes(H, idx), "h33.g6", header=False)

# Theorem 4, step 6 (uniqueness of the equality case): with B = the lines of AG(2,3) and all probabilities 1/12, a
# C-event whose A-neighbourhood is the triple c needs  sum_{b ∋ a, b ∩ c = ∅} <g_c|g_b> = 1/2  for every point a ∉ c.
# Exact check over all 84 triples: solvable iff c is a line, and then the solution is 1/2 on its two parallels.
Ls = [L for _, L in lines]
for c in itertools.combinations(pts, 3):
    c = frozenset(c); D = [b for b in Ls if not b & c]; T = sp.symbols(f"t0:{len(D)}")
    sol = sp.linsolve([sum(T[k] for k, b in enumerate(D) if a in b) - half for a in pts if a not in c], *T)
    assert (sol != sp.EmptySet) == (c in Ls)
    if c in Ls: assert sol == sp.FiniteSet((half, half))
print("equality case: for all 84 point triples c the equations are solvable iff c is a line (unique solution 1/2, 1/2)")

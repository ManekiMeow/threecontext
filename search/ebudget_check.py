"""Lemma 39 and Corollary 40 (RESULTS.md §2.4) checked on the H33 realization.

Lemma 39: for every event x of context k and every context l,
    sum_{y in I_l} |<g_y|E|g_x>|^2 = ||Pi_l E g_x||^2 <= (1 - 9 c_x)/4,
with E = F - (3/2)(1 + |psi><psi|) and explicit entries:
    l != k:  <g_y|E|g_x> = <g_y|g_x>/2 + <g_y|Pi_m|g_x> - (3/2) sqrt(c_x c_y)
    l == k:  <g_y|E|g_x> = <g_y|Pi_l'|g_x> + <g_y|Pi_m|g_x> - (3/2) sqrt(c_x c_y)   (y != x)
Corollary 40 (covering): c_x * c(B_0) <= 1/9 - c_x and c_x * c(D_2) <= 1/9 - c_x.
"""
import itertools

import numpy as np

pts = [(x, y) for x in range(3) for y in range(3)]
lines = []
for dx, dy in [(1, 0), (0, 1), (1, 1), (1, 2)]:
    nx, ny = (-dy) % 3, dx % 3
    for off in range(3):
        lines.append(frozenset(p for p in pts if (nx * p[0] + ny * p[1]) % 3 == off))
V = [("A", p) for p in pts] + [("B", l) for l in lines] + [("C", l) for l in lines]
n = len(V)
G = np.eye(n)
for i, j in itertools.combinations(range(n), 2):
    (ki, xi), (kj, xj) = V[i], V[j]
    if {ki, kj} in ({"A", "B"}, {"A", "C"}):
        p, l = (xi, xj) if ki == "A" else (xj, xi)
        if p in l:
            G[i, j] = G[j, i] = 1 / np.sqrt(12)
    if {ki, kj} == {"B", "C"} and xi != xj and not (xi & xj):
        G[i, j] = G[j, i] = 0.5
w, U = np.linalg.eigh(G)
g = (U[:, w > 1e-9] * np.sqrt(w[w > 1e-9])).T                 # columns: event vectors, d = 20
d = g.shape[0]
ctx = [list(range(0, 9)), list(range(9, 21)), list(range(21, 33))]
psi = g[:, :9].sum(1) / 3
for i in range(n):                                            # phases <g|psi> > 0
    if g[:, i] @ psi < 0:
        g[:, i] *= -1
Pi = [g[:, I] @ g[:, I].T for I in ctx]
E = sum(Pi) - 1.5 * (np.eye(d) + np.outer(psi, psi))
c = (g.T @ psi) ** 2
adj = np.abs(g.T @ g) > 1e-9
worst, err, viol = -1.0, 0.0, 0
for k in range(3):
    for x in ctx[k]:
        Eg = E @ g[:, x]
        budget = (1 - 9 * c[x]) / 4
        assert abs(Eg @ Eg - budget) < 1e-12                 # Lemma 3
        for l in range(3):
            worst = max(worst, sum((g[:, y] @ Eg) ** 2 for y in ctx[l]) - budget)
            others = [j for j in range(3) if j != l] if l == k else [3 - k - l]
            for y in ctx[l]:
                if y == x:
                    continue
                f = sum(g[:, y] @ Pi[j] @ g[:, x] for j in others) - 1.5 * np.sqrt(c[x] * c[y])
                if l != k:
                    f += g[:, y] @ g[:, x] / 2
                err = max(err, abs(f - g[:, y] @ Eg))
            if l != k:
                m = 3 - k - l
                B0 = [y for y in ctx[l] if not adj[x, y] and not any(adj[x, z] and adj[z, y] for z in ctx[m])]
                viol += c[x] * c[B0].sum() > 1 / 9 - c[x] + 1e-12
        D2 = [y for y in ctx[k] if y != x and not any(adj[x, z] and adj[z, y] for z in range(n) if z not in ctx[k])]
        viol += c[x] * c[D2].sum() > 1 / 9 - c[x] + 1e-12
print(f"H33: max_x,l [sum_(y in I_l) |<g_y|E g_x>|^2 - (1 - 9c_x)/4] = {worst:.1e} (Lemma 39: <= 0)")
print(f"     entry formulas, max error {err:.1e}; covering violations (Corollary 40): {viol}")

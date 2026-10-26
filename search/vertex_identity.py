# Numerical sanity check of the per-vertex identities behind RESULTS.md 2 (Lemma 3, Theorem 4) on known realizations.
# For a realization (unit vectors g_i, contexts I_1..I_3, state psi), F = sum_i g_i g_i^dag and E = F - 3/2 - 3/2 psi psi^dag:
#   (i)  ||E g_i||^2 = (1 - 9 c_i)/4 for every event i (c_i = |<psi|g_i>|^2), hence c_i <= 1/9;
#   (ii) c_i (1/alpha + 1/beta + 1/rho) <= 1, alpha/beta = probability of the neighbourhood of i in the other two
#        contexts, rho = 1 - alpha - beta.
# Realizations are rebuilt from the exact certificates (B = Y + diag(c) = Gram matrix of w_i = <g_i|psi> g_i).
import json, numpy as np
from verify_all import tricayley, tricayley_z2z6

def realization(V, E, Y, c):
    idx = {v: i for i, v in enumerate(V)}; B = np.diag([float(c[v]) for v in V])
    for (u, v), y in Y.items(): B[idx[u], idx[v]] = B[idx[v], idx[u]] = float(y)
    ev, U = np.linalg.eigh(B); k = ev > 1e-10; W = U[:, k] * np.sqrt(ev[k])      # rows w_i
    cc = (W ** 2).sum(1); g = W / np.sqrt(cc)[:, None]
    psi = W[[idx[v] for v in V if v[0] == 0]].sum(0); psi /= np.linalg.norm(psi)
    return idx, g, psi

def check(name, V, E, Y, c):
    idx, g, psi = realization(V, E, Y, c); d = g.shape[1]
    F = g.T @ g; Em = F - 1.5 * np.eye(d) - 1.5 * np.outer(psi, psi)
    p = (g @ psi) ** 2
    lhs = ((g @ Em.T) ** 2).sum(1); rhs = (1 - 9 * p) / 4
    adj = {v: set() for v in V}
    for u, v in E: adj[u].add(v); adj[v].add(u)
    worst = 0
    for v in V:
        a, b = [sum(p[idx[u]] for u in adj[v] if u[0] == k) for k in range(3) if k != v[0]]
        worst = max(worst, p[idx[v]] * (1 / a + 1 / b + 1 / (1 - a - b)))
    ctx = [round(float(sum(p[idx[v]] for v in V if v[0] == k)), 9) for k in range(3)]
    print(f"{name}: dim {d}, context sums {ctx}, max c_i = {p.max():.4f} (<= 1/9 = 0.1111), "
          f"max |(i) lhs - rhs| = {abs(lhs - rhs).max():.1e}, max c(1/a+1/b+1/r) = {worst:.4f} (<= 1)")

from fractions import Fraction as Fr
S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
V, E = tricayley_z2z6(S36); check("H36", V, E, {e: Fr(1, 48) for e in E}, {v: Fr(1, 12) for v in V})
V, E = tricayley(13, ((0, 1, 3, 9), (0, 1, 10), (1, 6, 8))); w = [Fr(1, 4), Fr(1, 3), Fr(1, 3)]
check("H39", V, E, {(u, v): w[u[0]] * Fr(1, 13) for u, v in E}, {v: Fr(1, 13) for v in V})
dj = json.load(open("shr_q13_1_0139_014_1257_cert.json"))
V, E = tricayley(13, tuple(tuple(s) for s in dj["S"]), removed={(0, 5), (2, 11)})
Y = {}
for a, b, y in dj["edges"]:
    a, b = tuple(a), tuple(b); Y[(a, b) if (a, b) in E else (b, a)] = Fr(y)
check("H37", V, E, Y, {tuple(v): Fr(x) for v, x in zip(dj["nodes"], dj["c"])})
import networkx as nx, itertools
H = nx.read_graph6("h33.g6"); lab = lambda i: (0 if i < 9 else 1 if i < 21 else 2, i)
V = [lab(i) for i in range(33)]; E = {(lab(u), lab(v)) for u, v in H.edges()}
cv = {v: Fr(1, 9) if v[0] == 0 else Fr(1, 12) for v in V}
# Y_uv = sqrt(c_u c_v) <g_u|g_v>: incidences sqrt(1/108)*sqrt(1/12)... use floats via the Gram matrix directly
Yf = {(u, v): (np.sqrt(float(cv[u] * cv[v])) * (1 / np.sqrt(12) if 0 in (u[0], v[0]) else 0.5)) for u, v in E}
check("H33", V, E, Yf, cv)

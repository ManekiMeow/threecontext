import sys, numpy as np
sys.path.insert(0, "../search")
from ghzstates import realisation, completions, symmetry_unitary, triangle_tensors
from ghzmilp import extra_ghz_state
from verify_all import tricayley_z2z6, tricayley
from perkel import perkel
def graph(name):
    if name == "H36":
        S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
        V, E = tricayley_z2z6(S36); idx = {v: i for i, v in enumerate(V)}; P = np.eye(36)
        for u, v in E: P[idx[u], idx[v]] = P[idx[v], idx[u]] = 0.25
        shift = [idx[(k, (g[0], (g[1] + 1) % 6))] for (k, g) in V]
    elif name == "H39":
        V, E = tricayley(13, ((0, 1, 3, 9), (0, 1, 10), (1, 6, 8))); idx = {v: i for i, v in enumerate(V)}
        P = np.eye(39); w = [0.25, 1 / 3, 1 / 3]
        for u, v in E: P[idx[u], idx[v]] = P[idx[v], idx[u]] = w[u[0]]
        shift = [idx[(k, (x + 1) % 13)] for (k, x) in V]
    else:
        G = perkel(); V = sorted(G); idx = {v: i for i, v in enumerate(V)}; P = np.eye(57)
        for u, v in G.edges(): P[idx[u], idx[v]] = P[idx[v], idx[u]] = 1 / 3
        shift = [idx[(k, (x + 1) % 19)] for (k, x) in V]
    return P, [v[0] for v in V], shift
if __name__ == "__main__":
  name, mode = sys.argv[1], sys.argv[2]
  P, cls, shift = graph(name)
  X, psi0 = realisation(P, cls); d = X.shape[1]
  B = completions(X, cls, mode, np.random.default_rng(0), symmetry_unitary(X, shift))
  E12, E23, E13 = triangle_tensors(B)
  print(f"{name} d={d} completion={mode}: rays {3*d}, triangles {int((E12 @ E23 * E13).sum())}", flush=True)
  found, stats = extra_ghz_state(B, [psi0], time_limit=900)
  names = {0: "feasible", 1: "time limit", 2: "INFEASIBLE"}
  print(f"{name} completion={mode}: GHZ state orthogonal to psi0 exists: {found is not None};"
        f" per-r* status counts: { {names.get(s, s): stats.count(s) for s in set(stats)} } over {len(stats)} of {d}", flush=True)

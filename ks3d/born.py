"""Degree-2 ('Born') realizability test for orthogonality graphs in C^3 (RESULTS.md §7.3).

For rays v_i put Y_ij = |<v_i|v_j>|^2 and X = (3Y - J)/2. Then X is the Gram matrix of the unit
vectors (traceless parts of the projectors, normalised), so
    X psd,  X_ii = 1,  X_ij = -1/2 on edges,  X_ij >= -1/2,  X_ia + X_ib <= 1/2 (edge ab, any i),
the last because two orthonormal rays carry total probability <= 1 in any state v_i. These are the
only constraints used: everything here is a semidefinite (degree-2 in Y) test.

Propagation: if every feasible X has X_ij = -1/2 for a non-edge ij, then v_i is orthogonal to v_j
in every embedding, so the edge ij can be added. If the enlarged graph contains a 4-cycle the
graph is not embeddable (two distinct rays have exactly one common orthogonal ray in C^3). If
some X_ij is forced to 1, two rays coincide, which is also excluded.
"""
import itertools
import json
import os
import sys
import warnings

import cvxpy as cp
import networkx as nx
import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)
HERE = os.path.dirname(os.path.abspath(__file__))
EPS = 1e-5


def born_problem(n, E):
    X = cp.Variable((n, n), symmetric=True)
    A = np.zeros((n, n), dtype=bool)
    for i, j in E:
        A[i, j] = A[j, i] = True
    cons = [X >> 0, cp.diag(X) == 1]
    for i in range(n):
        for j in range(i + 1, n):
            cons.append(X[i, j] == -0.5 if A[i, j] else X[i, j] >= -0.5)
    for a, b in E:
        for i in range(n):
            if i not in (a, b):
                cons.append(X[i, a] + X[i, b] <= 0.5)
    return X, cons, A


def extreme(n, E, pair, sense, T=None):
    X, cons, _ = born_problem(n, E) if T is None else born3_problem(n, E, T)
    i, j = pair
    p = cp.Problem(cp.Maximize(X[i, j]) if sense > 0 else cp.Minimize(X[i, j]), cons)
    p.solve(solver=cp.CLARABEL)
    if p.status not in ("optimal", "optimal_inaccurate"):
        return None
    return p.value


def born3_problem(n, E, T):
    """Degree-3 strengthening. With Z_ijk = Re <i|j><j|k><k|i> (fully symmetric), for each ray i
    the Gram matrix of |i> and the vectors Pi_j|i> is [[1, Y_i.], [Y_i., W_i]] with
    Re W_i[j,k] = Z_ijk and W_i[j,j] = Y_ij: it must be psd, vanish on edges jk, and satisfy
    sum_{t in T} Z_ijt = Y_ij for every basis T. (This is the degree-2 moment matrix of the state
    v_i; the symmetry of Z couples the states.)"""
    X, cons, A = born_problem(n, E)
    Y = (2 * X + 1) / 3
    trip = {S: k for k, S in enumerate(itertools.combinations(range(n), 3))}
    Z = cp.Variable(len(trip))

    def z(i, j, k):
        if j == k:
            return Y[i, j]
        if i == j:
            return Y[i, k]
        if i == k:
            return Y[i, j]
        return Z[trip[tuple(sorted((i, j, k)))]]

    for i in range(n):
        others = [j for j in range(n) if j != i]
        M = cp.Variable((n, n), symmetric=True)
        cons += [M >> 0, M[0, 0] == 1]
        for a, j in enumerate(others, 1):
            cons += [M[0, a] == Y[i, j], M[a, a] == Y[i, j]]
            for b in range(a + 1, n):
                cons.append(M[a, b] == z(i, j, others[b - 1]))
    for (i, j, k), idx in trip.items():
        if A[i, j] or A[j, k] or A[i, k]:
            cons.append(Z[idx] == 0)
    for t3 in T:
        for i in range(n):
            for j in range(n):
                if i < j:
                    cons.append(sum(z(i, j, t) for t in t3) == Y[i, j])
    return X, cons, A


def has_square(G):
    """A 4-cycle as a subgraph (not necessarily induced)."""
    for u, v in itertools.combinations(G.nodes(), 2):
        if len(set(G[u]) & set(G[v])) >= 2:
            return (u, v, *sorted(set(G[u]) & set(G[v]))[:2])
    return None


def propagate(G, degree=2):
    """Returns (verdict, log). verdict: 'refuted' or 'passes'."""
    G = nx.convert_node_labels_to_integers(G)
    n = G.number_of_nodes()
    log = []
    while True:
        E = sorted(tuple(sorted(e)) for e in G.edges())
        T = [t for t in itertools.combinations(range(n), 3)
             if G.has_edge(t[0], t[1]) and G.has_edge(t[0], t[2]) and G.has_edge(t[1], t[2])]
        T = None if degree == 2 else T
        X, cons, A = born_problem(n, E) if T is None else born3_problem(n, E, T)
        pairs = [(i, j) for i, j in itertools.combinations(range(n), 2) if not A[i, j]]
        # one solve that pushes every non-edge away from -1/2 and from 1 (capped); a pair that ends
        # up away from both cannot be forced, so only the others need individual solves
        obj = sum(cp.minimum(X[i, j] + 0.5, 0.05) + cp.minimum(1 - X[i, j], 0.05) for i, j in pairs)
        p = cp.Problem(cp.Maximize(obj), cons)
        p.solve(solver=cp.CLARABEL)
        if p.status not in ("optimal", "optimal_inaccurate"):
            log.append("Born SDP infeasible")
            return "refuted", log
        Xv = X.value
        suspects = [(i, j) for i, j in pairs if Xv[i, j] < -0.5 + 1e-3 or Xv[i, j] > 1 - 1e-3]
        new = []
        for i, j in suspects:
            hi = extreme(n, E, (i, j), +1, T)
            if hi is not None and hi < -0.5 + EPS:
                new.append((i, j))
            lo = extreme(n, E, (i, j), -1, T)
            if lo is not None and lo > 1 - EPS:
                log.append(f"rays {i},{j} forced equal")
                return "refuted", log
        if not new:
            return "passes", log
        log.append(f"forced orthogonal: {new}")
        G.add_edges_from(new)
        sq = has_square(G)
        if sq:
            log.append(f"4-cycle {sq} after adding forced edges")
            return "refuted", log


def _deg3(s):
    verdict, log = propagate(nx.from_graph6_bytes(s.encode()), degree=3)
    return s, verdict, log


if __name__ == "__main__":
    import spectral
    from multiprocessing import Pool
    # mugs.json: the 17 minimal unembeddable graphs of Li-Bright-Ganesh (orders 10-12), copied from
    # github.com/BrianLi009/SAT-nauty, log_archive/min_nonembed_graph_10-12.txt. The 13 listed in
    # that solver's cadical-rcl/src/unembeddable_graphs.h are marked with *.
    in_solver = {"ICOedPKL?", "I?qa``eeO", "J@`ICC[Gv_?", "JsOH?`JH_i_", "J`hSA?hD_F_", "JRPCSGoAgJ_",
                 "K@QCAGdA_c~?", "KAICH@@DOT^?", "Kko_GC`C?b`q", "K`?LAQOP@KbK", "K`?KA@se_YCX",
                 "Ks?I@?XDaTCi", "KqCaC?[I_J?Z"}
    mugs = json.load(open(os.path.join(HERE, "mugs.json")))
    code = lambda s: f"`{s}`" if "``" in s or "`" not in s else f"`` {s} ``"   # markdown code span
    if "--degree3" in sys.argv:
        with Pool(4) as pool:
            res = pool.map(_deg3, [s for gs in mugs.values() for s in gs], chunksize=1)
        print("| graph6 | degree-3 Born propagation |")
        print("|---|---|")
        for s, verdict, log in res:
            print(f"| {code(s)}{'*' if s in in_solver else ''} | {verdict}: {'; '.join(log)} |", flush=True)
        sys.exit()
    print("| graph6 | n | edges | triangles | Hoffman | inertia | p=2 | theta(Gc) | Born slack |"
          " Born propagation |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for order, gs in mugs.items():
        for s in gs:
            G = nx.from_graph6_bytes(s.encode())
            g = spectral.from_graph(G.number_of_nodes(), list(G.edges()))
            sb = spectral.spectral_bounds(g.A.astype(float))
            verdict, log = propagate(G)
            star = "*" if s in in_solver else ""
            print(f"| {code(s)}{star} | {g.n} | {len(g.E)} | {g.m} | {sb['hoffman']:.3f} | {sb['energy0']:.3f} |"
                  f" {sb['energy2']:.3f} | {spectral.theta_complement(g):.3f} |"
                  f" {spectral.born_sdp_subgraph(g):.3f} | {verdict}: {'; '.join(log)} |", flush=True)

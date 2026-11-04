"""RESULTS.md Section 2.5: the E^2-scalar case as a weighted cubic graph problem.

A uniform 30-event GHZ realization with E^2 scalar is the same as a Hermitian weighting Y of a triangle-free
(10,10,10) support with zero diagonal blocks, row sums 1 into each other context, and

    Y^3 = (7/9) Y + (2/9)(J - I).

This script
  1. calibrates the least-squares solver on H36, whose weights Y = A/4 (spectrum {2, 1/2, 0, -1}) satisfy the
     analogous cubic (Y - 1/2) Y (Y + 1) = J/4; it must recover them from random starts on H36's support;
  2. runs the solver on supports from the SAT stream of event30_uniform_search.build (triangle-free,
     Proposition 37 degree pairs) and reports the smallest residual found.
Heuristic: a large residual on a support is evidence, not proof, that the support has no such weighting.

usage: python3 search/e2scalar_cubic.py [SECONDS] [SEED]      (default 300 s)
"""
import sys
import time

import numpy as np
from scipy.optimize import least_squares


def calibrate_h36():
    from verify_all import tricayley_z2z6
    S = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
    V, E = tricayley_z2z6(S)
    idx = {v: i for i, v in enumerate(V)}
    N = len(V)
    ctx = np.array([V[i][0] for i in range(N)])
    edges = np.array([(idx[u], idx[v]) for u, v in E])
    J = np.ones((N, N))
    kappa = (2 - 0.5) * 2 * 3 / N                    # (Y - 1/2) Y (Y + 1) on the all-ones vector, divided by N
    iu = np.triu_indices(N)

    def res(v):
        Y = np.zeros((N, N))
        Y[edges[:, 0], edges[:, 1]] = v
        Y[edges[:, 1], edges[:, 0]] = v
        C = (Y - 0.5 * np.eye(N)) @ Y @ (Y + np.eye(N)) - kappa * J
        rs = [Y[x, ctx == k].sum() - 1 for x in range(N) for k in range(3) if k != ctx[x]]
        return np.concatenate([C[iu], rs])
    rng = np.random.default_rng(0)
    out = []
    for _ in range(3):
        r = least_squares(res, 0.25 * (1 + 0.3 * rng.normal(size=len(edges))), method="lm", max_nfev=3000,
                          xtol=1e-14, ftol=1e-14)
        out.append((r.cost, np.abs(r.x - 0.25).max()))
    print("H36 control (cubic (Y-1/2)Y(Y+1) = kappa J): residual costs", [f"{c:.1e}" for c, _ in out],
          "| max |w - 1/4|", [f"{d:.1e}" for _, d in out])


def run_stream(seconds, seed):
    import networkx as nx
    from event30_uniform_search import build
    N = 30
    ctx = np.array([0] * 10 + [1] * 10 + [2] * 10)
    J = np.ones((N, N))
    iu = np.triu_indices(N)
    rng = np.random.default_rng(seed)
    solver, edge, pairs, parts = build(seed)
    t0, n, best = time.time(), 0, np.inf
    while time.time() - t0 < seconds:
        solver.conf_budget(5000)
        sat = solver.solve_limited()
        if sat is None:
            continue
        if not sat:
            break
        mdl = {x for x in solver.get_model() if x > 0}
        es = np.array([(u, v) for u, v in pairs if edge(u, v) in mdl])
        solver.add_clause([-edge(u, v) if edge(u, v) in mdl else edge(u, v) for u, v in pairs])
        deg = np.zeros((N, 3))
        for u, v in es:
            deg[u, ctx[v]] += 1
            deg[v, ctx[u]] += 1

        def res(w):
            Y = np.zeros((N, N))
            Y[es[:, 0], es[:, 1]] = w
            Y[es[:, 1], es[:, 0]] = w
            C = Y @ Y @ Y - 7 / 9 * Y - 2 / 9 * (J - np.eye(N))
            rs = [Y[x, ctx == k].sum() - 1 for x in range(N) for k in range(3) if k != ctx[x]]
            return np.concatenate([C[iu], rs])
        cost = min(least_squares(res, np.array([1 / deg[u, ctx[v]] for u, v in es]) * (1 + 0.3 * rng.normal(size=len(es))),
                                 method="lm", max_nfev=3000, xtol=1e-14, ftol=1e-14).cost for _ in range(4))
        n += 1
        best = min(best, cost)
        if cost < 1e-16:
            print("SOLUTION FOUND on support", es.tolist(), flush=True)
            break
    print(f"SAT-stream supports tested: {n}; smallest residual cost {best:.3e} (a weighting needs 0)", flush=True)


if __name__ == "__main__":
    seconds = float(sys.argv[1]) if len(sys.argv) > 1 else 300
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 7
    calibrate_h36()
    run_stream(seconds, seed)

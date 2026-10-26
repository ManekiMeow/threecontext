"""Lovász theta versus triangle packing as a KS screen, on known 3D ray sets.

G = orthogonality graph (edge = orthogonal rays). For a 010-colouring I (an independent set
meeting every triangle) we need alpha(G) >= nu(G), nu = max number of vertex-disjoint triangles.
Hence theta(G) < nu(G)  ==>  alpha(G) < nu(G)  ==>  no 010-colouring (KS), a *sufficient* test.
Also every graph realised by rays in R^3 has theta(G) >= n/3 (maximally mixed state).
"""
import json
import sys

import numpy as np
from cvxopt import matrix, solvers
from ortools.sat.python import cp_model

solvers.options["show_progress"] = False


def graph(vs):
    n = len(vs)
    E = [(i, j) for i in range(n) for j in range(i + 1, n) if np.dot(vs[i], vs[j]) == 0]
    T = [(i, j, k) for i, j in E for k in range(j + 1, n)
         if np.dot(vs[i], vs[k]) == 0 and np.dot(vs[j], vs[k]) == 0]
    return n, E, T


def theta(n, E):
    """theta(G) = max <J, X>, tr X = 1, X_ij = 0 on edges, X psd (cvxopt SDP, dual form)."""
    # variables: t, and y_ij for edges;  minimise t s.t. t I - J + sum y_ij (E_ij + E_ji) >= 0
    m = 1 + len(E)
    c = matrix([1.0] + [0.0] * len(E))
    G = np.zeros((n * n, m))
    G[:, 0] = -np.eye(n).ravel()
    for k, (i, j) in enumerate(E):
        A = np.zeros((n, n)); A[i, j] = A[j, i] = -1.0
        G[:, k + 1] = A.ravel()
    h = matrix(-np.ones((n, n)).ravel())
    sol = solvers.sdp(c, Gs=[matrix(G)], hs=[matrix(-np.ones((n, n)))])
    return sol["primal objective"]


def cp_max(n, E, T, what):
    m = cp_model.CpModel()
    if what == "alpha":
        x = [m.NewBoolVar("") for _ in range(n)]
        for i, j in E:
            m.AddBoolOr([x[i].Not(), x[j].Not()])
        m.Maximize(sum(x))
    else:
        x = [m.NewBoolVar("") for _ in T]
        for v in range(n):
            m.Add(sum(x[k] for k, t in enumerate(T) if v in t) <= 1)
        m.Maximize(sum(x))
    s = cp_model.CpSolver(); s.Solve(m)
    return int(s.ObjectiveValue())


def colourable(n, E, T):
    from pysat.solvers import Cadical153
    with Cadical153(bootstrap_with=[[-(i + 1), -(j + 1)] for i, j in E] +
                    [[a + 1, b + 1, c + 1] for a, b, c in T]) as s:
        return s.solve()


def report(name, vs):
    n, E, T = graph(vs)
    th = theta(n, E)
    a, nu = cp_max(n, E, T, "alpha"), cp_max(n, E, T, "nu")
    print(f"{name:28s} n={n:2d} bases={len(T):2d} alpha={a:2d} theta={th:6.3f} n/3={n/3:5.2f} "
          f"nu={nu:2d} theta<nu={th < nu - 1e-6!s:5s} KS={not colourable(n, E, T)}", flush=True)


if __name__ == "__main__":
    sets = json.load(open(sys.argv[1]))
    for name, vs in sets.items():
        vs = [np.array(v) for v in vs]
        report(name, vs)
        if name == "min_height2":
            for i in range(0, len(vs), 10):
                report(f"  {name} minus ray {i}", vs[:i] + vs[i + 1:])

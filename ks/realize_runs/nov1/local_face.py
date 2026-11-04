"""Does the local theta optimum pin the moduli |<u|c>|^2 of a realization?

For a ray u, B_cc' = <u|c><c|c'><c'|u>/3 (c, c' in N(u)) is an optimal point of the local Lovasz SDP (value 3).
Solve the dual for a maximal-rank slack Z; the optimal face is {B = V W V^T : W >= 0, linear constraints}, V = ker Z.
Report: face dimension (real), and the dimension of its image under B -> diag(B) (= moduli/3).
usage: diagface.py pattern.json|.npy [u ...]
"""
import os
import sys
import numpy as np
import cvxpy as cp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))  # ks/
from realize_search import load_pattern

X = load_pattern(sys.argv[1])
n = len(X)
d = n // 4
us = [int(a) for a in sys.argv[2:]] or [0, d, 2 * d, 3 * d]
for u in us:
    N = np.where(X[u])[0]
    m = len(N)
    A = X[np.ix_(N, N)]
    orth = [(i, j) for i in range(m) for j in range(i + 1, m) if not A[i, j]]
    # dual: min lam s.t. lam I - J - sum y_ij (E_ij + E_ji) >= 0
    lam = cp.Variable()
    y = cp.Variable(len(orth))
    I_, J_ = zip(*orth)
    Ymat = np.zeros((m, m, len(orth)))
    E = np.zeros((m * m, len(orth)))
    for k, (i, j) in enumerate(orth):
        E[i * m + j, k] = 1
        E[j * m + i, k] = 1
    Z = lam * np.eye(m) - np.ones((m, m)) - cp.reshape(E @ y, (m, m), order="C")
    prob = cp.Problem(cp.Minimize(lam), [Z >> 0])
    prob.solve(solver="CLARABEL")
    Zv = lam.value * np.eye(m) - np.ones((m, m)) - (E @ y.value).reshape(m, m)
    w, Q = np.linalg.eigh((Zv + Zv.T) / 2)
    rkZ = int((w > 1e-6 * max(1, w.max())).sum())
    V = Q[:, w <= 1e-6 * max(1, w.max())]                  # kernel of Z
    r = V.shape[1]
    # primal too (to check strict complementarity)
    B = cp.Variable((m, m), symmetric=True)
    cons = [B >> 0, cp.trace(B) == 1] + [B[i, j] == 0 for (i, j) in orth]
    pv = cp.Problem(cp.Maximize(cp.sum(B)), cons)
    pv.solve(solver="CLARABEL")
    wb = np.linalg.eigvalsh(B.value)
    rkB = int((wb > 1e-6 * wb.max()).sum())
    # face directions: W sym r x r with V W V^T zero on orth pairs, trace 0, sum 0
    iu = np.triu_indices(r)
    def lift(wvec):
        W = np.zeros((r, r)); W[iu] = wvec; W = W + W.T - np.diag(np.diag(W)); return V @ W @ V.T
    nW = len(iu[0])
    cols = []
    for k in range(nW):
        e = np.zeros(nW); e[k] = 1
        cols.append(lift(e))
    L = np.array(cols)                                        # nW x m x m
    rows = [L[:, i, j] for (i, j) in orth] + [np.trace(L, axis1=1, axis2=2), L.sum(axis=(1, 2))]
    C = np.array(rows)                                        # constraints x nW
    sv = np.linalg.svd(C, compute_uv=False)
    tol = 1e-8 * sv.max()
    rkC = int((sv > tol).sum())
    _, s2, Vt = np.linalg.svd(C)
    null = Vt[rkC:].T                                         # nW x dimface
    Dg = np.array([L[:, i, i] for i in range(m)])             # m x nW
    DN = Dg @ null
    sd = np.linalg.svd(DN, compute_uv=False) if DN.size else np.array([0.0])
    rkD = int((sd > 1e-8 * max(1.0, sd.max() if sd.size else 1)).sum())
    print(f"u = {u} (basis {u // d}): |N(u)| = {m}, orth pairs {len(orth)}; theta = {pv.value:.6f} / dual {lam.value:.6f};"
          f" rank B = {rkB}, rank Z = {rkZ}, ker Z = {r}; face dim (real) = {null.shape[1]}; "
          f"diag-image dim = {rkD}  (singular values of diag map: {np.round(sd[:6], 6)})", flush=True)
    diagB = np.diag(B.value) * 3
    print(f"   moduli at the IPM optimum: min {diagB.min():.5f}, max {diagB.max():.5f}, 1/h = {1/ (m/3):.5f}", flush=True)

"""Level N = 2 of the rank-constrained hierarchy of Yu, Simnacher, Nguyen and Gühne (PRX Quantum 3, 010340 (2022);
references/) for orthonormal representations: does a graph H on n vertices have unit vectors in C^k with
<a_i|a_j> = 0 whenever i != j are non-adjacent? (§III B of the paper; for a four-context skeleton, H is the
non-orthogonality graph and k = d.) Level 1 is the Lovász bound k >= theta(H).

Complex case, with the U(k) symmetry of the auxiliary space used as in their Eq. (25): Phi_AB = Phi_I (x) 1 +
Phi_V (x) V on (C^n (x) C^k)^{(x)2}, with P = Phi_I Hermitian and swap-invariant on C^n (x) C^n and Phi_V = V P.
M = k P + V P plays the role of Gamma (x) Gamma. Constraints:
  P +- VP >= 0 (Phi_AB >= 0);  P^{T_A} >= 0 and P^{T_A} + k (VP)^{T_A} >= 0 (PPT);  k^2 Tr P + k Tr VP = 1;
  M[(i,j),(i,j')] = (1/n) sum_m M[(m,j),(m,j')]    (unit diagonal: equal A-diagonal blocks);
  M[(i,j),(i',j')] = 0 for non-adjacent i != i'     (orthogonality).
We maximise a common margin t (each PSD matrix minus t I stays PSD). Feasible problems have t = 0 (the feasible set
has no interior), infeasible ones t < 0. The feasible set grows with k (Appendix C), so k may be real.

Every matrix is an index permutation of the real and imaginary parts r, s of P (row-major vectors), and Hermitian
PSD constraints use the real embedding [[X, -Y], [Y, X]] >= 0. Solver: SCS.

usage: python3 rankhier.py shrikhande K...      (the Mermin / Shrikhande pattern of RESULTS.md §3.10, n = 16, d = 4)
       python3 rankhier.py c5 K...              (5-cycle: minimum dimension 3, Lovász bound 2.236)
"""
import sys
import time

import cvxpy as cp
import numpy as np
import scipy.sparse as sp


def feasible(A, k, eps=1e-7, max_iters=200000, verbose=False, bases=None):
    """bases: optional list of vertex lists that must be complete orthonormal bases of C^k (then k = their size).
    Each adds the degree-2 identity sum_{m in B} Gamma_im Gamma_mj = Gamma_ij, which is linear in M:
    sum_{m in B} M[(i,m),(m,j)] = (1/n) sum_l M[(i,l),(j,l)]."""
    n = len(A)
    N = n * n
    pair = lambda i, j: i * n + j
    idx = lambda a, b: a * N + b
    sw = np.array([pair(j, i) for i in range(n) for j in range(n)])            # V on C^n (x) C^n
    a_all, b_all = np.divmod(np.arange(N * N), N)
    pV = sw[a_all] * N + b_all                                                  # (VP)[a, b] = P[V a, b]
    ia, ja = np.divmod(a_all, n)
    ib, jb = np.divmod(b_all, n)
    pTA = (ib * n + ja) * N + (ia * n + jb)                                     # P^{T_A}[(i,j),(i',j')] = P[(i',j),(i,j')]
    pVTA = (ja * n + ib) * N + (ia * n + jb)                                    # (VP)^{T_A} = P[(j,i'),(i,j')]
    pT = b_all * N + a_all                                                      # transpose
    pVV = sw[a_all] * N + sw[b_all]                                             # P[V a, V b]

    r, s, t = cp.Variable(N * N), cp.Variable(N * N), cp.Variable()
    mat = lambda v: cp.reshape(v, (N, N), order="C")
    I = np.eye(N)

    def herm_psd(X, Y):
        return cp.bmat([[X, -Y], [Y, X]]) >> 0

    cons = [r == r[pT], s == -s[pT], r == r[pVV], s == s[pVV]]                  # Hermitian, swap-invariant
    cons += [herm_psd(mat(r + r[pV]) - t * I, mat(s + s[pV])),
             herm_psd(mat(r - r[pV]) - t * I, mat(s - s[pV])),
             herm_psd(mat(r[pTA]) - t * I, mat(s[pTA])),
             herm_psd(mat(r[pTA] + k * r[pVTA]) - t * I, mat(s[pTA] + k * s[pVTA]))]
    diag = np.arange(N) * (N + 1)
    cons.append(k * k * cp.sum(r[diag]) + k * cp.sum(r[pV][diag]) == 1)
    Mr, Ms = k * r + r[pV], k * s + s[pV]
    # equal A-diagonal blocks: M[(i,j),(i,j')] - (1/n) sum_m M[(m,j),(m,j')] = 0
    rows, cols, vals = [], [], []
    e = 0
    for i in range(n):
        for j in range(n):
            for jp in range(n):
                rows.append(e); cols.append(idx(pair(i, j), pair(i, jp))); vals.append(1.0)
                for m in range(n):
                    rows.append(e); cols.append(idx(pair(m, j), pair(m, jp))); vals.append(-1.0 / n)
                e += 1
    C = sp.csr_matrix((vals, (rows, cols)), shape=(e, N * N))
    cons += [C @ Mr == 0, C @ Ms == 0]
    zero = [idx(pair(i, j), pair(ip, jp)) for i in range(n) for ip in range(n) if i != ip and not A[i, ip]
            for j in range(n) for jp in range(n)]
    if zero:
        zero = np.array(zero)
        cons += [Mr[zero] == 0, Ms[zero] == 0]
    for B in bases or []:
        # M[(i,j),(i',j')] ~ Gamma_{i'i} Gamma_{j'j}, so sum_{m in B} M[(i,m),(m,j')] ~ (Gamma^T Gamma^T)_{ij'} over B
        rows, cols, vals = [], [], []
        e = 0
        for i in range(n):
            for jp in range(n):
                for m in B:
                    rows.append(e); cols.append(idx(pair(i, m), pair(m, jp))); vals.append(1.0)
                for l in range(n):
                    rows.append(e); cols.append(idx(pair(i, l), pair(jp, l))); vals.append(-1.0 / n)
                e += 1
        Cb = sp.csr_matrix((vals, (rows, cols)), shape=(e, N * N))
        cons += [Cb @ Mr == 0, Cb @ Ms == 0]
    prob = cp.Problem(cp.Maximize(t), cons)
    prob.solve(solver="SCS", eps_abs=eps, eps_rel=eps, max_iters=max_iters, verbose=verbose)
    return prob.status, prob.value


def shrikhande_four_bases():
    """Shrikhande graph on Z_4 x Z_4, connection set {±(1,0), ±(0,1), ±(1,1)}; the four bases are the cosets of
    2Z_4 x 2Z_4 (independent sets). Vertices are ordered basis by basis."""
    conn = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    pts = sorted(((x, y) for x in range(4) for y in range(4)), key=lambda v: (v[0] % 2, v[1] % 2, v))
    n = len(pts)
    A = np.zeros((n, n), dtype=bool)
    for a, u in enumerate(pts):
        for b, v in enumerate(pts):
            A[a, b] = ((v[0] - u[0]) % 4, (v[1] - u[1]) % 4) in conn
    return A


def cycle(n):
    A = np.zeros((n, n), dtype=bool)
    for i in range(n):
        A[i, (i + 1) % n] = A[(i + 1) % n, i] = True
    return A


if __name__ == "__main__":
    graph = sys.argv[1]
    complete = graph.endswith("+bases")
    graph = graph.replace("+bases", "")
    A = shrikhande_four_bases() if graph == "shrikhande" else cycle(5)
    bases = [list(range(4 * b, 4 * b + 4)) for b in range(4)] if complete else None
    for k in map(float, sys.argv[2:]):
        t0 = time.time()
        status, val = feasible(A, k, bases=bases)
        print(f"{graph}{' (complete bases)' if complete else ''}: k = {k:g}: status {status}, margin t = {val:.3e} "
              f"({time.time() - t0:.0f} s)", flush=True)

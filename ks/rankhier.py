"""Level N = 2 of the rank-constrained hierarchy of Yu, Simnacher, Nguyen and Gühne (PRX Quantum 3, 010340 (2022);
references/) for orthonormal representations: does a graph H on n vertices have unit vectors in C^k with
<a_i|a_j> = 0 whenever i != j are non-adjacent? (§III B of the paper; for a four-context skeleton, H is the
non-orthogonality graph and k = d.) Level 1 is the Lovász bound k >= theta(H).

This is Eq. (40) of the paper with SEP relaxed to PPT, complex case, reduced by the U(k) symmetry of the auxiliary
space as in their Eq. (25): Phi_AB = Phi_I (x) 1 + Phi_V (x) V on (C^n (x) C^k)^{(x)2}, with P = Phi_I Hermitian and
swap-invariant on C^n (x) C^n and Phi_V = V P (support on the symmetric subspace). The reduction is exact.
M = k P + V P plays the role of Gamma (x) Gamma. Constraints:
  P +- VP >= 0 (Phi_AB >= 0);  P^{T_A} >= 0 and P^{T_A} + k (VP)^{T_A} >= 0 (PPT);  k^2 Tr P + k Tr VP = n^2;
  M[(i,j),(i,j')] = (1/n) sum_m M[(m,j),(m,j')]    (unit diagonal: equal A-diagonal blocks);
  M[(i,j),(i',j')] = 0 for non-adjacent i != i'     (orthogonality).
By default we maximise a common margin t (each PSD matrix minus t I stays PSD). Feasible problems have t = 0 (the
feasible set has no interior), infeasible ones t < 0. The feasible set grows with k (Appendix C), so k may be real.
The normalisation n^2 (rather than 1) keeps the entries of P of order 1/k^2, so that a small negative margin is not
lost in the solver tolerance.

Every matrix is an index permutation of the real and imaginary parts r, s of P (row-major vectors), and Hermitian
PSD constraints use the real embedding [[X, -Y], [Y, X]] >= 0. Solver: SCS (Clarabel needs > 15 GB at n = 11).

Checks (RESULTS.md §3.13):
  C5: infeasible at k = 2, 2.2; feasible at 2.5, 3, 4.
  the paper's Fig. 3 graph (11 vertices, theta = 4): feasible for k >= 4.01 (8425 SCS iterations at 4.01, 500 at
      4.05); no feasible point found at k = 4 or 4.005 (residual stalls at 2-3e-5). The paper excludes dimension 4.
  Shrikhande pattern (16 rays, local theta 2, not realizable in C^4): feasible at k = 4, infeasible at 3.5.

usage: python3 rankhier.py shrikhande K...      (the Mermin / Shrikhande pattern of RESULTS.md §3.10, n = 16, d = 4)
       python3 rankhier.py shrikhande+bases K...  (with the completeness identities of the four bases)
       python3 rankhier.py fig3 K...            (the 11-vertex graph of the paper's Fig. 3)
       python3 rankhier.py c5 K...              (5-cycle: minimum dimension 3, Lovász bound 2.236)
"""
import sys
import time

import cvxpy as cp
import numpy as np
import scipy.sparse as sp


def feasible(A, k, eps=1e-7, max_iters=200000, verbose=False, bases=None, scale=None, margin=True, solver="SCS"):
    """bases: optional list of vertex lists that must be complete orthonormal bases of C^k (then k = their size).
    Each adds the degree-2 identity sum_{m in B} Gamma_im Gamma_mj = Gamma_ij, which is linear in M:
    sum_{m in B} M[(i,m),(m,j)] = (1/n) sum_l M[(i,l),(j,l)].
    scale: right-hand side of the normalisation (default N = n^2). margin: maximise the common PSD margin t
    (t < 0 means infeasible); margin=False solves the pure feasibility problem. solver: "SCS" (first order, eps,
    max_iters) or "CLARABEL" (interior point; accurate margins on small graphs, n <= 11 or so)."""
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

    r, s = cp.Variable(N * N), cp.Variable(N * N)
    t = cp.Variable() if margin else 0.0              # margin=False: pure feasibility (SCS certifies infeasibility)
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
    # normalisation k^2 Tr P + k Tr VP = scale; scale = N keeps the entries of P of order 1/k^2, so that a negative
    # margin is not lost in the solver tolerance (with scale = 1 they are of order 1/(k^2 N))
    cons.append(k * k * cp.sum(r[diag]) + k * cp.sum(r[pV][diag]) == (N if scale is None else scale))
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
    prob = cp.Problem(cp.Maximize(t) if margin else cp.Minimize(0), cons)
    if solver == "CLARABEL":
        prob.solve(solver="CLARABEL", verbose=verbose)
    else:
        prob.solve(solver="SCS", eps_abs=eps, eps_rel=eps, max_iters=max_iters, verbose=verbose)
    return prob.status, prob.value


def feasible_real(A, k, eps=1e-7, max_iters=200000, verbose=False, scale=None, margin=True):
    """Real orthonormal representations (the paper's Eq. (40) with the extra constraint rho^{T_A} = rho), reduced
    by O(k): rho = Phi_1 (x) 1 + Phi_V (x) V + Phi_W (x) W with W = sum |aa><bb| = V^{T_2}. With P = Phi_1:
    Phi_V = V P, P = P^T = V P V = P^{T_1}, Phi_W = Q = (VP)^{T_1} = V Q. PSD on the three isotypic blocks of
    {1, V, W}: P + VP + kQ (on the maximally entangled vector), P + VP, P - VP; PPT is implied by rho^{T_A} = rho.
    M = kP + VP + Q plays the role of Gamma (x) Gamma."""
    n = len(A)
    N = n * n
    pair = lambda i, j: i * n + j
    idx = lambda a, b: a * N + b
    sw = np.array([pair(j, i) for i in range(n) for j in range(n)])
    a_all, b_all = np.divmod(np.arange(N * N), N)
    pV = sw[a_all] * N + b_all
    ia, ja = np.divmod(a_all, n)
    ib, jb = np.divmod(b_all, n)
    pTA = (ib * n + ja) * N + (ia * n + jb)
    pVTA = (ja * n + ib) * N + (ia * n + jb)
    pT = b_all * N + a_all
    pVV = sw[a_all] * N + sw[b_all]
    r = cp.Variable(N * N)
    t = cp.Variable() if margin else 0.0
    mat = lambda v: cp.reshape(v, (N, N), order="C")
    I = np.eye(N)
    q = r[pVTA]
    cons = [r == r[pT], r == r[pVV], r == r[pTA], q == q[pV]]
    cons += [mat(r + r[pV] + k * q) - t * I >> 0, mat(r + r[pV]) - t * I >> 0, mat(r - r[pV]) - t * I >> 0]
    diag = np.arange(N) * (N + 1)
    cons.append(k * k * cp.sum(r[diag]) + k * cp.sum(r[pV][diag]) + k * cp.sum(q[diag]) == (N if scale is None else scale))
    Mr = k * r + r[pV] + q
    rows, cols, vals = [], [], []
    e = 0
    for i in range(n):
        for j in range(n):
            for jp in range(n):
                rows.append(e); cols.append(idx(pair(i, j), pair(i, jp))); vals.append(1.0)
                for m in range(n):
                    rows.append(e); cols.append(idx(pair(m, j), pair(m, jp))); vals.append(-1.0 / n)
                e += 1
    cons.append(sp.csr_matrix((vals, (rows, cols)), shape=(e, N * N)) @ Mr == 0)
    zero = [idx(pair(i, j), pair(ip, jp)) for i in range(n) for ip in range(n) if i != ip and not A[i, ip]
            for j in range(n) for jp in range(n)]
    if zero:
        cons.append(Mr[np.array(zero)] == 0)
    prob = cp.Problem(cp.Maximize(t) if margin else cp.Minimize(0), cons)
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


FIG3_EDGES = [(1, 2), (1, 3), (1, 5), (1, 7), (1, 8), (1, 10), (2, 3), (2, 4), (2, 6), (2, 11), (3, 4), (3, 5), (3, 9),
              (4, 5), (4, 6), (4, 7), (4, 9), (5, 6), (5, 8), (6, 7), (6, 8), (7, 9), (7, 11), (8, 10), (8, 11), (9, 10)]


def fig3():
    """The 11-vertex graph of Fig. 3 of the paper (vertices 1..11 as drawn there), theta = 4."""
    A = np.zeros((11, 11), dtype=bool)
    for i, j in FIG3_EDGES:
        A[i - 1, j - 1] = A[j - 1, i - 1] = True
    return A


if __name__ == "__main__":
    graph = sys.argv[1]
    complete = graph.endswith("+bases")
    graph = graph.replace("+bases", "")
    A = {"shrikhande": shrikhande_four_bases, "fig3": fig3}.get(graph, lambda: cycle(5))()
    bases = [list(range(4 * b, 4 * b + 4)) for b in range(4)] if complete else None
    for k in map(float, sys.argv[2:]):
        t0 = time.time()
        status, val = feasible(A, k, bases=bases)
        print(f"{graph}{' (complete bases)' if complete else ''}: k = {k:g}: status {status}, margin t = {val:.3e} "
              f"({time.time() - t0:.0f} s)", flush=True)

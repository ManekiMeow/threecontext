# Exact upper-bound certificate theta(Hbar) < 3.
# Lovasz dual: theta(G) <= lambda_max(A) for any symmetric A with A_ii = 1 and A_ij = 1
# whenever i != j are NOT adjacent in G.  For G = Hbar this means A_ij = 1 on edges of H
# (and on the diagonal), and A_ij free on non-edges of H.  We round a numerical optimum to
# rationals and verify t*I - A is positive definite exactly for a rational t < 3.
import numpy as np, cvxpy as cp
from fractions import Fraction as Fr
from gcert import exact_psd


def dual_upper(H):
    nodes = sorted(H); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
    A = cp.Variable((n, n), symmetric=True)
    t = cp.Variable()
    fixed = np.eye(n)
    for u, v in H.edges(): fixed[idx[u], idx[v]] = fixed[idx[v], idx[u]] = 1
    cons = [cp.multiply(fixed, A) == fixed, t * np.eye(n) - A >> 0]
    cp.Problem(cp.Minimize(t), cons).solve(solver="CLARABEL")
    return nodes, fixed, A.value, float(t.value)


def exact_upper(H, den=10**4):
    nodes, fixed, A, t = dual_upper(H)
    n = len(nodes)
    Aq = [[Fr(1) if fixed[i][j] else Fr(round(A[i][j] * den), den) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i): Aq[i][j] = Aq[j][i]
    # smallest rational t with a little slack that certifies t I - Aq PSD
    for slack in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2):
        tq = Fr(round((t + slack) * den), den)
        M = [[(tq if i == j else 0) - Aq[i][j] for j in range(n)] for i in range(n)]
        if exact_psd(M):
            return tq, t
    return None, t

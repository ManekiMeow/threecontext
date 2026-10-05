# Exact certificate that theta(complement of H) >= 3 for a 3-coloured triangle-free graph H.
#
# Certificate: rational weights c_i > 0 and a rational symmetric Y, zero diagonal,
# supported on E(H), such that
#   (i)  sum_{j in N(i), colour(j)=k} Y_ij = c_i  for every vertex i and colour k != colour(i),
#   (ii) B = Y + diag(c) is positive semidefinite.
# Then B / sum(c) is feasible for the Lovasz SDP of Hbar (B_ij = 0 on non-edges of H)
# with objective sum(B)/tr(B) = (sum c + 2 sum c)/sum c = 3, so theta(Hbar) >= 3, and
# theta(Hbar) <= chi(H) = 3.  (Equivalently W_ij = Y_ij/sqrt(c_i c_j) is a Hoffman-type
# matrix with W sqrt(c) = 2 sqrt(c), W + I PSD.)
#
# B always has the 2-dim kernel spanned by 1 - 3*1_{I_k}; we maximise the smallest
# eigenvalue of B on the complement of span{1_{I_k}} so that rational rounding inside
# the affine space defined by (i) keeps (ii) exactly.
import numpy as np, cvxpy as cp, sympy as sp
from fractions import Fraction as Fr


def margin_W(H, colour):
    nodes = sorted(H); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
    E = [(idx[u], idx[v]) for u, v in H.edges()]
    w = cp.Variable(len(E))
    c = cp.Variable(n)
    tau = cp.Variable()
    # constraint matrix for (i)
    rows = []
    for v in nodes:
        for k in range(3):
            if k == colour[v]: continue
            r = np.zeros(len(E))
            for e, (a, b) in enumerate(E):
                if (a == idx[v] and colour[nodes[b]] == k) or (b == idx[v] and colour[nodes[a]] == k):
                    r[e] = 1
            rows.append(r)
    A = np.array(rows)
    Wm = 0
    for e, (a, b) in enumerate(E):
        M = np.zeros((n, n)); M[a, b] = M[b, a] = 1
        Wm = Wm + w[e] * M
    C = np.zeros((n, 3))
    for v in nodes: C[idx[v], colour[v]] = 1
    Q, _ = np.linalg.qr(C)
    P = np.eye(n) - Q @ Q.T
    # rows of A correspond to (vertex, colour) pairs, in order; rhs is c_vertex
    R = np.zeros((len(rows), n)); r = 0
    for v in nodes:
        for k in range(3):
            if k == colour[v]: continue
            R[r, idx[v]] = 1; r += 1
    cons = [A @ w == R @ c, cp.sum(c) == 3, c >= 1e-3,
            P @ (Wm + cp.diag(c)) @ P - tau * P >> 0]
    prob = cp.Problem(cp.Maximize(tau), cons)
    prob.solve(solver="CLARABEL")
    return nodes, E, A, R, w.value, c.value, tau.value


def exact_psd(A):
    n = len(A); A = [row[:] for row in A]
    for i in range(n):
        if A[i][i] < 0: return False
        if A[i][i] == 0:
            if any(A[i][j] != 0 for j in range(i + 1, n)): return False
            continue
        for j in range(i + 1, n):
            if A[j][i] != 0:
                r = A[j][i] / A[i][i]
                for k in range(i, n): A[j][k] -= r * A[i][k]
    return True


def rational_certificate(H, colour, dens=(12, 60, 360, 2520, 27720, 10**6)):
    nodes, E, A, R, w, c, tau = margin_W(H, colour)
    if tau is None or tau <= 1e-7:
        return None, tau
    n = len(nodes)
    As = sp.Matrix(np.hstack([A, -R]).astype(int)).col_join(sp.Matrix([[0] * len(E) + [1] * n]))
    rhs = sp.Matrix([0] * (As.rows - 1) + [3])
    z = list(w) + list(c)
    for den in dens:
        zr = sp.Matrix([sp.Rational(round(x * den), den) for x in z])
        res = As * zr - rhs
        sol, params = As.gauss_jordan_solve(res)
        sol = sol.subs({p: 0 for p in params})
        zq = zr - sol
        assert As * zq == rhs
        wq, cq = zq[:len(E)], zq[len(E):]
        if min(cq) <= 0: continue
        W = [[Fr(0)] * n for _ in range(n)]
        for (a, b), x in zip(E, wq):
            fx = Fr(int(sp.fraction(x)[0]), int(sp.fraction(x)[1]))
            W[a][b] = W[b][a] = fx
        for i in range(n):
            x = cq[i]; W[i][i] = Fr(int(sp.fraction(x)[0]), int(sp.fraction(x)[1]))
        if exact_psd(W):
            return (nodes, E, [str(x) for x in wq], [str(x) for x in cq]), tau
    return None, tau

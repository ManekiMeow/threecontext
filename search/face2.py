# Optimal face of the Lovasz SDP for theta(Hbar):  is the optimal Gram matrix unique, and what
# ranks occur?  B0 = interior-point optimum (maximal rank in the face).  Every optimal B whose range
# lies in range(B0) is U M U^T with M >= 0 satisfying linear constraints (zeros on non-edges,
# trace 1, sum 3).  If that affine space is a point, the Gram matrix (hence the dimension) is
# determined by the graph.  Otherwise we walk to the boundary repeatedly: each step lowers the rank.
import sys, numpy as np, scipy.linalg as sl, networkx as nx
from lowrank import load
from theta import theta_comp


def constraints(G, n, U):
    k = U.shape[1]; iu = np.triu_indices(k)
    def lin(x, y):
        O = np.outer(x, y); O = O + O.T; O[np.diag_indices(k)] /= 2; return O[iu]
    rows = [lin(U[i], U[j]) for i in range(n) for j in range(i + 1, n) if not G.has_edge(i, j)]
    rows.append(np.eye(k)[iu])                     # trace (homogeneous part handled by B itself)
    s = U.T @ np.ones(n); rows.append(lin(s, s) / 2)
    return np.array(rows), iu


H = load(sys.argv[1]); nodes = sorted(H); n = len(nodes)
G0 = nx.relabel_nodes(H, {v: i for i, v in enumerate(nodes)}); G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(G0.edges())
t, B0 = theta_comp(G)
ev, V = np.linalg.eigh(B0); U = V[:, ev > 1e-6 * ev.max()]
A, iu = constraints(G, n, U)
null = A.shape[1] - np.linalg.matrix_rank(A, tol=1e-7 * np.linalg.norm(A, 2))
print(sys.argv[1], "theta=%.7f  max rank in optimal face=%d  face dimension=%d" % (t, U.shape[1], null), flush=True)
B = B0; rng = np.random.default_rng(0)
while True:
    ev, V = np.linalg.eigh(B); keep = ev > 1e-8 * ev.max(); U = V[:, keep]; k = U.shape[1]
    A, iu = constraints(G, n, U)
    N = sl.null_space(A, rcond=1e-9)
    if N.shape[1] == 0: break
    x = N @ rng.standard_normal(N.shape[1]); D = np.zeros((k, k)); D[iu] = x; D = D + D.T - np.diag(np.diag(D))
    M = U.T @ B @ U
    L = np.linalg.cholesky(M); Li = np.linalg.inv(L)
    mu = np.linalg.eigvalsh(Li @ D @ Li.T)
    step = -1 / mu.min() if mu.min() < 0 else -1 / mu.max()   # D traceless-like: one side hits boundary
    B = U @ (M + step * D) @ U.T; B = (B + B.T) / 2
mask_ok = max(abs(B[i, j]) for i in range(n) for j in range(i + 1, n) if not G.has_edge(i, j))
w = np.linalg.eigvalsh(B)
print("   extreme point reached: rank %d, min eig %.1e, tr %.6f, sum %.6f, max |B_ij| on non-edges %.1e"
      % (int((w > 1e-8 * w.max()).sum()), w.min(), np.trace(B), B.sum(), mask_ok), flush=True)

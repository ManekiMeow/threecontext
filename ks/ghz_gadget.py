# Four covering bases built from a three-context GHZ paradox (H36 or Perkel):
#   B_k = I_k  ∪  C_k     (C_k an orthonormal basis of S_k^perp),  k = 1,2,3
#   B_4 = {psi} ∪ C_4     (C_4 an orthonormal basis of psi^perp)
# Completions are chosen to maximise linked subspaces (internal bases):
#   * C_k contains an ONB of T_k = (S_k + S_{k+1})^perp, shared with C_{k+1} where possible,
#   * C_4 = ONB(S_1 ∩ psi^perp) ∪ C_1   (links S_1 between B_1 and B_4).
# Every orthonormal basis inside the union counts as a context; KS-colourability via SAT.
import sys, itertools, numpy as np, networkx as nx
from pysat.solvers import Cadical153
from ghzgeom import realisation_from_gram
sys.path.insert(0, "../search")

def onb(M, tol=1e-9):          # orthonormal basis of the row space of M
    if len(M) == 0: return np.zeros((0, M.shape[1]))
    u, s, vt = np.linalg.svd(M); return vt[: (s > tol * s[0]).sum()]
def perp(M, d):                # orthonormal basis of the orthogonal complement of rowspace(M)
    if len(M) == 0: return np.eye(d)
    u, s, vt = np.linalg.svd(M); r = (s > 1e-9 * s[0]).sum(); return vt[r:]
def complete(inside, space):   # extend ONB 'inside' (subset of 'space') to an ONB of 'space'
    if len(inside) == 0: return space
    rest = space - (space @ inside.T) @ inside
    return np.vstack([inside, onb(rest)])

def build(X, psi, cls):
    d = X.shape[1]
    I = [X[[i for i, c in enumerate(cls) if c == k]] for k in range(3)]
    I = [Ik / np.linalg.norm(Ik, axis=1, keepdims=True) for Ik in I]
    Sp = [perp(I[k], d) for k in range(3)]                       # S_k^perp
    T = [perp(np.vstack([I[k], I[(k + 1) % 3]]), d) for k in range(3)]   # (S_k+S_{k+1})^perp
    C = []
    for k in range(3):
        shared = [T[k]]
        prev = T[(k - 1) % 3]           # (S_{k-1}+S_k)^perp also lies in S_k^perp; include if orthogonal to T[k]
        if len(prev) and len(T[k]) and np.abs(prev @ T[k].T).max() < 1e-9: shared.append(prev)
        C.append(complete(onb(np.vstack(shared)) if sum(len(s) for s in shared) else np.zeros((0, d)), Sp[k]))
    A1 = perp(np.vstack([psi, Sp[0]]), d)                          # S_1 ∩ psi^perp
    B4 = np.vstack([psi, A1, C[0]])
    bases = [np.vstack([I[k], C[k]]) for k in range(3)] + [B4]
    for b in bases: assert np.allclose(b @ b.T, np.eye(d), atol=1e-8)
    return bases

def ks_analysis(bases, tol=1e-8, cap=10**6):
    d = bases[0].shape[1]
    vecs = []                                    # merge identical rays
    for b in bases:
        for v in b:
            if not any(abs(abs(v @ w) - 1) < 1e-8 for w in vecs): vecs.append(v)
    V = np.array(vecs); n = len(V)
    O = np.abs(V @ V.T) < tol
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from((i, j) for i in range(n) for j in range(i + 1, n) if O[i, j])
    ctx = [c for c in nx.find_cliques(G) if len(c) == d]
    s = Cadical153(); var = lambda i: i + 1
    for c in ctx:
        s.add_clause([var(i) for i in c])
        for a, b in itertools.combinations(c, 2): s.add_clause([-var(a), -var(b)])
    for i, j in G.edges(): s.add_clause([-var(i), -var(j)])
    sat = s.solve(); model = s.get_model() if sat else None
    return n, len(ctx), sat, model, V

if __name__ == "__main__":
    from verify_all import tricayley_z2z6
    S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
    Vv, E = tricayley_z2z6(S36); idx = {v: i for i, v in enumerate(Vv)}
    P = np.eye(36)
    for u, v in E: P[idx[u], idx[v]] = P[idx[v], idx[u]] = 0.25
    cls = [v[0] for v in Vv]
    X, psi = realisation_from_gram(P, cls)
    bases = build(X, psi, cls)
    n, nctx, sat, model, V = ks_analysis(bases)
    print("H36 construction: rays", n, " internal orthonormal bases", nctx, " KS-colourable:", sat)
    if sat:
        ones = [i for i in range(n) if model[i] > 0]
        print("  a surviving colouring has 1s on", len(ones), "rays; overlap with psi of each 1-ray:",
              [round(abs(V[i] @ psi), 4) for i in ones])

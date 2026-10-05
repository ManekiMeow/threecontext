# Geometry of a three-context GHZ realisation: S_k = span(I_k), their intersections and sums.
import sys, numpy as np, networkx as nx
sys.path.insert(0, "../search")

def realisation_from_gram(P, classes):
    """P = Gram matrix (unit diagonal, zeros on non-edges); returns vectors (rows) in R^r and psi."""
    w, U = np.linalg.eigh(P); keep = w > 1e-9
    X = U[:, keep] * np.sqrt(w[keep])            # rows = vectors, X X^T = P
    k0 = [i for i, c in enumerate(classes) if c == 0]
    # psi = normalised sum of h_i = sqrt(c_i) g_i; with uniform c within the class: psi ∝ sum g_i
    psi = X[k0].sum(0); psi /= np.linalg.norm(psi)
    return X, psi

def dim(M, tol=1e-8):
    if len(M) == 0: return 0
    s = np.linalg.svd(np.atleast_2d(M), compute_uv=False); return int((s > tol * s[0]).sum())

def report(name, X, psi, classes):
    d = X.shape[1]
    S = [X[[i for i, c in enumerate(classes) if c == k]] for k in range(3)]
    for k in range(3):
        a = np.abs(S[k] @ psi) ** 2
        print(f"  class {k}: |I|={len(S[k])}  sum |<psi|g>|^2 = {a.sum():.6f}")
    print(f"{name}: d = {d}, dims S_k = {[dim(s) for s in S]}")
    for k, l in [(0, 1), (0, 2), (1, 2)]:
        ds = dim(np.vstack([S[k], S[l]]))
        print(f"  dim(S{k}+S{l}) = {ds},  dim(S{k}∩S{l}) = {dim(S[k]) + dim(S[l]) - ds},  dim((S{k}+S{l})^perp) = {d - ds}")
    print(f"  dim(S0+S1+S2) = {dim(np.vstack(S))}")

if __name__ == "__main__":
    from verify_all import tricayley_z2z6
    S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
    V, E = tricayley_z2z6(S36); idx = {v: i for i, v in enumerate(V)}
    P = np.eye(36)
    for u, v in E: P[idx[u], idx[v]] = P[idx[v], idx[u]] = 0.25
    cls = [v[0] for v in V]
    X, psi = realisation_from_gram(P, cls); report("H36", X, psi, cls)
    sys.path.insert(0, "../search"); from perkel import perkel
    G = perkel(); Vp = sorted(G); ip = {v: i for i, v in enumerate(Vp)}
    P = np.eye(57)
    for u, v in G.edges(): P[ip[u], ip[v]] = P[ip[v], ip[u]] = 1 / 3
    cls = [v[0] for v in Vp]
    X, psi = realisation_from_gram(P, cls); report("Perkel", X, psi, cls)
    np.save("h36_vectors.npy", realisation_from_gram(*(lambda: (np.eye(1), [0]))()) [0]) if False else None

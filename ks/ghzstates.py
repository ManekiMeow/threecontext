# Largest family of mutually orthogonal GHZ states for three completed contexts.
#
# Given a three-context GHZ realisation (vectors I_1, I_2, I_3 in C^d, state psi0), complete each I_k to an
# orthonormal basis B_k of C^d.  A unit vector psi is a GHZ state for (B_1, B_2, B_3) if its support
# {r in B_1 ∪ B_2 ∪ B_3 : <r|psi> != 0} contains no triangle of pairwise non-orthogonal rays (one per basis).
# A four-context KS set with B_1, B_2, B_3 (Route A) needs d mutually orthogonal GHZ states (the fourth basis).
#
# Search: minimise f(psi) = sum over triangles (a,b,c) of |<a|psi>|^2 |<b|psi>|^2 |<c|psi>|^2 over unit psi
# orthogonal to the states found so far; snap the minimiser to the subspace orthogonal to its (numerically)
# zero rays and verify the triangle-free support exactly; repeat greedily with restarts.
import sys, numpy as np, torch, networkx as nx
sys.path.insert(0, "../search")
torch.set_default_dtype(torch.float64)
TOL = 1e-8


def realisation(P, cls):
    w, U = np.linalg.eigh(P); keep = w > 1e-9
    X = (U[:, keep] * np.sqrt(w[keep])).astype(complex)
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    k0 = [i for i, c in enumerate(cls) if c == 0]
    psi = X[k0].T @ (X[k0] @ np.ones(X.shape[1]))     # placeholder, replaced below
    # psi0 = top eigenvector of F = sum_k Pi_k (eigenvalue 3)
    F = X.conj().T @ X
    ev, V = np.linalg.eigh(F); psi = V[:, -1]
    assert abs(ev[-1] - 3) < 1e-6, ev[-1]
    return X, psi


def perp(M, d):
    if len(M) == 0: return np.eye(d, dtype=complex)
    u, s, vh = np.linalg.svd(M); r = (s > 1e-9 * s[0]).sum(); return vh[r:].conj()


def symmetry_unitary(X, perm):
    """Unitary U on C^d with U g_i = g_perm(i) (exists since the Gram matrix is invariant)."""
    Y = X[perm]
    U = (np.linalg.pinv(X) @ Y).T                    # solve X U^T = Y in least squares
    assert np.allclose(X @ U.T, Y, atol=1e-6)
    return U


def completions(X, cls, mode, rng, U=None):
    d = X.shape[1]; B = []
    for k in range(3):
        I = X[[i for i, c in enumerate(cls) if c == k]]
        Sp = perp(I, d)                              # rows: ONB of S_k^perp
        if mode == "random":
            Q, _ = np.linalg.qr(rng.standard_normal((len(Sp), len(Sp))) + 1j * rng.standard_normal((len(Sp), len(Sp))))
            C = Q.T @ Sp
        else:                                        # eigenvectors of the symmetry restricted to S_k^perp
            M = Sp.conj() @ U @ Sp.T                  # matrix of U on S_k^perp (rows of Sp as basis)
            ev, W = np.linalg.eig(M)
            W, _ = np.linalg.qr(W)                   # orthonormalise inside degenerate eigenspaces
            C = W.T @ Sp
        Bk = np.vstack([I, C])
        assert np.allclose(Bk.conj() @ Bk.T, np.eye(d), atol=1e-7), "not an ONB"
        B.append(Bk)
    return B


def triangle_tensors(B):
    O = lambda a, b: (np.abs(a.conj() @ b.T) > TOL).astype(float)
    return O(B[0], B[1]), O(B[1], B[2]), O(B[0], B[2])


def ghz_check(psi, B, E12, E23, E13):
    s = [np.abs(b.conj() @ psi) > 1e-7 for b in B]
    sub12 = E12[np.ix_(s[0], s[1])]; sub23 = E23[np.ix_(s[1], s[2])]; sub13 = E13[np.ix_(s[0], s[2])]
    tri = (sub12 @ sub23 * sub13).sum() if sub12.size and sub23.size and sub13.size else 0
    return tri == 0, [int(x.sum()) for x in s]


def find_family(B, psi0, restarts=40, steps=1500, seed=0):
    d = B[0].shape[1]; E12, E23, E13 = triangle_tensors(B)
    Bt = [torch.tensor(b) for b in B]; E = [torch.tensor(e) for e in (E12, E23, E13)]
    fam = [psi0]
    assert ghz_check(psi0, B, E12, E23, E13)[0]
    gen = torch.Generator().manual_seed(seed)
    while len(fam) < d:
        Pm = torch.tensor(perp(np.array(fam).conj(), d))    # rows: ONB of complement of the family
        best = None
        for r in range(restarts):
            z = torch.randn(2, Pm.shape[0], generator=gen, requires_grad=True)
            opt = torch.optim.Adam([z], lr=0.05)
            for it in range(steps):
                c = torch.complex(z[0], z[1]); psi = Pm.T @ c; psi = psi / psi.norm()
                p = [(b.conj() @ psi).abs() ** 2 for b in Bt]
                f = p[0] @ ((E[0] * (E[2] @ torch.diag(p[2]) @ E[1].T)) @ p[1])
                opt.zero_grad(); f.backward(); opt.step()
            psi = (Pm.T @ torch.complex(z[0], z[1])).detach().numpy(); psi /= np.linalg.norm(psi)
            # snap: project onto the complement of the (near-)zero rays and of the family
            zero = np.vstack([b[np.abs(b.conj() @ psi) < 1e-3] for b in B] + [np.array(fam).conj()])
            V = perp(zero.conj(), d)
            if len(V) == 0: continue
            cand = V.T @ (V.conj() @ psi); n = np.linalg.norm(cand)
            if n < 1e-6: continue
            cand /= n
            ok, sizes = ghz_check(cand, B, E12, E23, E13)
            if ok: best = (cand, sizes); break
        if best is None: break
        fam.append(best[0])
        print(f"    GHZ state #{len(fam)} found, support sizes {best[1]}", flush=True)
    return len(fam)


if __name__ == "__main__":
    from verify_all import tricayley_z2z6, tricayley
    from perkel import perkel
    graphs = {}
    S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
    V, E = tricayley_z2z6(S36); idx = {v: i for i, v in enumerate(V)}
    P = np.eye(36)
    for u, v in E: P[idx[u], idx[v]] = P[idx[v], idx[u]] = 0.25
    shift = [idx[(k, ((g[0]), (g[1] + 1) % 6))] for (k, g) in V]          # translation by (0,1) in Z2 x Z6
    graphs["H36"] = (P, [v[0] for v in V], shift)
    V, E = tricayley(13, ((0, 1, 3, 9), (0, 1, 10), (1, 6, 8))); idx = {v: i for i, v in enumerate(V)}
    P = np.eye(39); w = [0.25, 1 / 3, 1 / 3]
    for u, v in E: P[idx[u], idx[v]] = P[idx[v], idx[u]] = w[u[0]]
    shift = [idx[(k, (x + 1) % 13)] for (k, x) in V]
    graphs["H39"] = (P, [v[0] for v in V], shift)
    G = perkel(); V = sorted(G); idx = {v: i for i, v in enumerate(V)}
    P = np.eye(57)
    for u, v in G.edges(): P[idx[u], idx[v]] = P[idx[v], idx[u]] = 1 / 3
    shift = [idx[(k, (x + 1) % 19)] for (k, x) in V]
    graphs["Perkel"] = (P, [v[0] for v in V], shift)
    rng = np.random.default_rng(0)
    for name in sys.argv[1].split(","):
        P, cls, shift = graphs[name]
        X, psi0 = realisation(P, cls); d = X.shape[1]
        U = symmetry_unitary(X, shift)
        for mode in sys.argv[2].split(","):
            B = completions(X, cls, mode, rng, U)
            E12, E23, E13 = triangle_tensors(B)
            ntri = int((E12 @ E23 * E13).sum())
            print(f"{name} (d={d}) completion={mode}: triangles among the 3d rays = {ntri}", flush=True)
            m = find_family(B, psi0)
            print(f"{name} completion={mode}: largest orthogonal GHZ family found = {m} of d = {d}", flush=True)

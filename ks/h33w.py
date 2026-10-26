# W-symmetric four-basis configurations built from H33 (RESULTS.md 3.11).
# B2 = B-lines + {beta_k}, B3 = C-lines + {gamma_k} (block-pure completions); W = the unitary involution carrying the
# H33 copy of psi onto the dual H33 copy of the point a0 = (0,0); it fixes B2 and B3 as sets of rays.
# Choose any ONB A' of S1-perp, put B1 = points + A' and B4 = W(B1).  Then psi in B4 and a0 = W(psi) in B1 are GHZ
# states (triangle-free neighbourhoods isomorphic to H33), and W swaps B1 <-> B4.
import itertools, numpy as np
from scipy.optimize import linprog
import h33model as M, h33dual as Dl
from fourbasis import ortho_graph, all_bases

def build_W(a0=(0, 0)):
    a = M.point[a0]; Y = Dl.dual_points(a0)
    ph = lambda v, ref: v * abs(np.vdot(ref, v)) / np.vdot(ref, v)
    src, tgt = [M.psi], [a]
    for x in M.pts: src.append(M.point[x]); tgt.append(ph(Y[x], a))
    for i in range(4):
        kp, km = M.K[i], M.neg(M.K[i]); c0 = M.dot(M.K[i], a0)
        Bi = {0: M.Bline[(i, c0)], 1: M.beta[kp], 2: M.beta[km]}; Ci = {0: M.Cline[(i, c0)], 1: M.gamma[km], 2: M.gamma[kp]}
        for c in range(3):
            off = (c - c0) % 3
            src += [M.Bline[(i, c)], M.Cline[(i, c)]]; tgt += [ph(Bi[off], a), ph(Ci[off], a)]
    S = np.array(src).T; T = np.array(tgt).T; W = T @ np.linalg.pinv(S)
    assert np.allclose(W.conj().T @ W, np.eye(20), atol=1e-10) and np.allclose(W @ S, T, atol=1e-10)
    return W

W = build_W()
B2 = M.B + [M.beta[k] for k in M.nz]; B3 = M.C + [M.gamma[k] for k in M.nz]
H3 = [(M.u[0] + s1 * M.u[1] + s2 * M.u[2] + s1 * s2 * M.u[3]) / 2 for s1, s2 in [(1, -1), (-1, 1), (-1, -1)]]
A_SYM = H3 + [M.eps[k] for k in M.nz]

def configuration(Ap):
    B1 = M.A + list(Ap); B4 = [W @ v for v in B1]
    vecs = B1 + B2 + B3 + B4
    return vecs, [list(range(20 * j, 20 * j + 20)) for j in range(4)]

def k4_list(vecs, cover, tol=1e-9):
    V = np.array(vecs); ov = np.abs(V.conj() @ V.T) > tol
    return [t for t in itertools.product(*cover) if all(ov[x, y] for x, y in itertools.combinations(t, 2))]

def report(vecs, cover, tol=1e-9):
    V = np.array(vecs); Gm = np.abs(V.conj() @ V.T)
    for c in cover: assert np.allclose(Gm[np.ix_(c, c)], np.eye(20), atol=1e-8)
    shared = (int(np.sum(np.abs(Gm - 1) < 1e-8)) - 80) // 2
    G = ortho_graph(vecs, tol); bs = all_bases(G, 20); cov = {frozenset(c) for c in cover}
    internal = [b for b in bs if b not in cov]
    T = k4_list(vecs, cover, tol)
    good = [t for t in T if all(len(set(t) & b) == 1 for b in internal)]
    in_k4 = {r for t in T for r in t}
    return {"shared rays": shared, "internal": len(internal),
            "internal>=3 pieces": sum(1 for b in internal if sum(1 for c in cover if b & set(c)) >= 3),
            "K4": len(T), "rays in no K4": [len([r for r in c if r not in in_k4]) for c in cover],
            "KS colourings": len(good)}, T, internal

if __name__ == "__main__":
    vecs, cover = configuration(A_SYM)
    res, T, internal = report(vecs, cover); print("symmetric A':", res)
    rng = np.random.default_rng(5)
    for t in range(3):
        Q, _ = np.linalg.qr(np.array(M.A).T); P = np.eye(20) - Q @ Q.conj().T
        X = P @ (rng.normal(size=(20, 11)) + 1j * rng.normal(size=(20, 11))); Qc, _ = np.linalg.qr(X)
        res, T, internal = report(*configuration([Qc[:, j] for j in range(11)])); print("random A':", res)

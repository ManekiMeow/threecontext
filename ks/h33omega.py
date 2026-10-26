# The orbit Omega = Gamma.psi (10000 rays of C^20, RESULTS.md 3.11) with, for every ray s = g.psi, its own H33 copy:
# ctx9(s) = g(points) (9 rays of Omega), ctxB(s) = g(B-lines) (12 of the 20 rays of B2), ctxC(s) = g(C-lines) in B3.
# s is a GHZ state for (B1, B2, B3) iff ctx9(s) is contained in B1.
import numpy as np, pickle, os
import h33model as M
from h33group import gens
from h33orbit import rkey

def build():
    B2 = M.B + [M.beta[k] for k in M.nz]; B3 = M.C + [M.gamma[k] for k in M.nz]
    P = np.array([M.point[x] for x in M.pts]); Bl = np.array(M.B); Cl = np.array(M.C)
    seen = {rkey(M.psi): 0}; vecs = [M.psi]; mats = [np.eye(20, dtype=complex)]; frontier = [0]
    while frontier:
        nxt = []
        for j in frontier:
            for g in gens:
                U = g @ mats[j]; v = U @ M.psi; k = rkey(v)
                if k not in seen:
                    seen[k] = len(vecs); vecs.append(v); mats.append(U); nxt.append(len(vecs) - 1)
        frontier = nxt
    O = np.array(vecs)
    def find(v, pool):
        ov = np.abs(pool.conj() @ v); j = int(np.argmax(ov)); assert abs(ov[j] - 1) < 1e-8; return j
    B2a, B3a = np.array(B2), np.array(B3)
    ctx9 = np.array([[find(U @ p, O) for p in P] for U in mats])
    ctxB = np.array([[find(U @ b, B2a) for b in Bl] for U in mats])
    ctxC = np.array([[find(U @ c, B3a) for c in Cl] for U in mats])
    return O, ctx9, ctxB, ctxC, B2a, B3a

if __name__ == "__main__":
    O, ctx9, ctxB, ctxC, B2a, B3a = build()
    np.savez_compressed("h33_omega.npz", O=O, ctx9=ctx9, ctxB=ctxB, ctxC=ctxC, B2=B2a, B3=B3a)
    n = len(O); orth = np.abs(O.conj() @ O.T) < 1e-9
    print("Omega:", n, "rays; each orthogonal to", int(orth[0].sum()), "others")
    # sanity: ctx9(psi) = points, ctx9 rays are orthonormal and have overlap 1/3 with s
    s = 1234; c = ctx9[s]
    print("ctx9 orthonormal:", np.allclose(np.abs(O[c].conj() @ O[c].T), np.eye(9)), " |<s|ctx9>| =",
          sorted(set(np.round(np.abs(O[c].conj() @ O[s]), 6))))
    # how many distinct 9-contexts, and how often does a ray of Omega lie in some ray's 9-context
    keys = {tuple(sorted(r)) for r in ctx9}; print("distinct 9-contexts:", len(keys))
    psi_perp = np.where(orth[0])[0]; S1perp = np.where(orth[ctx9[0]].all(axis=0))[0]
    print("rays of Omega orthogonal to psi:", len(psi_perp), "; orthogonal to all 9 points:", len(S1perp))

# Block-aligned completions of H33 (RESULTS.md 3.11).  In the coordinates of h33model, C^20 = K0 + 8 blocks; each
# block C^2_k carries the six directions chi, gamma, e, eps, f, beta (0, 30, ..., 150 degrees).  A design chooses
#   * B4 inside each block: one orthogonal pair {chi,eps}, {gamma,f} or {e,beta};
#   * the completion rays of B1 (eps), B2 (beta), B3 (gamma) for the block pair {k, -k}: pure, or the real/imaginary
#     mixtures (x_k + x_-k)/sqrt2, (x_k - x_-k)/sqrt2  resp.  (x_k + i x_-k)/sqrt2, (x_k - i x_-k)/sqrt2;
#   * orthonormal bases of R = K0 - psi for B1 and B4.
# A KS colouring of the union is exactly a K4 (one pairwise non-orthogonal ray per covering basis) meeting every
# internal basis once, so the colourings are counted directly.
import itertools, numpy as np, networkx as nx
import h33model as M
from fourbasis import ortho_graph, all_bases

u = M.u; s2 = np.sqrt(2)
H3 = [(u[0] + s1 * u[1] + s2_ * u[2] + s1 * s2_ * u[3]) / 2 for s1, s2_ in [(1, -1), (-1, 1), (-1, -1)]]
def pairing(j):                       # {(u_a-u_b)/sqrt2, (u_c-u_d)/sqrt2, (u_a+u_b-u_c-u_d)/2}
    a, b, c, d = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)][j]
    return [(u[a] - u[b]) / s2, (u[c] - u[d]) / s2, (u[a] + u[b] - u[c] - u[d]) / 2]
RBASES = {"H": H3, "P0": pairing(0), "P1": pairing(1), "P2": pairing(2)}
PAIRS = {"chi": (M.chi, M.eps), "gam": (M.gamma, M.f), "e": (M.e, M.beta)}
def comp(dirs, kp, mode):
    k, m = kp; x, y = dirs[k], dirs[m]
    if mode == "pure": return [x, y]
    if mode == "re": return [(x + y) / s2, (x - y) / s2]
    return [(x + 1j * y) / s2, (x - 1j * y) / s2]
KPAIRS = [(M.K[i], M.neg(M.K[i])) for i in range(4)]

def build(design):
    """design: dict with 'B4' -> {k: pair name}, 'A','B','C' -> {class i: mode}, 'RA', 'RD' -> R-basis names."""
    Ap = list(RBASES[design["RA"]]); Bp = []; Cp = []; Dp = list(RBASES[design["RD"]])
    for i, kp in enumerate(KPAIRS):
        Ap += comp(M.eps, kp, design["A"][i]); Bp += comp(M.beta, kp, design["B"][i]); Cp += comp(M.gamma, kp, design["C"][i])
        for k in kp:
            a, b = PAIRS[design["B4"][k]]; Dp += [a[k], b[k]]
    bases = [M.A + Ap, M.B + Bp, M.C + Cp, [M.psi] + Dp]
    vecs = [v for b in bases for v in b]
    return vecs, [list(range(20 * j, 20 * j + 20)) for j in range(4)]

def analyse(vecs, cover, tol=1e-9):
    V = np.array(vecs); Gm = np.abs(V.conj() @ V.T)
    for c in cover: assert np.allclose(Gm[np.ix_(c, c)], np.eye(20), atol=1e-8)
    shared = int(np.sum(np.abs(Gm - 1) < 1e-8)) - 80                  # equal rays in different bases
    if shared: return {"shared": shared // 2}
    G = ortho_graph(vecs, tol); ov = Gm > tol
    bs = all_bases(G, 20); cov = {frozenset(c) for c in cover}
    internal = [b for b in bs if b not in cov]
    pcs = [sum(1 for c in cover if b & set(c)) for b in internal]
    T = [t for t in itertools.product(*cover) if all(ov[a, b] for a, b in itertools.combinations(t, 2))]
    good = [t for t in T if all(len(set(t) & b) == 1 for b in internal)]
    return {"internal": len(internal), "pieces>=3": sum(p >= 3 for p in pcs), "K4": len(T), "KS colourings": len(good),
            "B4 rays in a colouring": len({t[3] for t in good})}

if __name__ == "__main__":
    import sys
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    best = None; seen = 0
    for trial in range(N):
        d = {"B4": {k: rng.choice(list(PAIRS)) for kp in KPAIRS for k in kp},
             "A": [rng.choice(["pure", "re", "im"]) for _ in range(4)], "B": [rng.choice(["pure", "re", "im"]) for _ in range(4)],
             "C": [rng.choice(["pure", "re", "im"]) for _ in range(4)]}
        d["RA"], d["RD"] = rng.choice(list(RBASES), 2, replace=False)
        res = analyse(*build(d))
        if "shared" in res: continue
        seen += 1
        if best is None or res["KS colourings"] < best[0]["KS colourings"]: best = (res, d); print(trial, res, flush=True)
    print("designs without shared rays:", seen, "best:", best[0])

# Reproduces every number of RESULTS.md 3.11 (building a four-context KS set from H33).  Runtime ~3 minutes.
import itertools, numpy as np
from scipy.optimize import linprog
import h33model as M, h33ctx as X, h33design as D, h33dual as Dl, h33w as HW
from h33group import gens
from h33orbit import orbit
rng = np.random.default_rng(2026)
cover = [list(range(20 * j, 20 * j + 20)) for j in range(4)]

def haar_perp(S):
    Q, _ = np.linalg.qr(np.array(S).T); P = np.eye(20) - Q @ Q.conj().T
    Z = P @ (rng.normal(size=(20, 20 - len(S))) + 1j * rng.normal(size=(20, 20 - len(S))))
    Qc, _ = np.linalg.qr(Z); return [Qc[:, j] for j in range(20 - len(S))]

def k4_mixture(vecs, rho):
    """Is the Born vector of rho a convex combination of K4 indicators?  (If yes, the union is not KS.)"""
    T = HW.k4_list(vecs, cover); V = np.array(vecs)
    p = np.real(np.einsum("ri,ij,rj->r", V.conj(), rho, V)); A = np.zeros((80, len(T)))
    for j, t in enumerate(T): A[list(t), j] = 1
    return len(T), linprog(np.zeros(len(T)), A_eq=A, b_eq=p, bounds=(0, None), method="highs").status == 0

print("1. Route A is impossible with H33: overlaps between parallel distinct B- and C-lines")
ov = sorted({round(float(abs(np.vdot(M.Bline[l], M.Cline[m]))), 9) for l in M.Bline for m in M.Cline if l[0] == m[0] and l[1] != m[1]})
print("   |<b|c>| =", ov, "> 1/3, so by Lemma 14 every such pair lies in a K4 of any four-basis completion")

print("2. Structure: class subspaces V_i, block-diagonal completions, equiangular lines, duality W, orbit")
idx = {i: [M.U0[i]] + [M.CHI[k] for k in (M.K[i], M.neg(M.K[i]))] + [M.EPS[k] for k in (M.K[i], M.neg(M.K[i]))] for i in range(4)}
for name, Bs in [("B2 = B-lines + betas", HW.B2), ("B3 = C-lines + gammas", HW.B3)]:
    print(f"   {name}: rays per class subspace", [sum(1 for v in Bs if np.linalg.norm(v[idx[i]]) > 1e-9) for i in range(4)],
          "(each ray inside one V_i)" if all(sum(np.linalg.norm(v[idx[i]]) > 1e-9 for i in range(4)) == 1 for v in Bs) else "")
O, ok = orbit(M.psi, gens); O = np.array(O)
G6 = np.abs(O.conj() @ O.T) * 6
print(f"   orbit of psi under <AGL(2,3), W>: {len(O)} rays (closed: {ok}); 6|<u|v>| in",
      sorted({int(x) for x in np.round(G6.ravel())}), "; max deviation", f"{float(np.abs(G6 - np.round(G6)).max()):.1e}")
comp = O[:, idx[0]] * 2; reps = []
for v in comp:
    if not any(abs(abs(np.vdot(r, v)) - 1) < 1e-8 for r in reps): reps.append(v)
Gr = np.abs(np.array(reps).conj() @ np.array(reps).T)
print(f"   components in V_0: all of norm 1/2, {len(reps)} distinct lines, mutual |overlap| {sorted({round(float(x), 6) for x in Gr[np.triu_indices(len(reps), 1)]})}")
W = HW.W
print("   W unitary involution:", np.allclose(W.conj().T @ W, np.eye(20)) and np.allclose(W @ W, (W @ W)[0, 0] * np.eye(20)),
      "; W.psi = point (0,0):", abs(abs(np.vdot(W @ M.psi, M.point[(0, 0)])) - 1) < 1e-9)
Ys = {a: np.array([v for x, v in Dl.dual_points(a).items() if x != (0, 0)]) for a in M.pts}
nd = len({tuple(np.round(v, 6)) for a in M.pts for v in Ys[a]})
pairs = [p for p in itertools.combinations(M.pts, 2) if np.abs(Ys[p[0]].conj() @ Ys[p[1]].T).max() < 1e-9]
print(f"   dual 8-sets of the 9 points: {nd} distinct vectors; mutually orthogonal pairs of 8-sets: {len(pairs)}")

print("3. Per-ray test (NCF) and the state-averaging test")
th = X.Three(*X.completions("random", rng)); vecs = M.A + list(th.bases[0][9:]) + th.bases[1] + th.bases[2]
Dp = haar_perp([M.psi]); vals = [th.ncf(y)[0] for y in Dp]
print(f"   random completion + random B4: NCF of the 19 rays other than psi in [{min(vals):.6f}, {max(vals):.6f}] (1 = colourable with that ray = 1)")
for kind in ["symmetric", "random"]:
    th = X.Three(*X.completions(kind, rng)); v, lam, y = X.min_ratio(th, restrict=True, starts=10)
    print(f"   {kind} completion: a state orthogonal to psi with NCF = {th.ncf(y)[0]:.4f} exists (contextual)")
rho = (np.eye(20) - np.outer(M.psi, M.psi.conj())) / 19
for kind in ["symmetric", "random"]:
    Ap, Bp, Cp = X.completions(kind, rng)
    nT, mix = k4_mixture(M.A + list(Ap) + M.B + list(Bp) + M.C + list(Cp) + [M.psi] + haar_perp([M.psi]), rho)
    print(f"   {kind} completion + random B4: {nT} K4s; uniform state on psi-perp is a K4 mixture: {mix}")
d = {"B4": {k: "e" for kp in D.KPAIRS for k in kp}, "A": ["pure"] * 4, "B": ["re"] * 4, "C": ["pure"] * 4, "RA": "P0", "RD": "P1"}
vecs, _ = D.build(d); print("   block-aligned design:", D.analyse(vecs, cover), "; uniform state on psi-perp is a K4 mixture:", k4_mixture(vecs, rho)[1])

print("4. W-symmetric designs B4 = W(B1): psi and the point (0,0) are GHZ states")
res, T, internal = HW.report(*HW.configuration(HW.A_SYM)); print("   symmetric A':", res)
two = np.array([M.psi, M.point[(0, 0)]]); Q = np.linalg.qr(two.T, mode="complete")[0][:, 2:]
rho2 = Q @ Q.conj().T / 18
print("   uniform state on {psi, a0}-perp is a K4 mixture:", k4_mixture(HW.configuration(HW.A_SYM)[0], rho2)[1])

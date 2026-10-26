# Can a second copy of the H33 configuration (state y orthogonal to psi) live in the SAME three bases B1, B2, B3?
# Compatibility: in every basis, the two contexts placed there form an orthonormal set up to shared rays, i.e. every
# overlap |<u|Uv>| is 0 or 1.  Search: real orthogonal U = expm(X), penalty
#   sum_k sum_{u,v} m^2 (1 - m^2) + <psi|U psi>^2,  m = <u|U v>,
# minimised from random starts for each placement of the second copy's 9-context (B1, B2 or B3).
import sys, numpy as np
from scipy.linalg import expm, expm_frechet
from scipy.optimize import minimize
import h33model as Mo

V = Mo.A + Mo.B + Mo.C + [Mo.psi]
G = np.real(np.array([[np.vdot(x, y) for y in V] for x in V]))
ev, Uv = np.linalg.eigh(G); keep = ev > 1e-9
R = (Uv[:, keep] * np.sqrt(ev[keep])).T          # columns: real vectors in R^20 with the same Gram matrix
A, B, C, psi = R[:, :9], R[:, 9:21], R[:, 21:33], R[:, 33]
ctx = [A, B, C]

def make(place):
    # place = basis index (0,1,2) receiving the second copy's 9-context; its 12-contexts fill the other two
    others = [k for k in range(3) if k != place]
    Q = {place: A, others[0]: B, others[1]: C}
    def f(x):
        X = np.zeros((20, 20)); X[np.triu_indices(20, 1)] = x; X = X - X.T
        U = expm(X); val = 0.0; Gr = np.zeros((20, 20))
        for k in range(3):
            Mk = ctx[k].T @ U @ Q[k]
            val += np.sum(Mk ** 2 - Mk ** 4); Gr += ctx[k] @ (2 * Mk - 4 * Mk ** 3) @ Q[k].T
        s = psi @ U @ psi; val += s ** 2; Gr += 2 * s * np.outer(psi, psi)
        GX = expm_frechet(X.T, Gr, compute_expm=False); GX = GX - GX.T
        return val, GX[np.triu_indices(20, 1)]
    return f, Q

if __name__ == "__main__":
    restarts = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    rng = np.random.default_rng(0)
    for place in range(3):
        f, Q = make(place); best = None
        for r in range(restarts):
            x0 = rng.normal(scale=1.0, size=190)
            res = minimize(f, x0, jac=True, method="L-BFGS-B", options={"maxiter": 3000})
            if best is None or res.fun < best.fun: best = res
        print(f"9-context of the second copy in B{place + 1}: best penalty {best.fun:.3e} over {restarts} starts", flush=True)

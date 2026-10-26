# Contextual states for three completed H33 bases (RESULTS.md 3.11).  For a state y let p_r = |<r|y>|^2 on
# B1 u B2 u B3 and NCF(y) = max{ sum q_T : sum_{T contains r} q_T <= p_r }, T = pairwise non-orthogonal transversals.
# A ray y of B4 with NCF(y) = 1 gives a KS colouring with y = 1, so in a KS set every ray of B4 has NCF < 1.
# NCF(y) = min{ lambda.p(y) : lambda >= 0, lambda(T) >= 1 for all T }; min_ratio() minimises the smoothed ratio
# lambda_min(sum_r lambda_r |r><r| on psi-perp) / min_T lambda(T) to find contextual states orthogonal to psi.
import itertools, numpy as np
from scipy.optimize import linprog
import h33model as M

def completions(kind, rng=None):
    if kind == "symmetric":                                      # the Aut(H33)-invariant completion
        h = [(M.u[0] + s1 * M.u[1] + s2 * M.u[2] + s1 * s2 * M.u[3]) / 2 for s1, s2 in [(1, -1), (-1, 1), (-1, -1)]]
        return h + [M.eps[k] for k in M.nz], [M.beta[k] for k in M.nz], [M.gamma[k] for k in M.nz]
    out = []
    for S in (M.A, M.B, M.C):                                    # Haar-random ONB of each complement
        Q, _ = np.linalg.qr(np.array(S).T); P = np.eye(20) - Q @ Q.conj().T
        X = P @ (rng.normal(size=(20, 20 - len(S))) + 1j * rng.normal(size=(20, 20 - len(S))))
        Qc, _ = np.linalg.qr(X); out.append([Qc[:, j] for j in range(20 - len(S))])
    return out

class Three:
    def __init__(self, Ap, Bp, Cp, tol=1e-9):
        self.bases = [M.A + list(Ap), M.B + list(Bp), M.C + list(Cp)]
        self.V = np.array([v for b in self.bases for v in b])
        ov = np.abs(self.V.conj() @ self.V.T) > tol
        self.T = [t for t in itertools.product(range(20), range(20, 40), range(40, 60))
                  if ov[t[0], t[1]] and ov[t[0], t[2]] and ov[t[1], t[2]]]
        self.Aub = np.zeros((len(self.T), 60))
        for j, t in enumerate(self.T): self.Aub[j, list(t)] = 1
        Q = np.linalg.qr(np.array([M.psi]).T, mode="complete")[0]; self.Pperp = Q[:, 1:]   # basis of psi-perp
    def ncf(self, y):
        p = np.abs(self.V.conj() @ y) ** 2
        res = linprog(p, A_ub=-self.Aub, b_ub=-np.ones(len(self.T)), bounds=(0, None), method="highs")
        return res.fun, res.x

def min_ratio(th, restrict=True, starts=20, tau=0.02, seed=0):
    """min over lambda >= 0 of  lambda_min(O_lambda on psi-perp) / min_T lambda(T); a value < 1 certifies a contextual
    state orthogonal to psi (the bottom eigenvector).  Smooth surrogate: soft-min over triangles, lambda = theta^2."""
    from scipy.optimize import minimize
    P = th.Pperp if restrict else np.eye(20)
    Vc = th.V; A = th.Aub; rng = np.random.default_rng(seed); best = (np.inf, None, None)
    def f(t):
        lam = t ** 2
        O = (Vc.T * lam) @ Vc.conj(); Op = P.conj().T @ O @ P
        w, U = np.linalg.eigh((Op + Op.conj().T) / 2); y = P @ U[:, 0]
        g_eig = np.abs(Vc.conj() @ y) ** 2                       # d lambda_min / d lambda_r
        s = A @ lam; m = s.min(); e = np.exp(-(s - m) / tau); Z = e.sum()
        soft = m - tau * np.log(Z); g_soft = (e / Z) @ A          # d softmin / d lambda
        val = w[0] / soft; g = (g_eig * soft - w[0] * g_soft) / soft ** 2
        return val, g * 2 * t
    for _ in range(starts):
        res = minimize(f, rng.uniform(0.5, 1.5, 60), jac=True, method="L-BFGS-B", options={"maxiter": 2000})
        lam = res.x ** 2; lam = lam / (A @ lam).min()             # exact feasibility
        O = (Vc.T * lam) @ Vc.conj(); Op = P.conj().T @ O @ P
        w, U = np.linalg.eigh((Op + Op.conj().T) / 2)
        if w[0] < best[0]: best = (w[0], lam, P @ U[:, 0])
    return best

if __name__ == "__main__":
    import sys
    rng = np.random.default_rng(7)
    for kind in ["symmetric", "random"]:
        th = Three(*completions(kind, rng)); v, lam, y = min_ratio(th, restrict=True, starts=10)
        print(f"{kind} completion: {len(th.T)} triangles; contextual state orthogonal to psi with NCF = {th.ncf(y)[0]:.4f}")

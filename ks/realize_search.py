"""Numerical realizability tests for four-basis zero patterns (RESULTS.md §3.12): can a 4d-vertex non-orthogonality
graph X (four parts of size d) be realized by four orthonormal bases of C^d?

  rrr      relaxed reflect-reflect (Elser): x <- x + beta (P_A(2 P_B(x) - x) - P_B(x)), beta = 1 is Douglas–Rachford.
           A = {4Q : Q a rank-d projector}, B = {identity diagonal blocks, zeros where X has no edge}; with --flat h
           every allowed overlap is also forced to modulus 1/sqrt(h).
  general  any realization: basis 3 standard, bases 0-2 the columns of free matrices V_k; L-BFGS on orthonormality
           plus the forbidden overlaps.
  flat     flat realizations only: G = I + S/sqrt(h) with S Hermitian, |S_ij| = 1 on the edges of X and 0 elsewhere.
           G is the Gram matrix of four bases iff S^2 = 2 sqrt(h) S + 3h I; minimise ||S^2 - 2 sqrt(h) S - 3h I||^2
           over the edge phases with L-BFGS, from random starts.
  control  a planted realizable pattern of the same density (calibration): four bases of (C^2)^m, basis k acting by a
           Hadamard on the qubits in Q_k, with |Q_k delta Q_l| = m - 2 for all k != l. It is flat with h = 2^(m-2).

Patterns are read from ks/cayley4_candidate_p*.json ("sets") or from a .npy adjacency matrix. A solution with residual
< 1e-9 is saved next to the log.

usage: python3 realize_search.py rrr PATTERN --beta 0.5 --iters 20000 --seed 1 [--flat H] [--out DIR]
       python3 realize_search.py general PATTERN --restarts 20 --seed 1 [--out DIR]
       python3 realize_search.py flat PATTERN --restarts 50 --seed 1 [--out DIR]
       python3 realize_search.py control --qubits 6 --out DIR       (writes control_m6.npy)
"""
import argparse
import json
import os
import time

import numpy as np
from scipy.optimize import minimize


def load_pattern(path):
    """Adjacency matrix (bool, 4d x 4d) of the non-orthogonality graph."""
    if path.endswith(".npy"):
        return np.load(path).astype(bool)
    rec = json.load(open(path))
    p = rec["p"]
    X = np.zeros((4 * p, 4 * p), dtype=bool)
    for key, S in rec["sets"].items():
        k, l = int(key[0]), int(key[1])
        for s in range(p):
            for t in S:
                X[k * p + s, l * p + (s + t) % p] = True
    return X | X.T


def control(m):
    """Planted pattern: bases B_k = H^{Q_k} (computational basis) of (C^2)^m, |Q_k delta Q_l| = m - 2 (k != l).
    Four such sets exist only for m = 4, 6 (density 1/4, like the candidates)."""
    assert m in (4, 6)
    q = (m - 2) // 2
    Q = [set(), set(range(0, 2 * q)), set(range(0, q)) | set(range(2 * q, 3 * q)), set(range(q, 3 * q))]
    d = 2 ** m
    H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
    B = []
    for s in Q:
        M = np.ones((1, 1))
        for j in range(m):
            M = np.kron(M, H if j in s else np.eye(2))
        B.append(M)
    V = np.hstack(B)
    G = V.conj().T @ V
    X = (np.abs(G) > 1e-9) & ~np.kron(np.eye(4, dtype=bool), np.ones((d, d), dtype=bool))
    return X, G


def rrr(X, beta, iters, seed, flat=0.0, log=print, tol=1e-9):
    d = X.shape[0] // 4
    n = 4 * d
    part = np.arange(n) // d
    diag = part[:, None] == part[None, :]
    zero = ~X & ~diag
    nz = X & ~diag
    eye = np.eye(n)
    rng = np.random.default_rng(seed)

    def PA(M):
        w, U = np.linalg.eigh((M + M.conj().T) / 2)
        U = U[:, -d:]
        return 4 * U @ U.conj().T

    def PB(M):
        M = M.copy()
        M[zero] = 0
        M[diag] = eye[diag]
        if flat:
            M[nz] = M[nz] / np.maximum(np.abs(M[nz]), 1e-15) / np.sqrt(flat)
        return M

    V = np.hstack([np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))[0] for _ in range(4)])
    x = V.conj().T @ V
    best, t0 = np.inf, time.time()
    for it in range(iters):
        pb = PB(x)
        x = x + beta * (PA(2 * pb - x) - pb)
        if it % 100 == 99:
            pa = PA(PB(x))
            r = np.linalg.norm(PB(pa) - pa)
            best = min(best, r)
            if it % 1000 == 999:
                log(f"seed {seed} beta {beta} flat {flat} it {it + 1}: residual {r:.5f} best {best:.5f} "
                    f"({time.time() - t0:.0f} s)")
            if r < tol:
                return r, PB(pa)
    return best, None


def flat_search(X, restarts, seed, log=print, h=None, maxiter=4000):
    n = X.shape[0]
    d = n // 4
    iu, ju = np.where(np.triu(X, 1))
    if h is None:
        h = int(round(X[0, d:2 * d].sum())) if X[0, d:2 * d].sum() else int(round(X[0].sum() / 3))
    sh = np.sqrt(h)
    rng = np.random.default_rng(seed)

    def build(theta):
        S = np.zeros((n, n), dtype=complex)
        z = np.exp(1j * theta)
        S[iu, ju] = z
        S[ju, iu] = z.conj()
        return S

    def fg(theta):
        S = build(theta)
        R = S @ S - 2 * sh * S
        R[np.diag_indices(n)] -= 3 * h
        M = S @ R + R @ S - 2 * sh * R
        return float(np.real(np.vdot(R, R))), -4 * np.imag(M[ju, iu] * S[iu, ju])

    best = (np.inf, None)
    for r in range(restarts):
        t0 = time.time()
        res = minimize(fg, rng.uniform(0, 2 * np.pi, len(iu)), jac=True, method="L-BFGS-B",
                       options=dict(maxiter=maxiter, maxcor=30, ftol=1e-15, gtol=1e-10))
        # residual on the scale of G = I + S/sqrt(h): ||G^2 - 4G|| = ||R|| / h
        rel = np.sqrt(res.fun) / h
        log(f"seed {seed} restart {r}: ||G^2 - 4G|| = {rel:.6f} after {res.nit} iterations ({time.time() - t0:.0f} s)")
        if rel < best[0]:
            best = (rel, res.x)
        if rel < 1e-9:
            break
    return best[0], (np.eye(n) + build(best[1]) / sh) if best[1] is not None else None


def general_search(X, restarts, seed, log=print, maxiter=20000):
    """Any realization: basis 3 is the standard basis (gauge), bases 0-2 are the columns of V_0, V_1, V_2. Minimise
    sum_k ||V_k^† V_k - I||^2 + sum over pairs of the squared overlaps that X forbids, with L-BFGS."""
    d = X.shape[0] // 4
    blk = lambda k, l: X[k * d:(k + 1) * d, l * d:(l + 1) * d]
    Z3 = [~blk(3, k) for k in range(3)]                       # (V_k)_{st} = <e_s|v^k_t> must vanish
    Zkl = {(k, l): ~blk(k, l) for k in range(3) for l in range(k + 1, 3)}
    rng = np.random.default_rng(seed)
    eye = np.eye(d)

    def unpack(x):
        return (x[:3 * d * d] + 1j * x[3 * d * d:]).reshape(3, d, d)

    def fg(x):
        V = unpack(x)
        G = np.zeros_like(V)
        f = 0.0
        for k in range(3):
            R = V[k].conj().T @ V[k] - eye
            f += np.real(np.vdot(R, R))
            G[k] += 2 * V[k] @ R
            E = np.where(Z3[k], V[k], 0)
            f += np.real(np.vdot(E, E))
            G[k] += E
        for (k, l), Z in Zkl.items():
            E = np.where(Z, V[k].conj().T @ V[l], 0)
            f += np.real(np.vdot(E, E))
            G[k] += V[l] @ E.conj().T
            G[l] += V[k] @ E
        G = 2 * G.reshape(-1)
        return f, np.concatenate([G.real, G.imag])

    best = (np.inf, None)
    for r in range(restarts):
        t0 = time.time()
        V0 = np.stack([np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))[0] for _ in range(3)])
        res = minimize(fg, np.concatenate([V0.real.ravel(), V0.imag.ravel()]), jac=True, method="L-BFGS-B",
                       options=dict(maxiter=maxiter, maxcor=30, ftol=1e-16, gtol=1e-12))
        log(f"seed {seed} restart {r}: objective {res.fun:.3e} after {res.nit} iterations ({time.time() - t0:.0f} s)")
        if res.fun < best[0]:
            best = (res.fun, res.x)
        if res.fun < 1e-12:                                   # the d = 64 control stops at ~1e-15 when solved
            break
    V = unpack(best[1])
    W = np.hstack([V[0], V[1], V[2], eye])                    # columns: bases 0, 1, 2, 3
    Gram = W.conj().T @ W
    return np.sqrt(best[0]), Gram


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["rrr", "flat", "general", "control"])
    ap.add_argument("pattern", nargs="?")
    ap.add_argument("--beta", type=float, default=0.5)
    ap.add_argument("--iters", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--flat", type=float, default=0.0)
    ap.add_argument("--restarts", type=int, default=20)
    ap.add_argument("--qubits", type=int, default=6)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.mode == "control":
        X, G = control(a.qubits)
        np.save(os.path.join(a.out, f"control_m{a.qubits}.npy"), X)
        d = X.shape[0] // 4
        print(f"control m = {a.qubits}: d = {d}, degree per other basis {X[0, d:2 * d].sum()}, "
              f"planted residual {np.linalg.norm(G @ G - 4 * G):.2e}")
        return
    X = load_pattern(a.pattern)
    tag = os.path.basename(a.pattern).split(".")[0]
    logf = open(os.path.join(a.out, f"{a.mode}_{tag}_s{a.seed}_b{a.beta}_f{a.flat:g}.log"), "a")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n")
        logf.flush()
    if a.mode == "rrr":
        r, G = rrr(X, a.beta, a.iters, a.seed, flat=a.flat, log=log)
    elif a.mode == "general":
        r, G = general_search(X, a.restarts, a.seed, log=log)
    else:
        r, G = flat_search(X, a.restarts, a.seed, log=log)
    log(f"done: best residual {r:.3e}")
    if r < 1e-6 and G is not None:
        log(f"check: ||G^2 - 4G|| = {np.linalg.norm(G @ G - 4 * G):.2e}, forbidden overlaps "
            f"{np.abs(G[~X & ~np.kron(np.eye(4, dtype=bool), np.ones((len(G) // 4,) * 2, dtype=bool))]).max():.2e}")
        np.save(os.path.join(a.out, f"SOLUTION_{a.mode}_{tag}_s{a.seed}.npy"), G)
        log("*** solution saved ***")
    if False:
        np.save(os.path.join(a.out, f"SOLUTION_{a.mode}_{tag}_s{a.seed}.npy"), G)
        log("*** solution saved ***")


if __name__ == "__main__":
    main()

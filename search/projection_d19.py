"""Checks for RESULTS.md, Sections 2.1-2.2 (projection bound and the dimension-19 case).

  python3 search/projection_d19.py              H33 equality structure, and the two d = 19 shapes of Lemma 33
  python3 search/projection_d19.py --adversarial  also an adversarial numerical test of Lemma 29
                                                  (all rank triples, D <= 6; a few minutes)
Lemma 29: 8 Re tr(ABC) >= rank A + rank B + rank C - 2D for orthogonal projections on C^D.
"""
import sys


def h33_structure():
    import itertools
    import numpy as np
    pts = [(x, y) for x in range(3) for y in range(3)]
    dirs = [(1, 0), (0, 1), (1, 1), (1, 2)]
    lines = []
    for (dx, dy) in dirs:
        for off in range(3):
            # line: points p with (normal . p) = off, normal = (-dy, dx)
            nx, ny = (-dy) % 3, dx % 3
            lines.append(frozenset(p for p in pts if (nx * p[0] + ny * p[1]) % 3 == off))
    assert len(set(lines)) == 12
    A = [("A", p) for p in pts]; B = [("B", l) for l in lines]; C = [("C", l) for l in lines]
    V = A + B + C
    n = len(V)
    G = np.eye(n)
    for i, j in itertools.combinations(range(n), 2):
        (ki, xi), (kj, xj) = V[i], V[j]
        if {ki, kj} == {"A", "B"} or {ki, kj} == {"A", "C"}:
            p, l = (xi, xj) if ki == "A" else (xj, xi)
            if p in l:
                G[i, j] = G[j, i] = 1 / np.sqrt(12)
        if {ki, kj} == {"B", "C"} and xi != xj and len(xi & xj) == 0:
            G[i, j] = G[j, i] = 0.5
    w, U = np.linalg.eigh(G)
    print("min eig", w.min().round(10), "rank", int((w > 1e-9).sum()))
    vecs = (U[:, w > 1e-9] * np.sqrt(w[w > 1e-9])).T          # columns = event vectors in R^20
    d = vecs.shape[0]
    psi = vecs[:, :9].sum(1) / 3
    print("|psi| =", np.linalg.norm(psi).round(12))
    proj = []
    for idx in (range(0, 9), range(9, 21), range(21, 33)):
        M = vecs[:, list(idx)]
        proj.append(M @ M.T)
    P = np.outer(psi, psi)
    for k, Pi in enumerate(proj):
        print(f"context {k}: Pi psi = psi: {np.allclose(Pi @ psi, psi)}; probs {np.round(sorted(set(np.round((vecs[:, [i for i in range(n) if V[i][0]=='ABC'[k]]].T @ psi) ** 2, 6))), 6)}")
    Ak = [Pi - P for Pi in proj]
    tr = np.trace(Ak[0] @ Ak[1] @ Ak[2])
    ranks = [int(round(np.trace(a))) for a in Ak]
    print("tr(A1A2A3) =", round(tr, 10), "ranks", ranks, "D =", d - 1, "Lemma 29 slack:", round(8 * tr - (sum(ranks) - 2 * (d - 1)), 10))
    S = sum(Ak); E = S - 1.5 * (np.eye(d) - P)
    ev = np.round(np.linalg.eigvalsh(E), 8)
    print("E eigenvalues (on psi-perp):", sorted(set(ev)), "multiplicities", {v: int((ev == v).sum()) for v in sorted(set(ev))})
    # principal angles between S_B and S_C
    for (a, b, name) in ((0, 1, "AB"), (0, 2, "AC"), (1, 2, "BC")):
        Qa = np.linalg.svd(vecs[:, [i for i in range(n) if V[i][0] == "ABC"[a]]], full_matrices=False)[0]
        Qb = np.linalg.svd(vecs[:, [i for i in range(n) if V[i][0] == "ABC"[b]]], full_matrices=False)[0]
        sv = np.round(np.linalg.svd(Qa.T @ Qb, compute_uv=False), 8)
        print(name, "cosines of principal angles:", {v: int((sv == v).sum()) for v in sorted(set(sv), reverse=True)})


def d19_shapes():
    import numpy as np
    from scipy.linalg import block_diag

    def lines_block(angles):
        return [np.outer(v, v) for v in (np.array([np.cos(a), np.sin(a)]) for a in angles)]

    def analyse(name, A):
        A1, A2, A3 = A
        D = len(A1)
        tr = np.trace(A1 @ A2 @ A3)
        ranks = [int(round(np.trace(a).real)) for a in A]
        slack = 8 * tr.real - (sum(ranks) - 2 * D)
        out = [f"{name}: D = {D}, ranks {ranks}, tr(A1A2A3) = {tr.real:.6f}{'+%.1ei' % tr.imag if abs(tr.imag) > 1e-12 else ''}, Lemma 29 slack {slack:.6f}"]
        for (k, l, m) in ((1, 2, 0), (0, 2, 1), (0, 1, 2)):
            X = np.eye(D) + 4 * (A[k] @ A[l] + A[l] @ A[k])           # 1 + 8 Re(A_k A_l)
            ev = np.sort(np.linalg.eigvalsh(X))
            out.append(f"   pair ({k+1},{l+1}): smallest eigenvalues of 1+8T {np.round(ev[:10], 4).tolist()}; tr(A_{m+1}(1+8T)) = {np.trace(A[m] @ X).real:.6f}")
        print("\n".join(out))

    # case (ii): nine blocks, lines at 0 (A2), arccos(1/3) (A3), and A1 along beta - gamma
    th = np.arccos(1 / 3)
    b, g = np.array([1.0, 0.0]), np.array([np.cos(th), np.sin(th)])
    al = (b - g) / np.linalg.norm(b - g)
    blk = [np.outer(al, al), np.outer(b, b), np.outer(g, g)]
    A = [block_diag(*[blk[k]] * 9) for k in range(3)]
    analyse("case (ii), uniform cosines (1/sqrt3, 1/3, 1/sqrt3)", A)
    E = A[0] + A[1] + A[2] - 1.5 * np.eye(18)
    print("   E eigenvalues:", sorted(set(np.round(np.linalg.eigvalsh(E), 6))))

    # case (i): eight trines + tail (e2, e3 orthogonal; v at angle phi, phase chi)
    tri = lines_block([2 * np.pi / 3, 0.0, np.pi / 3])            # A1, A2, A3 lines at 120 deg spacing (trine)
    for phi, chi in ((0.0, 0.0), (0.3, 0.0), (0.7, 1.1)):
        v = np.array([np.cos(phi), np.exp(1j * chi) * np.sin(phi)])
        tail = [np.outer(v, v.conj()), np.diag([1.0, 0.0]).astype(complex), np.diag([0.0, 1.0]).astype(complex)]
        A = [block_diag(*([tri[k].astype(complex)] * 8 + [tail[k]])) for k in range(3)]
        analyse(f"case (i), phi = {phi}, chi = {chi}", A)


def adversarial():
    import numpy as np
    from scipy.optimize import minimize
    rng = np.random.default_rng(1)

    def proj(x, D, r):
        M = (x[:D * r] + 1j * x[D * r:]).reshape(D, r)
        Qm, _ = np.linalg.qr(M)
        return Qm @ Qm.conj().T

    def f(x, D, ranks):
        k = 2 * D
        A = proj(x[:k * ranks[0]], D, ranks[0])
        B = proj(x[k * ranks[0]:k * (ranks[0] + ranks[1])], D, ranks[1])
        C = proj(x[k * (ranks[0] + ranks[1]):], D, ranks[2])
        return 8 * np.trace(A @ B @ C).real - (sum(ranks) - 2 * D)

    worst = {}
    for D in range(2, 7):
        for ra in range(1, D):
            for rb in range(ra, D):
                for rc in range(rb, D):
                    best = np.inf
                    for s in range(12):
                        x0 = rng.normal(size=2 * D * (ra + rb + rc))
                        res = minimize(f, x0, args=(D, (ra, rb, rc)), method="L-BFGS-B", options={"maxiter": 3000})
                        best = min(best, res.fun)
                    worst[(D, ra, rb, rc)] = best
    vals = sorted(worst.items(), key=lambda kv: kv[1])
    print("smallest slack found:", [(k, round(v, 6)) for k, v in vals[:6]])
    print("any negative?", any(v < -1e-7 for v in worst.values()), "cases", len(worst))


if __name__ == "__main__":
    print("== H33: equality in Lemma 29")
    h33_structure()
    print("== dimension 19: the two shapes of Lemma 33 (slack exactly 1)")
    d19_shapes()
    if "--adversarial" in sys.argv:
        print("== adversarial test of Lemma 29")
        adversarial()

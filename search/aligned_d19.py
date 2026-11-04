"""RESULTS.md Section 2.5, Proposition 43: when E^2 commutes with the three projections, the triple is aligned.

Write J_k = 2A_k - 1 and X = J_A + J_B + J_C = 2E on W = psi^perp (dimension 18). If X^2 commutes with every J_k, then
each anticommutator {J_k, J_l} is central, so every irreducible block has dimension <= 2 (Clifford relations).
In dimension 19 (shape (ii)) the triple is then a sum of nine 2-dimensional blocks carrying three lines each.
If E^2 is scalar, then E^2 = 1/36, every probability is 1/10, and the event Gram matrix G = I + Y satisfies
    Y^3 = (7/9) Y + (2/9)(J - I),        Re(Y_AB Y_BC Y_CA) = (1/9)(J - I) on the A-block.

Checks:
 1. random rank-9 triples on C^18 with S = A+B+C two-valued {4/3, 5/3} (found by least squares from random starts):
    the commutant has dimension 9 (nine pairwise non-isomorphic 2-dimensional blocks), the anticommutators are central;
 2. the identities for aligned data with all block defects 1/9 and event bases with uniform psi-components
    (any such data, triangle-free or not, satisfies them).
usage: python3 search/aligned_d19.py
"""
import numpy as np
from scipy.optimize import least_squares

rng = np.random.default_rng(0)
D, r = 18, 9


def projs(x):
    out = []
    for k in range(3):
        a = x[k * 2 * D * r:(k + 1) * 2 * D * r]
        V = np.linalg.qr((a[:D * r] + 1j * a[D * r:]).reshape(D, r))[0]
        out.append(V @ V.conj().T)
    return out


def two_valued_residual(x):
    S = sum(projs(x))
    M = (S - 4 / 3 * np.eye(D)) @ (S - 5 / 3 * np.eye(D))
    return np.concatenate([M.real.ravel(), M.imag.ravel()])


def commutant_dim(Ps):
    I = np.eye(D)
    M = np.vstack([np.kron(P, I) - np.kron(I, P.T) for P in Ps])
    s = np.linalg.svd(M, compute_uv=False)
    return int((s < 1e-7 * s.max()).sum())


print("1. two-valued S = A+B+C on C^18 (ranks 9,9,9), random solutions:")
for trial in range(3):
    s = least_squares(two_valued_residual, rng.normal(size=3 * 2 * D * r), method="trf", max_nfev=400,
                      xtol=1e-15, ftol=1e-15, gtol=1e-15)
    A, B, C = projs(s.x)
    J = [2 * P - np.eye(D) for P in (A, B, C)]
    central = all(np.allclose(K @ Jm, Jm @ K, atol=1e-8)
                  for K in (J[0] @ J[1] + J[1] @ J[0], J[1] @ J[2] + J[2] @ J[1], J[0] @ J[2] + J[2] @ J[0]) for Jm in J)
    tr = np.trace(A @ B @ C)
    print(f"   residual {np.abs(s.fun).max():.1e}; commutant dimension {commutant_dim([A, B, C])}; "
          f"anticommutators central: {central}; tr(ABC) = {tr.real:.6f}{tr.imag:+.4f}i")

print("2. aligned data with every block defect 1/9 (|n_A + n_B + n_C| = 1/3 on the Bloch sphere), uniform probabilities:")
sig = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1.0, -1.0])]


def line(n):
    P = (np.eye(2) + sum(c * s for c, s in zip(n, sig))) / 2
    w, V = np.linalg.eigh(P)
    return V[:, 1]


def block():
    while True:   # three unit vectors whose sum has length 1/3
        nA, nB = (v / np.linalg.norm(v) for v in rng.normal(size=(2, 3)))
        m = nA + nB
        # choose nC on the unit sphere with |m + nC| = 1/3:  m.nC = (1/9 - |m|^2 - 1)/2
        t = (1 / 9 - m @ m - 1) / 2
        if abs(t) < np.linalg.norm(m):
            u = m / np.linalg.norm(m)
            perp = np.cross(u, rng.normal(size=3)); perp /= np.linalg.norm(perp)
            a = t / np.linalg.norm(m)
            nC = a * u + np.sqrt(1 - a * a) * perp
            return [line(nA), line(nB), line(nC)]


blocks = [block() for _ in range(9)]
d = 19
psi = np.zeros(d, complex); psi[0] = 1


def context(k):
    Bk = np.zeros((d, 10), complex); Bk[0, 0] = 1
    for j in range(9):
        Bk[1 + 2 * j:3 + 2 * j, 1 + j] = blocks[j][k]
    # event basis with uniform psi-components: first row of U = 1/sqrt(10)
    M = np.linalg.qr(np.column_stack([np.ones(10) / np.sqrt(10), rng.normal(size=(10, 9)) + 1j * rng.normal(size=(10, 9))]))[0]
    M[:, 0] *= np.sign(M[0, 0].real)
    U = M.T.conj()          # rows: block coordinates (psi first), columns: events; first row uniform
    U[0] = np.abs(U[0])
    return Bk @ U


G_ = np.hstack([context(k) for k in range(3)])
G = G_.conj().T @ G_
Y = G - np.eye(30)
S_ops = [Gk @ Gk.conj().T - np.outer(psi, psi) for Gk in (G_[:, :10], G_[:, 10:20], G_[:, 20:])]
E = sum(S_ops) - 1.5 * (np.eye(d) - np.outer(psi, psi))
Jm = np.ones((30, 30))
print("   E^2 = 1/36 on psi-perp:", np.allclose(E @ E, (np.eye(d) - np.outer(psi, psi)) / 36),
      "| probabilities all 1/10:", np.allclose(np.abs(G_[0]) ** 2, 0.1))
print("   Y^3 = (7/9)Y + (2/9)(J - I):", np.allclose(Y @ Y @ Y, 7 / 9 * Y + 2 / 9 * (Jm - np.eye(30))))
YAB, YBC, YCA = Y[:10, 10:20], Y[10:20, 20:], Y[20:, :10]
M = YAB @ YBC @ YCA
print("   (M + M^dagger)/2 = (1/9)(J - I) for M = Y_AB Y_BC Y_CA:",
      np.allclose((M + M.conj().T) / 2, (np.ones((10, 10)) - np.eye(10)) / 9))
YAC = Y[:10, 20:]
Pm, Qm = YAB @ YAB.conj().T, YAC @ YAC.conj().T
print("3. aligned design identities: [P,Q] = 0", np.allclose(Pm @ Qm, Qm @ Pm), "| [P,M] = 0", np.allclose(Pm @ M, M @ Pm),
      "| P Q (I + M + M^dag - P - Q) = M M^dag", np.allclose(Pm @ Qm @ (np.eye(10) + M + M.conj().T - Pm - Qm), M @ M.conj().T))
print("   (diag(M) = 0, i.e. the per-event identity, holds only for triangle-free supports; these bases are random)")

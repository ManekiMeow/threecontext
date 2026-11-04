"""RESULTS.md Section 2.5, Lemma 42: identities of the transfer matrix Omega = Y_AB^{-1} Y_AC Y_CB in dimension 19.

Checked on random data in the generic shape (three random 10-dimensional context spans through psi in C^19, random
event bases with <g|psi> > 0); these are identities, so they need no triangle-freeness:
  * Z = Y Omega with Z = Y_AC Y_CB, Omega p_B = p_B, Y p_B = p_A;
  * Omega is the B-event matrix of T_B = (Pi_A restricted to S_B)^{-1} Pi_A Pi_C, and on ran A_B
    T_B = Q_B E - 1/2, with Q_B the oblique projection onto ran A_B along (ran A_A)^perp and E = A_A + A_B + A_C - 3/2;
  * aligned case (simultaneous 2-dimensional blocks): Omega is normal with spectrum {1, conj(tau_j)/|<a_j|b_j>|^2}.
usage: python3 search/transfer_d19.py
"""
import numpy as np

rng = np.random.default_rng(3)
d = 19
psi = np.zeros(d, complex)
psi[0] = 1


def events(Bk):
    """random orthonormal event basis of span(Bk) (Bk: d x 10 isometry with first column psi), <g|psi> > 0"""
    U = np.linalg.qr(rng.normal(size=(10, 10)) + 1j * rng.normal(size=(10, 10)))[0]
    U = U / (U[:, 0] / np.abs(U[:, 0]))[:, None]
    return Bk @ U.T


def generic_span():
    Q = np.linalg.qr(rng.normal(size=(d - 1, 9)) + 1j * rng.normal(size=(d - 1, 9)))[0]
    B = np.zeros((d, 10), complex)
    B[0, 0] = 1
    B[1:, 1:] = Q
    return B


spans = [generic_span() for _ in range(3)]
GA, GB, GC = (events(B) for B in spans)
Y, Z = GA.conj().T @ GB, GA.conj().T @ GC @ GC.conj().T @ GB
Om = np.linalg.solve(Y, Z)
pA, pB = (GA.conj().T @ psi).real, (GB.conj().T @ psi).real
print("generic: Z = Y Omega", np.allclose(Y @ Om, Z), "| Omega p_B = p_B", np.allclose(Om @ pB, pB),
      "| Y p_B = p_A", np.allclose(Y @ pB, pA))
# T_B on S_B and the formula T_B' = Q_B E - 1/2 on ran A_B
P = np.outer(psi, psi.conj())
Ak = [B @ B.conj().T - P for B in spans]
E = Ak[0] + Ak[1] + Ak[2] - 1.5 * (np.eye(d) - P)
VA, VB = spans[0][:, 1:], spans[1][:, 1:]                       # orthonormal bases of ran A_A, ran A_B
TB = GB @ Om @ GB.conj().T                                     # operator on S_B
M = VA.conj().T @ VB                                           # A_A restricted to ran A_B, in these bases (invertible)
QB = VB @ np.linalg.solve(M, VA.conj().T)                       # oblique projection onto ran A_B along (ran A_A)^perp
lhs = VB.conj().T @ TB @ VB
rhs = VB.conj().T @ (QB @ E - 0.5 * np.eye(d)) @ VB
print("generic: T_B = (Pi_A|S_B)^-1 Pi_A Pi_C on S_B",
      np.allclose(GA.conj().T @ (GA @ GA.conj().T) @ TB @ GB, GA.conj().T @ (GA @ GA.conj().T) @ (GC @ GC.conj().T) @ GB),
      "| T_B' = Q_B E - 1/2 on ran A_B", np.allclose(lhs, rhs))
# aligned case
blocks = [[np.array([np.cos(a), np.exp(1j * b) * np.sin(a)]) for a, b in rng.uniform(0, np.pi, size=(3, 2))]
          for _ in range(9)]


def aligned_span(k):
    B = np.zeros((d, 10), complex)
    B[0, 0] = 1
    for j in range(9):
        B[1 + 2 * j:3 + 2 * j, 1 + j] = blocks[j][k]
    return B


GA, GB, GC = (events(aligned_span(k)) for k in range(3))
Y, Z = GA.conj().T @ GB, GA.conj().T @ GC @ GC.conj().T @ GB
Om = np.linalg.solve(Y, Z)
tau = np.array([np.vdot(b[0], b[1]) * np.vdot(b[1], b[2]) * np.vdot(b[2], b[0]) for b in blocks])
lam = np.concatenate([[1], np.conj(tau) / np.abs([np.vdot(b[0], b[1]) for b in blocks]) ** 2])
print("aligned: Omega normal", np.allclose(Om @ Om.conj().T, Om.conj().T @ Om),
      "| spectrum {1, conj(tau_j)/|<a_j|b_j>|^2}", np.allclose(np.sort_complex(np.linalg.eigvals(Om)), np.sort_complex(lam)))

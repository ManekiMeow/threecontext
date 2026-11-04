"""RESULTS.md Section 2.5: group-covariant 30-event configurations in dimension 19 and why they fail.

A finite group acting unitarily, fixing psi up to a phase and transitive on each context forces equal probabilities
1/10. If the three context spans carry the same multiplicity-free permutation module 1 + sum(rho), then in d = 19 each
rho occurs with multiplicity 2, the projections are block diagonal (aligned) with blocks rho (x) C^2, and
Y_kl = E_1 + sum_rho a_kl(rho) E_rho with E_rho the primitive idempotents. Each block needs three lines in C^2, i.e. a
singular 3x3 Gram matrix of the a_kl(rho).

1. S5 on the 2-subsets of {1..5} (module 1 + V4 + V5): the only support with degrees 3-4 is disjointness for every
   pair (the Petersen relation). It forces a(V4) = -2/3 and a(V5) = 1/3, with Gram determinants -25/27 and 20/27.
2. Z10 (regular, module = all ten characters): the configuration is three functions y_AB, y_BC, y_AC on Z10, with sum 1
   and triangle-free supports. Each character chi != 0 needs
       Phi_chi = 1 + 2 Re(x y conj(z)) - |x|^2 - |y|^2 - |z|^2 = 0,
   where x, y, z are the Fourier coefficients. Triangle-freeness gives sum_chi Phi_chi = 10 (1 - sum ||y_kl||^2).
   The script reports the smallest residual over random starts for every support class whose degrees lie in {3, 4}
   (degrees 2 force |x(5)| in {0, 1}). This part is numerical illustration only. Section 2.5's theta certificates
   (search/symmetric_theta.py) exclude every Z10-invariant support in every dimension rigorously.

usage: python3 search/covariant_d19.py
"""
import itertools

import numpy as np
import sympy as sp
from scipy.optimize import least_squares


def petersen():
    P = list(itertools.combinations(range(5), 2))
    n = len(P)
    rel = lambda a, b: "E" if a == b else ("M" if len(set(a) & set(b)) == 1 else "D")
    D = sp.Matrix(n, n, lambda i, j: 1 if rel(P[i], P[j]) == "D" else 0)
    I = sp.eye(n)
    E1 = sp.ones(n, n) / 10
    E5 = (D - 3 * I) * (D + 2 * I) / ((1 - 3) * (1 + 2))          # Petersen eigenvalue 1
    E4 = (D - 3 * I) * (D - I) / ((-2 - 3) * (-2 - 1))             # Petersen eigenvalue -2
    assert E1 + E4 + E5 == I and E4 ** 2 == E4 and E5 ** 2 == E5 and E4.rank() == 4 and E5.rank() == 5
    a4, a5 = sp.symbols("a4 a5")
    Y = E1 + a4 * E4 + a5 * E5
    pos = {r: next((i, j) for i in range(n) for j in range(n) if rel(P[i], P[j]) == r) for r in "EMD"}
    sol = sp.solve([Y[pos["E"]], Y[pos["M"]]], [a4, a5], dict=True)[0]
    print("S5 / Petersen support: forced block overlaps a(V4) =", sol[a4], ", a(V5) =", sol[a5],
          "; overlap on disjoint pairs", sp.nsimplify(Y[pos["D"]].subs(sol)))
    for name, a in (("V4", sol[a4]), ("V5", sol[a5])):
        G = sp.Matrix([[1, a, a], [a, 1, a], [a, a, 1]])
        print(f"   block {name}: Gram determinant of the three lines = {G.det()} (needs 0)")


def z10():
    N = 10
    w = np.exp(2j * np.pi / N)
    F = np.array([[w ** (-c * s) for s in range(N)] for c in range(N)])

    def classes(sizes):
        seen, out = set(), []
        for SAB in itertools.combinations(range(N), sizes[0]):
            for SAC in itertools.combinations(range(N), sizes[1]):
                if 0 not in SAB or 0 not in SAC:
                    continue
                for SBC in itertools.combinations(range(N), sizes[2]):
                    if {(a + b) % N for a in SAB for b in SBC} & set(SAC):
                        continue
                    key = min((tuple(sorted((u * a + t1) % N for a in SAB)), tuple(sorted((u * a + t2) % N for a in SAC)),
                               tuple(sorted((u * a - t1 + t2) % N for a in SBC)))
                              for u in (1, 3, 7, 9) for t1 in range(N) for t2 in range(N)
                              if (-t1 * pow(u, -1, N)) % N in SAB and (-t2 * pow(u, -1, N)) % N in SAC)
                    if key not in seen:
                        seen.add(key)
                        out.append((SAB, SBC, SAC))
        return out

    def resid(v, sup):
        ys, p = [], 0
        for S in sup:
            k = len(S)
            z = v[p:p + k - 1] + 1j * v[p + k - 1:p + 2 * k - 2]
            p += 2 * k - 2
            y = np.zeros(N, complex)
            y[list(S[:-1])] = z
            y[S[-1]] = 1 - z.sum()
            ys.append(y)
        x, y, z = (F @ t for t in ys)
        x, y, z = x[1:], y[1:], z[1:]
        return 1 + 2 * np.real(x * y * np.conj(z)) - abs(x) ** 2 - abs(y) ** 2 - abs(z) ** 2

    rng = np.random.default_rng(0)
    worst = np.inf
    total = 0
    for sizes in itertools.product((3, 4), repeat=3):
        if sizes == (4, 4, 4):
            continue                                       # all events of type (4,4): Proposition 37
        for sup in classes(sizes):
            total += 1
            nv = sum(2 * len(S) - 2 for S in sup)
            best = min(least_squares(lambda v: resid(v, sup), rng.normal(size=nv) * 0.5, method="trf",
                                     max_nfev=2000).cost for _ in range(10))
            worst = min(worst, best)
    print(f"Z10-covariant, d = 19: {total} support classes with degrees in {{3,4}}; smallest residual "
          f"sum Phi^2/2 found = {worst:.4f} (a solution needs 0)")


if __name__ == "__main__":
    petersen()
    z10()

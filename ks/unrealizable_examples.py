"""KS-uncolourable four-basis skeletons with no quantum realization (RESULTS.md §3.14).

Each example is a 4-partite graph H on four classes of size d (the would-be orthonormal bases of C^d), with
H = non-orthogonality. Checks:
  * K4-free (exact). Then the union of the four bases is KS-uncolourable: a colouring picks one ray per basis,
    pairwise non-orthogonal, i.e. a K4 of H. The 0-1 colouring of the orthogonality graph is also run as SAT.
  * the combinatorial conditions of §3.12 ((DEG), (C2), (C4), (TRI)), for information.
  * local theta at one ray u per basis, with an exact rational dual certificate. If some u has
    theta(N(u)) < 3, there is no realization: completeness of the three other bases gives
    sum_{c in N(u)} |<u|c>|^2 = 3, while the Lovász bound for the orthonormal representation N(u) with handle u
    gives at most theta.

Examples:
  shrikhande  Cayley graph on Z4 x Z4, connection set {±(1,0), ±(0,1), ±(1,1)}, bases = cosets of 2Z4 x 2Z4
              (d = 4, 16 rays; every local graph is a 6-cycle, theta = 2).
  p41         the multiplier family F(41, 4) (§3.12): the SAT model with the largest local theta among all 113
              K4-free models satisfying (DEG), (C2), (C4), (TRI) and the quotient constraints (d = 41, 164 rays,
              40 neighbours per ray; theta = 2.7924 at every ray). Record: ks/unrealizable_p41.json.

usage: python3 unrealizable_examples.py [shrikhande] [p41] [path/to/record.json ...]
"""
import json
import sys
from fractions import Fraction

import cvxpy as cp
import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from cayley4 import check_conditions, colourable_complement  # noqa: E402
from realize_search import load_pattern  # noqa: E402


def shrikhande():
    conn = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    pts = sorted(((x, y) for x in range(4) for y in range(4)), key=lambda v: (v[0] % 2, v[1] % 2, v))
    X = np.array([[((v[0] - u[0]) % 4, (v[1] - u[1]) % 4) in conn for v in pts] for u in pts])
    return X


def k4_free(X, d):
    part = np.arange(len(X)) // d
    for v in range(d):
        for w in np.where(X[v] & (part == 1))[0]:
            for u in np.where(X[v] & X[w] & (part == 2))[0]:
                if (X[v] & X[w] & X[u] & (part == 3)).any():
                    return False
    return True


def local_theta(X, d, u):
    """theta of the exclusivity graph of N(u) (orthogonal = same basis or non-adjacent), primal and dual."""
    part = np.arange(len(X)) // d
    N = np.where(X[u])[0]
    m = len(N)
    orth = [(i, j) for i in range(m) for j in range(i + 1, m) if part[N[i]] == part[N[j]] or not X[N[i], N[j]]]
    lam = cp.Variable()
    y = cp.Variable(len(orth))
    E = np.zeros((m * m, len(orth)))
    for k, (i, j) in enumerate(orth):
        E[i * m + j, k] = E[j * m + i, k] = 1
    Z = lam * np.eye(m) - np.ones((m, m)) - cp.reshape(E @ y, (m, m), order="C")
    cp.Problem(cp.Minimize(lam), [(Z + Z.T) / 2 >> 0]).solve(solver="CLARABEL")
    return m, orth, float(lam.value), np.array(y.value)


def certify(m, orth, lam, y, slack=Fraction(1, 10 ** 4), digits=10 ** 9):
    """Exact check that (lam + slack) I - J - Y is positive definite with Y rationalised: an LDL^T with
    Fractions whose pivots are all positive. Returns the certified upper bound on theta, or None."""
    bound = Fraction(round(lam * digits), digits) + slack
    M = [[(bound if i == j else 0) - 1 for j in range(m)] for i in range(m)]
    for (i, j), v in zip(orth, y):
        q = Fraction(round(v * digits), digits)
        M[i][j] -= q
        M[j][i] -= q
    for k in range(m):                                   # in-place LDL^T (Gaussian elimination, symmetric)
        piv = M[k][k]
        if piv <= 0:
            return None
        for i in range(k + 1, m):
            if M[i][k] == 0:
                continue
            f = M[i][k] / piv
            for j in range(k + 1, m):
                if M[k][j] != 0:
                    M[i][j] -= f * M[k][j]
    return bound


def report(name, X, d):
    print(f"{name}: d = {d}, {len(X)} rays, degrees {sorted(set(X.sum(1).tolist()))}", flush=True)
    print(f"  K4-free: {k4_free(X, d)};  orthogonality graph 0-1 colourable (SAT): {colourable_complement(X, d)}",
          flush=True)
    issues = check_conditions(X, d)
    kinds = sorted({i[0] for i in issues})
    print(f"  combinatorial conditions (DEG, C2, C4, TRI): {'all hold' if not issues else 'fail: ' + ', '.join(kinds)}",
          flush=True)
    for k in range(4):
        u = k * d
        m, orth, lam, y = local_theta(X, d, u)
        b = certify(m, orth, lam, y)
        print(f"  ray {u} (basis {k}): |N(u)| = {m}, theta = {lam:.6f}; exact certificate theta <= "
              f"{float(b) if b is not None else 'FAILED'}{'  < 3: no realization' if b is not None and b < 3 else ''}",
              flush=True)


if __name__ == "__main__":
    here = __file__.rsplit("/", 1)[0]
    for arg in sys.argv[1:] or ["shrikhande", "p41"]:
        if arg == "shrikhande":
            report("Shrikhande pattern", shrikhande(), 4)
        else:
            path = f"{here}/unrealizable_p41.json" if arg == "p41" else arg
            rec = json.load(open(path))
            report(path.rsplit("/", 1)[-1], load_pattern(path), rec["p"])

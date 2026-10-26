# RESULTS.md 2, Theorem 21: a context with r rigid events (probability exactly 1/9) forces many events in the
# other contexts.  For rigid a, a' (same context A) and the events b of another context B (weights c_b):
#   (P1) sum_{b ~ a} c_b = 1/3,   (P2) sum_{b ~ a, a'} c_b = 1/12,   and every c_b <= 1/9.
# With k_b = number of rigid neighbours of b, (P1)/(P2) fix sum c_b k_b = r/3 and sum c_b k_b(k_b-1) = r(r-1)/12,
# and (P2) gives c_b <= 1/12 whenever k_b >= 2.  Hence |B| >= 12 - 3W, W = weight on events with k_b <= 1.
# (a) exact LP: max W over distributions of k in {0..r} (basic solutions, rational arithmetic);
# (b) MILP over the actual neighbourhood sets S of rigid events: minimum number of events of B.
import itertools
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

def lp_bound(r):
    if r <= 1: return Fr(1)                                     # every event has at most one rigid neighbour
    best = None
    for ks in itertools.combinations(range(r + 1), 3):          # basic feasible solutions have <= 3 nonzero weights
        A = [[Fr(1)] * 3, [Fr(k) for k in ks], [Fr(k * (k - 1)) for k in ks]]
        b = [Fr(1), Fr(r, 3), Fr(r * (r - 1), 12)]
        # solve 3x3 by Cramer's rule
        def det(M): return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        D = det(A)
        if D == 0: continue
        w = []
        for j in range(3):
            Mj = [row[:] for row in A]
            for i in range(3): Mj[i][j] = b[i]
            w.append(det(Mj) / D)
        if min(w) < 0: continue
        W = sum(wj for wj, k in zip(w, ks) if k <= 1)
        if best is None or W > best: best = W
    return best

def milp_bound(r):
    pts = range(r); subsets = [frozenset(s) for m in range(r + 1) for s in itertools.combinations(pts, m)]
    m = len(subsets); cap = [Fr(1, 12) if len(S) >= 2 else Fr(1, 9) for S in subsets]
    # variables: x_S (weights) then n_S (integer counts)
    rows, lo, hi = [], [], []
    def eq(coef_x, val):
        rows.append(list(coef_x) + [0] * m); lo.append(val); hi.append(val)
    eq([1] * m, 1)
    for a in pts: eq([1 if a in S else 0 for S in subsets], 1 / 3)
    for a, a2 in itertools.combinations(pts, 2): eq([1 if (a in S and a2 in S) else 0 for S in subsets], 1 / 12)
    for j in range(m):                                          # x_S <= cap_S * n_S
        row = [0] * (2 * m); row[j] = 1; row[m + j] = -float(cap[j]); rows.append(row); lo.append(-np.inf); hi.append(0)
    c = np.array([0] * m + [1] * m, float)
    res = milp(c, constraints=LinearConstraint(np.array(rows, float), lo, hi),
               integrality=np.array([0] * m + [1] * m), bounds=Bounds(0, np.inf), options={"time_limit": 600})
    return res.status, (round(res.fun) if res.x is not None else None)

if __name__ == "__main__":
    print(" r | max W (exact LP) | |B| >= 12 - 3W | MILP minimum of |B| (status 0 = optimal)")
    for r in range(1, 10):
        W = lp_bound(r); lb = 12 - 3 * W; ceil = -((-lb.numerator) // lb.denominator)
        st, val = milp_bound(r)
        print(f"{r:2d} | {str(W):>16} | {str(lb):>6} -> {ceil:2d}   | {val} (status {st})", flush=True)

def h33_sharpness():
    """In H33 any 3 points (rigid events) see exactly 12 B-lines: the bound for r = 3 is attained."""
    import networkx as nx
    H = nx.read_graph6("h33.g6"); P = range(9); B = range(9, 21)
    out = set()
    for T in itertools.combinations(P, 3):
        cls = {}
        for b in B: cls.setdefault(frozenset(a for a in T if H.has_edge(a, b)), []).append(b)
        sizes = tuple(sorted(len(S) for S in cls))
        out.add((tuple(sorted((len(S), len(v)) for S, v in cls.items())), len(B)))
    return out

def identity_check():
    """X2 = (1/4) I + (3/4) v v^dag + (R + J)/2 on context A, for any realization (no rigidity needed);
    v_a = <g_a|psi>, R = <g_a'|E|g_a>, J = <w_a'|(Pi2 - Pi3)|w_a>, w_a = g_a - psi <psi|g_a>."""
    import json
    from fractions import Fraction as F
    from verify_all import tricayley, tricayley_z2z6
    from vertex_identity import realization
    S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
    cases = []
    V, E = tricayley_z2z6(S36); cases.append(("H36", V, E, {e: F(1, 48) for e in E}, {v: F(1, 12) for v in V}))
    d = json.load(open("shr_q13_1_0139_014_1257_cert.json")); V, E = tricayley(13, tuple(tuple(s) for s in d["S"]), removed={(0, 5), (2, 11)})
    Y = {}
    for a, b, y in d["edges"]:
        a, b = tuple(a), tuple(b); Y[(a, b) if (a, b) in E else (b, a)] = F(y)
    cases.append(("H37", V, E, Y, {tuple(v): F(x) for v, x in zip(d["nodes"], d["c"])}))
    for name, V, E, Y, c in cases:
        idx, g, psi = realization(V, E, Y, c); dim = g.shape[1]
        ctx = [[idx[v] for v in V if v[0] == k] for k in range(3)]
        Pk = [g[ix].T @ g[ix].conj() for ix in ctx]                 # projectors onto the context spans
        F_ = sum(Pk); Em = F_ - 1.5 * np.eye(dim) - 1.5 * np.outer(psi, psi.conj())
        GA = g[ctx[0]].T                                            # columns g_a
        v = GA.conj().T @ psi; W = GA - np.outer(psi, v)            # w_a = g_a - psi <psi|g_a>
        X2 = GA.conj().T @ Pk[1] @ GA; R = GA.conj().T @ Em @ GA; J = W.conj().T @ (Pk[1] - Pk[2]) @ W
        rhs = 0.25 * np.eye(len(v)) + 0.75 * np.outer(v, v.conj()) + (R + J) / 2
        print(f"   {name}: max |X2 - [(1/4)I + (3/4)vv* + (R+J)/2]| = {np.abs(X2 - rhs).max():.1e}")

if __name__ == "__main__":
    print("H33 (r = 3): classes (rigid neighbours, number of B-lines) over all point triples:", h33_sharpness())
    print("exact identity for X2 on non-rigid realizations:"); identity_check()

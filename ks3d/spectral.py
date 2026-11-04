"""Graph spectra and the minimum Kochen–Specker problem in dimension 3 (RESULTS.md §7).

G  = orthogonality (exclusivity) graph of a set of rays, Gc = its complement (non-orthogonality
graph), T = orthogonal bases (triangles of G), m = |T|, w_T(v) = number of bases through v.

Theorem 7.1 says that for every set of rays in C^3 with m > 0:
  * theta(Gc) = xi_f(G) = 3, so no (weighted) Hoffman / inertia / p-energy bound exceeds 3;
  * theta(G, w_T) = theta'(G, w_T) = m, while alpha(G, w_T) = m iff the set is 010-colourable.
This script checks those numbers and tabulates the usual spectral parameters of G and Gc on known
KS sets and on non-KS controls.

usage: python3 spectral.py            (writes the tables printed in spectral.log)
"""
import itertools
import json
import os

import cvxpy as cp
import networkx as nx
import numpy as np
from pysat.solvers import Cadical153
from scipy.optimize import Bounds, LinearConstraint, milp

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-9


# ---------------------------------------------------------------------------------------------
# ray sets


def canon(v):
    v = np.asarray(v, dtype=complex)
    v = v / np.linalg.norm(v)
    k = next(i for i in range(len(v)) if abs(v[i]) > TOL)
    return v * (abs(v[k]) / v[k])


def ray_set(vs):
    out, keys = [], set()
    for v in vs:
        c = canon(v)
        key = tuple(np.round(np.concatenate([c.real, c.imag]), 8))
        if key not in keys:
            keys.add(key)
            out.append(c)
    return out


def peres33():
    r = np.sqrt(2)
    pats = [(0, 0, 1), (0, 1, 1), (0, -1, 1), (0, 1, r), (0, -1, r),
            (1, 1, r), (1, -1, r), (-1, 1, r), (-1, -1, r)]
    return ray_set(p for t in pats for p in itertools.permutations(t))


def cross_closure(rays, rounds=1):
    """Add the third vector of every orthogonal pair (cross product), `rounds` times."""
    for _ in range(rounds):
        new = list(rays)
        for a, b in itertools.combinations(rays, 2):
            if abs(np.vdot(a, b)) < TOL:
                new.append(np.conj(np.cross(a, b)))
        rays = ray_set(new)
    return rays


def known_sets():
    js = json.load(open(os.path.join(HERE, "theta_sets.json")))
    yo = ray_set(js["yu_oh_13"])
    return {
        "KS 31 (height 2)": ray_set(js["min_height2"]),
        "KS 31 (height 3)": ray_set(js["min_height3"]),
        "KS 33 (Peres)": peres33(),
        "Yu-Oh 13 (SI-C)": yo,
        "Yu-Oh closure 25": cross_closure(yo),
    }


# ---------------------------------------------------------------------------------------------
# graph data


class OG:
    """Orthogonality graph of a list of unit vectors."""

    def __init__(self, rays):
        self.rays = [np.asarray(r, dtype=complex) for r in rays]
        V = np.array(self.rays)
        self.gram = V.conj() @ V.T
        n = self.n = len(rays)
        self.A = (np.abs(self.gram) < TOL).astype(int)
        np.fill_diagonal(self.A, 0)
        self.E = [(i, j) for i in range(n) for j in range(i + 1, n) if self.A[i, j]]
        self.T = [t for t in itertools.combinations(range(n), 3)
                  if self.A[t[0], t[1]] and self.A[t[0], t[2]] and self.A[t[1], t[2]]]
        self.wT = np.array([sum(v in t for t in self.T) for v in range(n)])
        self.m = len(self.T)


def from_graph(n, E):
    """An OG-like object for an abstract graph (no vectors)."""
    g = OG.__new__(OG)
    g.n, g.rays, g.gram = n, None, None
    g.A = np.zeros((n, n), dtype=int)
    for i, j in E:
        g.A[i, j] = g.A[j, i] = 1
    g.E = sorted((min(i, j), max(i, j)) for i, j in E)
    g.T = [t for t in itertools.combinations(range(n), 3)
           if g.A[t[0], t[1]] and g.A[t[0], t[2]] and g.A[t[1], t[2]]]
    g.wT = np.array([sum(v in t for t in g.T) for v in range(n)])
    g.m = len(g.T)
    return g


# ---------------------------------------------------------------------------------------------
# spectral bounds (adjacency matrix, or any Hermitian weighting B of the same graph)


def spectral_bounds(B, p_list=(0, 2)):
    lam = np.linalg.eigvalsh(B)
    pos, neg = lam[lam > 1e-8], -lam[lam < -1e-8]
    out = {"l1": lam[-1], "ln": lam[0], "n+": len(pos), "n0": len(lam) - len(pos) - len(neg),
           "n-": len(neg), "hoffman": 1 + lam[-1] / -lam[0]}
    for p in p_list:
        ep, en = (pos ** p).sum(), (neg ** p).sum()
        out[f"energy{p}"] = 1 + max(ep / en, en / ep) ** (1 / abs(p - 1))
    return out


# ---------------------------------------------------------------------------------------------
# SDPs


def theta(n, E, w=None, nonneg=False):
    """Weighted Lovász theta (Schrijver's theta' if nonneg): max sum sqrt(w_i w_j) X_ij,
    X psd, tr X = 1, X_ij = 0 on edges."""
    w = np.ones(n) if w is None else np.asarray(w, float)
    s = np.sqrt(w)
    X = cp.Variable((n, n), symmetric=True)
    cons = [X >> 0, cp.trace(X) == 1] + [X[i, j] == 0 for i, j in E]
    if nonneg:
        cons.append(X >= 0)
    prob = cp.Problem(cp.Maximize(cp.sum(cp.multiply(np.outer(s, s), X))), cons)
    prob.solve(solver=cp.CLARABEL)
    return prob.value


def theta_complement(g):
    """theta(Gc) = strict vector chromatic number of G = best weighted Hoffman bound for chi(G)."""
    Ec = [(i, j) for i in range(g.n) for j in range(i + 1, g.n) if not g.A[i, j]]
    return theta(g.n, Ec)


def born_sdp_subgraph(g):
    """Degree-2 'Born' relaxation of an embedding in C^3 (born.py): x_i = normalised traceless part
    of the projector, <x_i,x_j> = (3|<v_i|v_j>|^2 - 1)/2 = -1/2 on edges, x_ia + x_ib <= 1/2 for
    every edge ab and every i. Distinct rays need <x_i,x_j> < 1. Returns the largest s with
    <x_i,x_j> <= 1 - s for all i != j (s <= 0: the graph has no embedding in C^3)."""
    n = g.n
    X = cp.Variable((n, n), symmetric=True)
    s = cp.Variable()
    cons = [X >> 0, cp.diag(X) == 1, s <= 1.5]
    for i in range(n):
        for j in range(i + 1, n):
            cons += [X[i, j] == -0.5] if g.A[i, j] else [X[i, j] >= -0.5, X[i, j] <= 1 - s]
    for a, b in g.E:
        for i in range(n):
            if i not in (a, b):
                cons.append(X[i, a] + X[i, b] <= 0.5)
    prob = cp.Problem(cp.Maximize(s), cons)
    prob.solve(solver=cp.CLARABEL)
    return prob.value if prob.status in ("optimal", "optimal_inaccurate") else -np.inf


# ---------------------------------------------------------------------------------------------
# exact combinatorial parameters


def alpha(g, w=None):
    w = np.ones(g.n) if w is None else np.asarray(w, float)
    if not g.E:
        return w.sum()
    M = np.zeros((len(g.E), g.n))
    for k, (i, j) in enumerate(g.E):
        M[k, i] = M[k, j] = 1
    r = milp(-w, constraints=LinearConstraint(M, -np.inf, 1), integrality=np.ones(g.n),
             bounds=Bounds(0, 1))
    return round(-r.fun)


def colourable_010(g):
    cl = [[-(i + 1), -(j + 1)] for i, j in g.E] + [[a + 1, b + 1, c + 1] for a, b, c in g.T]
    with Cadical153(bootstrap_with=cl) as s:
        return s.solve()


def k_colourable(g, k):
    var = lambda v, c: v * k + c + 1
    cl = [[var(v, c) for c in range(k)] for v in range(g.n)]
    cl += [[-var(i, c), -var(j, c)] for i, j in g.E for c in range(k)]
    cl.append([var(0, 0)])                       # symmetry: vertex 0 gets colour 0
    with Cadical153(bootstrap_with=cl) as s:
        return s.solve()


def chi_f(g):
    """Fractional chromatic number: LP over all maximal independent sets."""
    Gc = nx.complement(nx.from_numpy_array(g.A))
    I = list(nx.find_cliques(Gc))
    M = np.zeros((g.n, len(I)))
    for k, S in enumerate(I):
        M[S, k] = 1
    r = milp(np.ones(len(I)), constraints=LinearConstraint(M, 1, np.inf), bounds=Bounds(0, np.inf))
    return r.fun, len(I)


def clique_cover(g):
    """chi(Gc): fewest cliques of G (bases, orthogonal pairs, single rays) covering all rays."""
    cl = [list(t) for t in g.T] + [list(e) for e in g.E] + [[v] for v in range(g.n)]
    M = np.zeros((g.n, len(cl)))
    for k, S in enumerate(cl):
        M[S, k] = 1
    r = milp(np.ones(len(cl)), constraints=LinearConstraint(M, 1, np.inf),
             integrality=np.ones(len(cl)), bounds=Bounds(0, 1))
    return round(r.fun)


# ---------------------------------------------------------------------------------------------
# reports


def full_report(name, g, sdp=True):
    sb = spectral_bounds(g.A.astype(float))
    Ac = 1 - np.eye(g.n) - g.A
    sc = spectral_bounds(Ac)
    S = sum(np.outer(v, v.conj()) for v in g.rays) if g.rays is not None else None
    row = dict(name=name, n=g.n, edges=len(g.E), m=g.m, mindeg=int(g.A.sum(1).min()),
               maxdeg=int(g.A.sum(1).max()), KS=not colourable_010(g),
               chi3=k_colourable(g, 3), alpha=alpha(g), alpha_wT=alpha(g, g.wT),
               chif=chi_f(g)[0], cover=clique_cover(g), G=sb, Gc=sc)
    if S is not None:
        row["frame"] = np.linalg.eigvalsh(S)
    if sdp:
        row["theta_Gc"] = theta_complement(g)
        row["theta_G"] = theta(g.n, g.E)
        row["theta_wT"] = theta(g.n, g.E, g.wT)
        row["theta1_wT"] = theta(g.n, g.E, g.wT, nonneg=True)
    return row


def print_rows(rows):
    print("\n## Combinatorial parameters of G (orthogonality graph)\n")
    print("| set | n | edges | bases m | degree | KS | 3-colourable | alpha | alpha(G,w_T) | chi_f |"
          " clique cover = chi(Gc) |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['name']} | {r['n']} | {r['edges']} | {r['m']} | {r['mindeg']}–{r['maxdeg']} |"
              f" {'yes' if r['KS'] else 'no'} | {'yes' if r['chi3'] else 'no'} | {r['alpha']} |"
              f" {r['alpha_wT']} | {r['chif']:.4f} | {r['cover']} |")
    print("\n## Spectral bounds on chi(G) (all must be <= xi_f(G) = 3 for rays in C^3)\n")
    print("| set | lambda_1 | lambda_n | n+ / n0 / n- | Hoffman | inertia (p=0) | p=2 energy |"
          " theta(Gc) |")
    print("|---|---|---|---|---|---|---|---|")
    for r in rows:
        s = r["G"]
        print(f"| {r['name']} | {s['l1']:.4f} | {s['ln']:.4f} | {s['n+']} / {s['n0']} / {s['n-']} |"
              f" {s['hoffman']:.4f} | {s['energy0']:.4f} | {s['energy2']:.4f} |"
              f" {r.get('theta_Gc', float('nan')):.4f} |")
    print("\n## The complement Gc (non-orthogonality graph) and the theta numbers of G\n")
    print("| set | Gc: lambda_1 | Gc: lambda_n | Gc Hoffman (<= cover) | theta(G) | n/3 |"
          " theta(G,w_T) | theta'(G,w_T) | m | frame operator spectrum |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        s = r["Gc"]
        fr = " ".join(f"{x:.3f}" for x in r["frame"]) if "frame" in r else "-"
        print(f"| {r['name']} | {s['l1']:.4f} | {s['ln']:.4f} | {s['hoffman']:.4f} |"
              f" {r.get('theta_G', float('nan')):.4f} | {r['n'] / 3:.3f} |"
              f" {r.get('theta_wT', float('nan')):.4f} | {r.get('theta1_wT', float('nan')):.4f} |"
              f" {r['m']} | {fr} |")


def random_weightings(g, trials=2000, seed=1):
    """Max over random Hermitian edge weightings of the Hoffman / inertia / p=2 bounds."""
    rng = np.random.default_rng(seed)
    best = {"hoffman": 0, "energy0": 0, "energy2": 0}
    I, J = np.array(g.E).T
    for _ in range(trials):
        B = np.zeros((g.n, g.n), dtype=complex)
        z = rng.normal(size=len(g.E)) + 1j * rng.normal(size=len(g.E)) * rng.integers(0, 2)
        z *= rng.exponential(size=len(g.E)) ** rng.uniform(0, 3)
        B[I, J] = z
        B[J, I] = z.conj()
        sb = spectral_bounds(B)
        for k in best:
            best[k] = max(best[k], sb[k])
    return best


def main():
    sets = known_sets()
    rows = []
    for name, rays in sets.items():
        g = OG(rays)
        rows.append(full_report(name, g))
    g31 = OG(sets["KS 31 (height 2)"])
    for v in (0, 10, 20, 30):
        keep = [r for k, r in enumerate(g31.rays) if k != v]
        rows.append(full_report(f"KS 31 minus ray {v}", OG(keep)))
    print_rows(rows)

    print("\n## Barrier checks\n")
    for r in rows:
        print(f"{r['name']:24s} theta(G,w_T) - m = {r['theta_wT'] - r['m']:+.2e}   "
              f"theta'(G,w_T) - m = {r['theta1_wT'] - r['m']:+.2e}   "
              f"alpha(G,w_T) = {r['alpha_wT']} (m = {r['m']})   theta(Gc) - 3 = {r['theta_Gc'] - 3:+.1e}")
    print("\nRandom Hermitian edge weightings (2000 each): best bounds on chi(G)")
    for name in ("KS 31 (height 2)", "KS 33 (Peres)"):
        best = random_weightings(OG(sets[name]))
        print(f"  {name}: " + ", ".join(f"{k} {v:.4f}" for k, v in best.items()))


if __name__ == "__main__":
    main()

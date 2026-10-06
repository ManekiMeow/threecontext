# Tri-Cayley graphs over a finite group Gamma:
#   vertices (a,g), a in Z3, g in Gamma;  (a,g) ~ (a+1, g*s) for s in S_a.
# Triangle-free  <=>  e not in S0*S1*S2.  Left multiplication by Gamma acts as automorphisms.
# Normalisation: e in S0 and e in S1 (relabel (a,g) -> (a, g t_a)).
import sys, itertools, numpy as np, networkx as nx
from multiprocessing import Pool
from theta import theta_comp


def perm_group(gens):
    gens = [tuple(g) for g in gens]; e = tuple(range(len(gens[0])))
    G = {e}; frontier = [e]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = tuple(x[i] for i in g)
                if y not in G: G.add(y); new.append(y)
        frontier = new
    G = sorted(G); idx = {g: i for i, g in enumerate(G)}
    mul = [[idx[tuple(a[i] for i in b)] for b in G] for a in G]   # a*b = "apply b then a" composition table
    return mul, idx[e]


def groups():
    # name -> (multiplication table, identity)
    out = {}
    out["Z3xZ3"] = perm_group([[1, 2, 0, 3, 4, 5], [0, 1, 2, 4, 5, 3]])
    out["D5"] = perm_group([[1, 2, 3, 4, 0], [0, 4, 3, 2, 1]])
    out["Z2xZ6"] = perm_group([[1, 0, 2, 3, 4, 5, 6], [0, 1, 3, 2, 4, 5, 6], [0, 1, 2, 3, 5, 6, 4]])
    out["D6"] = perm_group([[1, 2, 3, 4, 5, 0], [0, 5, 4, 3, 2, 1]])
    out["A4"] = perm_group([[1, 2, 0, 3], [1, 0, 3, 2]])
    # Dic3 = <x,y | x^6 = 1, y^2 = x^3, y x y^-1 = x^-1> as a regular permutation group on 12 points
    els = [(i, j) for j in range(2) for i in range(6)]

    def m(p, q):
        (i1, j1), (i2, j2) = p, q
        if j1 == 0: i = (i1 + i2) % 6
        else: i = (i1 - i2) % 6
        j = j1 + j2
        if j == 2: j = 0; i = (i + 3) % 6
        return (i, j)
    ix = {g: k for k, g in enumerate(els)}
    out["Dic3"] = ([[ix[m(a, b)] for b in els] for a in els], ix[(0, 0)])
    for name, (mul, e) in out.items():
        n = len(mul)
        assert all(mul[e][x] == x for x in range(n))
        assert all(mul[mul[x][y]][z] == mul[x][mul[y][z]] for x in range(n) for y in range(n) for z in range(n))
    return out


def tricayley(mul, S):
    n = len(mul); H = nx.Graph()
    for a in range(3):
        for g in range(n):
            for s in S[a]: H.add_edge((a, g), ((a + 1) % 3, mul[g][s]))
    return H


def candidates(mul, e, sizes):
    n = len(mul); others = [x for x in range(n) if x != e]
    inv = {x: next(y for y in range(n) if mul[x][y] == e) for x in range(n)}
    seen = set()
    a, b, c = sizes
    for S0r in itertools.combinations(others, a - 1):
        S0 = (e,) + S0r
        for S1r in itertools.combinations(others, b - 1):
            S1 = (e,) + S1r
            forb = {inv[mul[x][y]] for x in S0 for y in S1}
            allowed = [z for z in range(n) if z not in forb]
            for S2 in itertools.combinations(allowed, c):
                H = tricayley(mul, (S0, S1, S2))
                if len(H) < 3 * n: continue
                key = tuple(np.round(np.linalg.eigvalsh(nx.to_numpy_array(H, nodelist=sorted(H))), 6))
                if key in seen: continue
                seen.add(key); yield (S0, S1, S2)


MUL = None


def work(S):
    t, _ = theta_comp(tricayley(MUL, S)); return t, S


if __name__ == "__main__":
    G = groups(); name = sys.argv[1]; MUL, e = G[name]
    n = len(MUL)
    for a in range(1, n + 1):
        for b in range(a, n + 1):
            for c in range(b, n + 1):
                if c > n - b or a + b + c > n + 3: continue
                cands = list(candidates(MUL, e, (a, b, c)))
                if not cands: continue
                with Pool(3) as P: res = P.map(work, cands)
                res.sort(reverse=True)
                for t, S in res:
                    if t > 3 - 1e-6: print("HIT", name, S, t, flush=True)
                print(name, (a, b, c), "tested", len(res), "best", res[0], flush=True)

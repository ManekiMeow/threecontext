"""RESULTS.md Section 2.6: no 30-event GHZ graph has an automorphism group that is transitive on a context.

Reduction (Section 2.6):
  * context-preserving automorphisms act faithfully on each context (Proposition 39 of the four-context branch);
  * a group transitive on context A has no fixed points on the other contexts (Corollary 40 there);
  * every transitive group of degree 10 contains one of six minimal transitive groups
    C10, D5 (regular), F20 (on 10 cosets of an involution), A5 (on 2-subsets), 2^4:5, 5^2:8.
For each minimal group M this script lists every faithful fixed-point-free action of M on B and C (with A transitive),
enumerates all maximal triangle-free M-invariant supports, and computes theta(Hbar); with --certify it also checks
an exact rational upper bound < 3 for one representative of every isomorphism class.

usage: python3 search/transitive_symmetry.py [--certify] [--only NAME]
"""
import itertools
import sys
import time

import numpy as np


# ------------------------------------------------------------------ groups as permutation groups on a base set
def closure(gens):
    n = len(gens[0])
    e = tuple(range(n))
    G, fr = {e}, [e]
    while fr:
        a = fr.pop()
        for g in gens:
            b = tuple(g[a[i]] for i in range(n))       # g o a
            if b not in G:
                G.add(b)
                fr.append(b)
    return sorted(G)


def comp(g, h):                                          # (g o h)(x) = g(h(x))
    return tuple(g[h[i]] for i in range(len(h)))


def subgroup(gens, n):
    return closure(list(gens)) if gens else [tuple(range(n))]


def groups():
    out = {}
    # C10 regular
    out["C10"] = dict(gens=[tuple((x + 1) % 10 for x in range(10))])
    # D5 on Z5 (natural, faithful)
    r = tuple((x + 1) % 5 for x in range(5)); s = tuple((-x) % 5 for x in range(5))
    out["D5"] = dict(gens=[r, s])
    # F20 = AGL(1,5) on Z5
    out["F20"] = dict(gens=[r, tuple((2 * x) % 5 for x in range(5))])
    # A5 on {0..4}
    out["A5"] = dict(gens=[(1, 2, 0, 3, 4), (1, 2, 3, 4, 0)])
    # 2^4:5 on Z5 x Z2 (index 2i+e): translation and an even sign change
    tr = tuple(2 * ((i + 1) % 5) + e for i in range(5) for e in range(2))
    fl = tuple(2 * i + ((e + (i < 2)) % 2) for i in range(5) for e in range(2))
    out["2^4:5"] = dict(gens=[tr, fl])
    # 5^2:8 on Z5 x {0,1} (index 5*eps + i): s1, s2 rotate the two blocks, t: (i,0)->(i,1), (i,1)->(2i,0)
    s1 = tuple(((i + 1) % 5) if eps == 0 else 5 + i for eps in range(2) for i in range(5))
    s2 = tuple(i if eps == 0 else 5 + (i + 1) % 5 for eps in range(2) for i in range(5))
    t = tuple(5 + i if eps == 0 else (2 * i) % 5 for eps in range(2) for i in range(5))
    out["5^2:8"] = dict(gens=[s1, s2, t])
    for name, d in out.items():
        d["elements"] = closure(d["gens"])
        d["n"] = len(d["gens"][0])
    return out


# ------------------------------------------------------------------ M-sets: unions of coset spaces M/S
def coset_action(G, S):
    """points = left cosets gS; returns list of permutations (one per element of G, same order as G)"""
    Sset = set(S)
    cosets, index = [], {}
    for g in G:
        key = frozenset(comp(g, s) for s in S)
        if key not in index:
            index[key] = len(cosets)
            cosets.append(key)
    perms = []
    for g in G:
        perms.append(tuple(index[frozenset(comp(g, x) for x in c)] for c in cosets))
    return perms


def mset(G, stabilizers):
    """disjoint union of coset spaces G/S for S in stabilizers; returns per-element permutations of 10 points"""
    parts = [coset_action(G, S) for S in stabilizers]
    sizes = [len(p[0]) for p in parts]
    off = np.cumsum([0] + sizes)
    perms = []
    for gi in range(len(G)):
        perm = []
        for k, p in enumerate(parts):
            perm += [off[k] + x for x in p[gi]]
        perms.append(tuple(int(x) for x in perm))
    assert sum(sizes) == 10, sizes
    # faithful and fixed-point-free checks
    assert len(set(perms)) == len(G), "not faithful"
    return perms


def action_list(name, G, n):
    """faithful fixed-point-free actions on 10 points needed for the transitive-symmetry reduction"""
    e = tuple(range(n))
    def sub(gens):
        return subgroup(gens, n)
    if name == "C10":
        return {"10": [sub([])]}
    if name == "D5":
        return {"10": [sub([])], "55": [sub([tuple((-x) % 5 for x in range(5))])] * 2}
    if name == "F20":
        return {"10": [sub([tuple((-x) % 5 for x in range(5))])], "55": [sub([tuple((2 * x) % 5 for x in range(5))])] * 2}
    if name == "A5":
        stab_pair = [g for g in G if {g[0], g[1]} == {0, 1}]           # setwise stabilizer of {0,1}: S3
        stab_pt = [g for g in G if g[4] == 4]                          # A4
        return {"10": [stab_pair], "55": [stab_pt] * 2}
    if name == "2^4:5":
        # translation subgroup Z2^4 = elements fixing every block {2i, 2i+1}; functionals f on it
        K = [g for g in G if all(g[2 * i] // 2 == i for i in range(5))]
        def vec(g):                                   # sign-change vector of g in K
            return tuple((g[2 * i] - 2 * i) % 2 for i in range(5))
        acts = {}
        for cls, f in (("10a", (1, 0, 0, 0, 0)), ("10b", (1, 1, 0, 0, 0)), ("10c", (1, 0, 1, 0, 0))):
            S = [g for g in K if sum(a * b for a, b in zip(f, vec(g))) % 2 == 0]
            assert len(S) == 8
            acts[cls] = [S]
        return acts
    if name == "5^2:8":
        P = [g for g in G if all((g[i] < 5) == (i < 5) for i in range(10)) and
             all((g[i] - i) % 5 == (g[0] - 0) % 5 for i in range(5)) and
             all((g[5 + i] - 5 - i) % 5 == (g[5] - 5) % 5 for i in range(5))]
        assert len(P) == 25
        t2 = tuple((2 * i) % 5 if eps == 0 else 5 + (2 * i) % 5 for eps in range(2) for i in range(5))
        def line(a, b):                                # subgroup generated by s1^a s2^b
            g = tuple(((i + a) % 5) if eps == 0 else 5 + (i + b) % 5 for eps in range(2) for i in range(5))
            return g
        acts = {}
        for cls, (a, b) in (("10a", (0, 1)), ("10b", (1, 1)), ("10c", (1, 2))):
            S = sub([line(a, b), t2])
            assert len(S) == 20, len(S)
            acts[cls] = [S]
        return acts
    raise ValueError(name)


# ------------------------------------------------------------------ orbitals, triangle table, maximal triples
def orbitals(PA, PB):
    seen, orbs = {}, []
    for x in range(10):
        for y in range(10):
            if (x, y) in seen:
                continue
            o = {(pa[x], pb[y]) for pa, pb in zip(PA, PB)}
            for p in o:
                seen[p] = len(orbs)
            orbs.append(o)
    lookup = np.zeros((10, 10), int)
    for (x, y), i in seen.items():
        lookup[x, y] = i
    return orbs, lookup


def maximal_triples(PA, PB, PC):
    oAB, lAB = orbitals(PA, PB)
    oBC, lBC = orbitals(PB, PC)
    oAC, lAC = orbitals(PA, PC)
    p, q, r = len(oAB), len(oBC), len(oAC)
    T = np.zeros((p, q), dtype=np.int64)            # bitmask of AC-orbitals closing a triangle
    for i, o in enumerate(oAB):
        a0, b0 = next(iter(o))
        for j in range(q):
            m = 0
            for c in range(10):
                if lBC[b0, c] == j:
                    m |= 1 << int(lAC[a0, c])
            T[i, j] = m
    fullr = (1 << r) - 1
    nq = 1 << q
    masks = np.arange(nq, dtype=np.int64)
    bit = [((masks >> j) & 1).astype(bool) for j in range(q)]
    G = np.zeros((p, nq), dtype=np.int64)           # G[i][S_BC] = OR_{j in S_BC} T[i, j]
    for i in range(p):
        for j in range(q):
            G[i] |= np.where(bit[j], T[i, j], 0)
    out = []
    for SAB in range(1, 1 << p):
        U = np.zeros(q, dtype=np.int64)
        for i in range(p):
            if SAB >> i & 1:
                U |= T[i]
        forb = np.zeros(nq, dtype=np.int64)
        for i in range(p):
            if SAB >> i & 1:
                forb |= G[i]
        SAC = fullr & ~forb
        ok = (SAC != 0) & (masks != 0)
        for i in range(p):                           # maximality of S_AB
            if not SAB >> i & 1:
                ok &= (G[i] & SAC) != 0
        for j in range(q):                           # maximality of S_BC
            ok &= bit[j] | ((U[j] & SAC) != 0)
        for SBC in np.nonzero(ok)[0]:
            out.append((SAB, int(SBC), int(SAC[SBC])))
    def edges(t3):
        SAB, SBC, SAC = t3
        E = set()
        for S, orbs, k, l in ((SAB, oAB, 0, 1), (SBC, oBC, 1, 2), (SAC, oAC, 0, 2)):
            for i, o in enumerate(orbs):
                if S >> i & 1:
                    for x, y in o:
                        E.add((10 * k + x, 10 * l + y))
        return sorted(E)
    return out, (oAB, oBC, oAC), edges


# ------------------------------------------------------------------ theta and certificates
def theta_solver(N=30):
    import cvxpy as cp
    B = cp.Variable((N, N), symmetric=True)
    mask = cp.Parameter((N, N))
    prob = cp.Problem(cp.Maximize(cp.sum(B)), [B >> 0, cp.trace(B) == 1, cp.multiply(mask, B) == 0])
    def th(E):
        M = np.ones((N, N))
        for u, v in E:
            M[u, v] = M[v, u] = 0
        np.fill_diagonal(M, 0)
        mask.value = M
        prob.solve(solver="CLARABEL")
        return prob.value
    return th


_AUT_CACHE = {}


def automorphisms(G, gens=None):
    """automorphisms of the abstract group G (elements as permutations), as dicts element -> element"""
    key = (len(G), G[0], G[-1])
    if key not in _AUT_CACHE:
        _AUT_CACHE[key] = _automorphisms(G, gens)
    return _AUT_CACHE[key]


def _automorphisms(G, given=None):
    gens = list(given) if given else []
    # pick a small generating set greedily
    H = [tuple(range(len(G[0])))]
    if gens:
        assert len(closure(gens)) == len(G)
    for g in ([] if gens else G):
        if g not in set(H):
            gens.append(g)
            H = closure(gens)
            if len(H) == len(G):
                break
    idx = {g: i for i, g in enumerate(G)}
    auts = []
    for imgs in itertools.product(G, repeat=len(gens)):
        f = {tuple(range(len(G[0]))): tuple(range(len(G[0])))}
        fr = [tuple(range(len(G[0])))]
        ok = True
        while fr and ok:
            x = fr.pop()
            for g, im in zip(gens, imgs):
                y, fy = comp(g, x), comp(im, f[x])
                if y in f:
                    if f[y] != fy:
                        ok = False
                        break
                else:
                    f[y] = fy
                    fr.append(y)
        if ok and len(f) == len(G) and len(set(f.values())) == len(G):
            auts.append(f)
    return auts


def intertwiners(PX, PX2):
    """bijections psi of 10 points with psi o PX[g] = PX2[g] o psi for all g (PX, PX2: per-element permutations)"""
    orbits, seen = [], set()
    for x in range(10):
        if x not in seen:
            o = {pg[x] for pg in PX}
            seen |= o
            orbits.append(min(o))
    out = []
    for imgs in itertools.product(range(10), repeat=len(orbits)):
        psi = [None] * 10
        ok = True
        for x0, y0 in zip(orbits, imgs):
            for pg, qg in zip(PX, PX2):
                x, y = pg[x0], qg[y0]
                if psi[x] is None:
                    psi[x] = y
                elif psi[x] != y:
                    ok = False
                    break
            if not ok:
                break
        if ok and None not in psi and len(set(psi)) == 10:
            out.append(tuple(psi))
    return out


def symmetry_group(G, PA, PB, PC, orbs, gens=None):
    """induced permutations of (AB, BC, AC) orbital indices from triples (psi_A, psi_B, psi_C) intertwining the
    actions twisted by a common automorphism of G"""
    oAB, oBC, oAC = orbs
    def index(orbits):
        d = {}
        for i, o in enumerate(orbits):
            for pr in o:
                d[pr] = i
        return d
    iAB, iBC, iAC = index(oAB), index(oBC), index(oAC)
    gidx = {g: i for i, g in enumerate(G)}
    perms = set()
    for f in automorphisms(G, gens):
        tw = [gidx[f[g]] for g in G]                  # element index -> index of its image
        PA2 = [PA[tw[i]] for i in range(len(G))]
        PB2 = [PB[tw[i]] for i in range(len(G))]
        PC2 = [PC[tw[i]] for i in range(len(G))]
        for a in intertwiners(PA, PA2):
            for b in intertwiners(PB, PB2):
                for c in intertwiners(PC, PC2):
                    pAB = tuple(iAB[(a[x], b[y])] for (x, y) in (next(iter(o)) for o in oAB))
                    pBC = tuple(iBC[(b[x], c[y])] for (x, y) in (next(iter(o)) for o in oBC))
                    pAC = tuple(iAC[(a[x], c[y])] for (x, y) in (next(iter(o)) for o in oAC))
                    perms.add((pAB, pBC, pAC))
    return sorted(perms)


def apply_mask(perm, m):
    out = 0
    i = 0
    while m:
        if m & 1:
            out |= 1 << perm[i]
        m >>= 1
        i += 1
    return out


def classes(triples, sym):
    seen, reps = set(), []
    for t3 in triples:
        if t3 in seen:
            continue
        orbit = {(apply_mask(pa, t3[0]), apply_mask(pb, t3[1]), apply_mask(pc, t3[2])) for pa, pb, pc in sym}
        seen |= orbit
        reps.append(t3)
    return reps


def main():
    certify = "--certify" in sys.argv
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    half = "--half" in sys.argv                       # also the intransitive (5,5)^3 actions of D5, F20, A5
    if certify:
        import networkx as nx
        from dualcert import exact_upper
    th = theta_solver()
    t0 = time.time()
    allok = True
    rows = []
    for name, d in groups().items():
        if only and name != only:
            continue
        G, n = d["elements"], d["n"]
        acts = action_list(name, G, n)
        Ps = {k: mset(G, v) for k, v in acts.items()}
        trans = [k for k in Ps if k.startswith("10")]
        combos = sorted({(a,) + tuple(sorted((b, c))) for a in trans for b in Ps for c in Ps})
        if half and "55" in Ps:
            combos.append(("55", "55", "55"))
        for combo in combos:
            PA, PB, PC = (Ps[k] for k in combo)
            triples, orbs, edges = maximal_triples(PA, PB, PC)
            if len(triples) > 300 and len(G) <= 80:
                sym = symmetry_group(G, PA, PB, PC, orbs, d["gens"])
            else:                                     # small families: no reduction needed
                sym = [(tuple(range(len(orbs[0]))), tuple(range(len(orbs[1]))), tuple(range(len(orbs[2]))))]
            reps = classes(triples, sym)
            vals = [th(edges(t3)) for t3 in reps]
            mx = max(vals) if vals else float("nan")
            line = (f"{name:6s} actions {combo}: orbitals {tuple(len(o) for o in orbs)}, maximal supports "
                    f"{len(triples)}, symmetry classes {len(reps)} (group of {len(sym)}), max theta {mx:.6f}")
            worst = None
            if certify:
                worst = 0
                for t3 in reps:
                    H = nx.Graph(); H.add_nodes_from(range(30)); H.add_edges_from(edges(t3))
                    bound, _ = exact_upper(H)
                    if bound is None or bound >= 3:
                        allok = False; worst = None; break
                    worst = max(worst, bound)
                line += f"; exact bounds < 3: {worst is not None}, largest {worst}"
            rows.append((name, combo, len(triples), len(reps), mx, worst))
            print(line, f"({time.time() - t0:.0f}s)", flush=True)
    if certify:
        print("ALL CERTIFIED" if allok else "SOME CLASS NOT CERTIFIED", flush=True)


if __name__ == "__main__":
    main()

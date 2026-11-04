"""Which single-class multiplier patterns are K4-free? (RESULTS.md §3.12, Theorem 23 and Proposition 24.)

Pattern: every S_kl (k < l) is one class C_j = g^j H_e of the e-th powers H_e in Z_p^* (g a primitive root), so
S_lk = -S_kl = C_{j + (p-1)/2}. A transversal K4 through (3, 0) is a triple a ∈ S_30, b ∈ S_31, c ∈ S_32 with
b - a ∈ S_01, c - a ∈ S_02, c - b ∈ S_12. Multiplying by H_e preserves everything, so a can be fixed to one
representative of S_30.

Proposition 24 says that every pattern contains such a K4 once p - (e-1)^2 (2e+1) sqrt(p) - 3e^2 > 0, so the scan up
to that threshold is complete. Union patterns (S_kl a union of classes, possibly with 0) contain a single-class
sub-pattern, so they can be K4-free only at primes where a single-class pattern is.

usage: python3 multiplier_classify.py E [PMAX]      (default PMAX = the Weil threshold for E)
"""
import itertools
import sys
import time

import numpy as np


def threshold(e):
    """Smallest p0 with p - (e-1)^2 (2e+1) sqrt(p) - 3e^2 > 0 for all p >= p0."""
    c = (e - 1) ** 2 * (2 * e + 1)
    s = (c + np.sqrt(c * c + 12 * e * e)) / 2
    return int(np.floor(s * s)) + 1


def primitive_root(p):
    fac = [q for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, int(q ** 0.5) + 1))]
    return next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in fac))


def k4_free_patterns(p, e):
    g = primitive_root(p)
    cls = np.full(p, -1)
    x = 1
    for i in range(p - 1):
        cls[x] = i % e
        x = x * g % p
    members = [np.where(cls == j)[0] for j in range(e)]
    neg = (p - 1) // 2 % e                     # -C_j = C_{j + neg}
    out = []
    for pat in itertools.product(range(e), repeat=6):
        if pat[0] != 0:                        # rotate classes so that S_01 = C_0 (multiplying by g)
            continue
        j01, j02, j03, j12, j13, j23 = pat
        a = int(members[(j03 + neg) % e][0])   # a ∈ S_30 = -S_03, one representative
        B = members[(j13 + neg) % e]
        B = B[cls[(B - a) % p] == j01]          # b ∈ S_31, b - a ∈ S_01
        C = members[(j23 + neg) % e]
        C = C[cls[(C - a) % p] == j02]          # c ∈ S_32, c - a ∈ S_02
        if len(B) and len(C) and (cls[(C[None, :] - B[:, None]) % p] == j12).any():
            continue
        out.append(pat)
    return out


def k4_free_patterns_fast(p, e, chunk=256):
    """Same result as k4_free_patterns, for large p: the candidate sets B (b ∈ S_31, b - a ∈ S_01) and C
    (c ∈ S_32, c - a ∈ S_02) depend on three class indices each and are precomputed; pairs (b, c) are tested in growing
    chunks and the test stops at the first K4, so only K4-free patterns cost a full |B| x |C| check."""
    g = primitive_root(p)
    cls = np.full(p, -1)
    x = 1
    for i in range(p - 1):
        cls[x] = i % e
        x = x * g % p
    members = [np.where(cls == j)[0] for j in range(e)]
    neg = (p - 1) // 2 % e
    reps = [int(members[j][0]) for j in range(e)]
    # sets[(a_class, z_class, diff_class)] = {z ∈ C_z : z - a ∈ C_diff}, with a = reps[a_class]
    sets = {}
    for ac in range(e):
        a = reps[ac]
        for zc in range(e):
            Z = members[zc]
            dz = cls[(Z - a) % p]
            for dc in range(e):
                sets[(ac, zc, dc)] = Z[dz == dc]
    out = []
    for pat in itertools.product(range(e), repeat=6):
        if pat[0] != 0:
            continue
        j01, j02, j03, j12, j13, j23 = pat
        ac = (j03 + neg) % e
        B = sets[(ac, (j13 + neg) % e, j01)]
        C = sets[(ac, (j23 + neg) % e, j02)]
        found = False
        if len(B) and len(C):
            for size in (64, 512, None):                  # small corner first: almost every pattern has a K4 there
                Bc, Cc = B[:size], C[:size]
                found = bool((cls[(Cc[None, :] - Bc[:, None]) % p] == j12).any())
                if found or size is None or (size >= len(B) and size >= len(C)):
                    break
            if not found and (len(B) > 512 or len(C) > 512):
                for lo in range(0, len(B), chunk):        # exhaustive check, in row blocks to bound memory
                    if (cls[(C[None, :] - B[lo:lo + chunk, None]) % p] == j12).any():
                        found = True
                        break
        if not found:
            out.append(pat)
    return out


def theta_orbits(p, e, sets, k):
    """theta of the exclusivity graph of N((k, 0)) for an H_e-invariant pattern (sets[(a, b)] for all ordered pairs).
    H_e = <g^e> ≅ Z_h acts on N((k, 0)) by multiplication: each class C_j inside S_kl is a free orbit (vertex
    c g^{e t} <-> (orbit, t)) and 0 ∈ S_kl is a fixed point. An invariant optimal X is circulant between free orbits
    and constant between a free orbit and a fixed point, so the PSD condition splits into the Fourier blocks
    X̂(w), w = 0..h/2 (free orbits only for w > 0)."""
    import cvxpy as cp
    g = primitive_root(p)
    h = (p - 1) // e
    gp = [pow(g, e * t, p) for t in range(h)]
    cls = {}
    x = 1
    for i in range(p - 1):
        cls[x] = i % e
        x = x * g % p
    others = [l for l in range(4) if l != k]
    orbits, fixed = [], []                                # orbits: (part, representative); fixed: part
    for l in others:
        S = sets[(k, l)]
        for j in sorted({cls[s] for s in S if s}):
            orbits.append((l, next(s for s in sorted(S) if s and cls[s] == j)))
        if 0 in S:
            fixed.append(l)
    adj = lambda l1, x1, l2, x2: l1 != l2 and (x2 - x1) % p in sets[(l1, l2)]
    q, r = len(orbits), len(fixed)
    xv = {}                                               # circulant first rows x_ab(delta), allowed delta only
    for a in range(q):
        for b in range(q):
            la, ca = orbits[a]
            lb, cb = orbits[b]
            if a == b:
                xv[(a, b)] = [0]
            else:
                xv[(a, b)] = [d for d in range(h) if adj(la, ca, lb, cb * gp[d] % p)]
    var = {key: cp.Variable(len(ds)) for key, ds in xv.items() if ds}
    c = {(a, i): cp.Variable() for a in range(q) for i in range(r)
         if adj(orbits[a][0], orbits[a][1], fixed[i], 0)}
    xf = cp.Variable((r, r), symmetric=True) if r else None
    cons = []
    for a in range(q):
        for b in range(a + 1, q):                          # X real symmetric: x_ba(d) = x_ab(-d)
            if (a, b) in var:
                db = {d: i for i, d in enumerate(xv[(b, a)])}
                cons += [var[(b, a)][db[(-d) % h]] == var[(a, b)][i] for i, d in enumerate(xv[(a, b)])]
    for i in range(r):
        for j in range(r):
            if i != j and not adj(fixed[i], 0, fixed[j], 0):
                cons.append(xf[i, j] == 0)
    hat = lambda a, b, w: (cp.sum(cp.multiply(np.exp(2j * np.pi * w * np.array(xv[(a, b)]) / h), var[(a, b)]))
                           if (a, b) in var else 0)
    for w in range(h // 2 + 1):
        n_blk = q + (r if w == 0 else 0)
        M = cp.Variable((n_blk, n_blk), hermitian=True)
        cons.append(M >> 0)
        for a in range(q):
            for b in range(q):
                cons.append(M[a, b] == hat(a, b, w))
        if w == 0:
            for a in range(q):
                for i in range(r):
                    cons.append(M[a, q + i] == (np.sqrt(h) * c[(a, i)] if (a, i) in c else 0))
            for i in range(r):
                for j in range(r):
                    cons.append(M[q + i, q + j] == xf[i, j])
    tr = h * sum(var[(a, a)][0] for a in range(q)) + (cp.trace(xf) if r else 0)
    obj = h * sum(cp.sum(v) for v in var.values()) + 2 * h * sum(c.values()) + (cp.sum(xf) if r else 0)
    cons.append(tr == 1)
    prob = cp.Problem(cp.Maximize(obj), cons)
    try:
        prob.solve(solver="CLARABEL")
    except cp.error.SolverError:                          # rare numerical failure: fall back to SCS
        prob.solve(solver="SCS", eps_abs=1e-9, eps_rel=1e-9, max_iters=200000)
    return prob.value


def union_thetas(p, e, T, z):
    """theta at all four bases for the union pattern with class sets T (pairs in order 01,02,03,12,13,23) and zero
    flags z."""
    import cyclo4
    C = cyclo4.classes(p, e)
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    sets = {}
    for (k, l), Tk, zk in zip(pairs, T, z):
        S = set(int(x) for j in Tk for x in C[j]) | ({0} if zk else set())
        sets[(k, l)] = S
        sets[(l, k)] = {(-x) % p for x in S}
    out = []
    for k in (3, 0, 1, 2):
        th = theta_orbits(p, e, sets, k)
        out.append(round(th, 6))
        if th < 3 - 1e-6:
            break
    return out


def full_check(p, e, pat):
    """The remaining necessary conditions for a single-class pattern: cyclo4.check ((DEG), (K4), (C2), (C4), (TRI))
    and theta at the neighbourhood of (k, 0) for all four bases. A neighbourhood is a tri-Cayley graph over
    H_e = <g^e> ≅ Z_h (index vertices of the coset x H_e by the exponent t of x g^{et}), so theta splits into
    h Hermitian 3x3 blocks (search/ftheta.py)."""
    import cyclo4
    sys.path.insert(0, __file__.rsplit("/", 2)[0] + "/search")
    from ftheta import ftheta
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    C = cyclo4.classes(p, e)                          # class i = g^i H_e with g the smallest primitive root
    S = {pr: np.isin(np.arange(p), C[c]) for pr, c in zip(pairs, pat)}
    why = cyclo4.check(p, S)
    if why is not None:
        return why, None
    g = cyclo4.primitive_root(p)
    h = (p - 1) // e
    sets = {}
    for (k, l), c in zip(pairs, pat):
        sets[(k, l)] = set(C[c])
        sets[(l, k)] = {(-x) % p for x in C[c]}
    powers = [pow(g, e * t, p) for t in range(h)]
    thetas = []
    for k in range(4):
        l1, l2, l3 = [l for l in range(4) if l != k]
        rep = {l: min(sets[(k, l)]) for l in (l1, l2, l3)}
        D = lambda a, b: [t for t in range(h) if (rep[b] * powers[t] - rep[a]) % p in sets[(a, b)]]
        thetas.append(round(ftheta(h, (D(l1, l2), D(l2, l3), D(l3, l1))), 6))
    return ("theta<3" if min(thetas) < 3 - 1e-6 else None), thetas


def union_survivors(p, e, log=print):
    """All union patterns S_kl = (union of classes T_kl) [+ {0}] passing cyclo4.check. A union is K4-free only if every
    single-class choice in T_01 x ... x T_23 is, so the T's form a box inside the K4-free single-class set F.
    Returns the survivors as (T tuple, zero flags), up to rotation of the classes."""
    import cyclo4
    canon = k4_free_patterns(p, e)
    F = {tuple((x + r) % e for x in pat) for pat in canon for r in range(e)}
    prefixes = [set(f[:m] for f in F) for m in range(7)]
    subsets = [frozenset(c) for r in range(1, e + 1) for c in itertools.combinations(range(e), r)]
    boxes = []

    def extend(Ts):
        m = len(Ts)
        if m == 6:
            boxes.append(tuple(Ts))
            return
        for T in subsets:
            if all(tuple(c) in prefixes[m + 1] for c in itertools.product(*(Ts + [T]))):
                extend(Ts + [T])
    extend([])
    seen, out = set(), []
    C = cyclo4.classes(p, e)
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for box in boxes:
        key = min(tuple(tuple(sorted((x + r) % e for x in T)) for T in box) for r in range(e))
        if key in seen:
            continue
        seen.add(key)
        for zflags in itertools.product((0, 1), repeat=6):
            S = {}
            for pr, T, z in zip(pairs, key, zflags):
                v = np.isin(np.arange(p), np.concatenate([C[j] for j in T]))
                v[0] = bool(z)
                S[pr] = v
            if cyclo4.check(p, S) is None:
                out.append((key, zflags))
    log(f"p={p}: {len(boxes)} boxes ({len(seen)} up to rotation), {len(out)} union patterns pass cyclo4.check")
    return out


def main():
    if sys.argv[1] == "unions":                       # python3 multiplier_classify.py unions E P...
        e = int(sys.argv[2])
        for p in map(int, sys.argv[3:]):
            for s in union_survivors(p, e):
                print("   ", s, flush=True)
        return
    if sys.argv[1] == "check":                        # python3 multiplier_classify.py check E P PATTERN(comma-separated)
        e, p, pat = int(sys.argv[2]), int(sys.argv[3]), tuple(int(x) for x in sys.argv[4].split(","))
        print(f"p={p} e={e} pattern {pat}: {full_check(p, e, pat)}", flush=True)
        return
    e = int(sys.argv[1])
    pmax = int(sys.argv[2]) if len(sys.argv) > 2 else threshold(e)
    print(f"# e = {e}: Weil threshold {threshold(e)}, scanning primes p = 1 mod {e} up to {pmax}", flush=True)
    for p in range(e + 1, pmax + 1):
        if (p - 1) % e or any(p % q == 0 for q in range(2, int(p ** 0.5) + 1)):
            continue
        t = time.time()
        free = k4_free_patterns_fast(p, e)
        print(f"p={p} h={(p - 1) // e}: {len(free)} K4-free single-class patterns (S_01 = C_0) {free} "
              f"({time.time() - t:.1f} s)", flush=True)


if __name__ == "__main__":
    main()

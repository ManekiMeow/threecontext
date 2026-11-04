"""Cyclotomic four-partite graphs over Z_p (RESULTS.md §3.12, family B).

Parts are copies of Z_p; (k, g) ~ (l, h) iff h - g ∈ S_kl, where every S_kl is a union of
cyclotomic classes of index e (cosets of the subgroup of e-th powers in Z_p^*), possibly with 0.
These graphs are Cayley graphs over Z_p (family A of cayley4.py) with the extra symmetry
x -> m x, m an e-th power. They are enumerated exhaustively, up to multiplying by a generator.
There is no normalisation 0 ∈ S_k3: a translation of one basis does not commute with x -> m x.
Every member is tested against the necessary conditions of RESULTS.md §2–3 ((K4), (DEG), (C2),
(C4), (TRI)), then the local GHZ test theta = 3, the 0-1 colouring of the complement, and
realizability (cayley4.realize).

usage: python3 cyclo4.py PMIN PMAX [EMAX]      (primes p in [PMIN, PMAX], indices e <= EMAX)
"""
import itertools
import sys
import time

import numpy as np

import cayley4 as c4

PAIRS = c4.PAIRS


def primitive_root(p):
    fac = [q for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, q))]
    return next(g for g in range(2, p) if all(pow(g, (p - 1) // q, p) != 1 for q in fac))


def classes(p, e):
    g = primitive_root(p)
    return [sorted({pow(g, i + e * j, p) for j in range((p - 1) // e)}) for i in range(e)]


def conv(a, b):
    """(a * b)(s) = #{c : a(c), b(s - c)} over Z_p, for 0/1 vectors."""
    return np.rint(np.real(np.fft.ifft(np.fft.fft(a) * np.fft.fft(b)))).astype(int)


def corr(a, b):
    """#{c : a(c), b(c - s)} as a function of s."""
    return np.rint(np.real(np.fft.ifft(np.fft.fft(a) * np.conj(np.fft.fft(b))))).astype(int)


def check(p, S):
    """Necessary conditions for the vectors S[(k, l)] (k < l). Returns the first failure or None."""
    ind = lambda k, l: S[(k, l)] if k < l else np.roll(S[(l, k)][::-1], 1)      # S_lk = -S_kl
    sizes = {pr: int(S[pr].sum()) for pr in PAIRS}
    for k in range(4):
        ss = [sizes[tuple(sorted((k, l)))] for l in range(4) if l != k]
        if min(ss) < 9 or sum(ss) < 30:
            return "DEG"
        for i, a in enumerate(ss):
            if a == 9 and min(ss[j] for j in range(3) if j != i) < 12:
                return "DEG(Thm5)"
    # (K4): the neighbourhood of (3, 0) is triangle-free
    N = [np.where(ind(3, k))[0] for k in range(3)]
    M01 = ind(0, 1)[(N[1][None, :] - N[0][:, None]) % p]
    M02 = ind(0, 2)[(N[2][None, :] - N[0][:, None]) % p]
    M12 = ind(1, 2)[(N[2][None, :] - N[1][:, None]) % p]
    if ((M01.astype(int) @ M12.astype(int)) * M02).sum() > 0:
        return "K4"
    for k, l in PAIRS:                                                # (C2), (C4)
        for m in range(4):
            if m in (k, l):
                continue
            cnt = conv(ind(k, m).astype(int), ind(m, l).astype(int))
            if (cnt[ind(k, l)] == 0).any():
                return "C2"
            if (cnt[~ind(k, l)] == 1).any():
                return "C4"
    for k in range(4):                                                # (C4), same basis
        for m in range(4):
            if m != k:
                a = ind(k, m).astype(int)
                if (corr(a, a)[1:] == 1).any():
                    return "C4same"
    for k, l, m in itertools.combinations(range(4), 3):               # (TRI)
        A = ind(l, m)[(np.arange(p)[None, :] - np.arange(p)[:, None]) % p]   # [s, c]: c - s ∈ S_lm
        if not (A[np.ix_(np.where(ind(k, l))[0], np.where(ind(k, m))[0])]).any():
            return "TRI"
    return None


def to_sets(p, S):
    """S[(k, l)] for all ordered pairs, as index lists (the format of Cayley4.sets)."""
    out = {}
    for k in range(4):
        for l in range(4):
            if k < l:
                out[(k, l)] = list(np.where(S[(k, l)])[0])
            elif k > l:
                out[(k, l)] = list(np.where(np.roll(S[(l, k)][::-1], 1))[0])
    return out


class Zp:
    """Z_p as a cayley4 group (identity 0)."""

    def __init__(self, p):
        self.n, self.name = p, f"Z{p}"
        self.mul = (np.arange(p)[:, None] + np.arange(p)[None, :]) % p
        self.inv = (-np.arange(p)) % p


def scan(p, e, log=print, realize_restarts=3):
    C = classes(p, e)
    blocks = [np.isin(np.arange(p), c) for c in C] + [np.arange(p) == 0]   # class 0..e-1, then {0}
    opts = []
    for mask in range(1, 1 << (e + 1)):
        opts.append((mask, np.any([blocks[i] for i in range(e + 1) if mask >> i & 1], axis=0)))
    stats = dict(p=p, e=e, total=0)
    enc = c4.Cayley4.__new__(c4.Cayley4)                 # only the helper methods are used
    enc.G, enc.d = Zp(p), p
    survivors = []
    rot = lambda mask: ((mask << 1) & ((1 << e) - 1)) | ((mask >> (e - 1)) & 1) | (mask & (1 << e))
    for combo in itertools.product(range(len(opts)), repeat=6):
        masks = [opts[i][0] for i in combo]
        # canonical up to multiplying everything by the generator (rotates the classes)
        if any(tuple(r) < tuple(masks) for r in _rotations(masks, rot, e)):
            continue
        stats["total"] += 1
        S = {pr: opts[i][1] for pr, i in zip(PAIRS, combo)}
        why = check(p, S)
        stats[why] = stats.get(why, 0) + 1
        if why is None:
            survivors.append(S)
    stats["survivors"] = len(survivors)
    results = []
    for S in survivors:
        sets = to_sets(p, S)
        thetas = []
        for k in (3, 0, 1, 2):
            V, A = enc.neighbourhood(sets, k)
            thetas.append(2.0 if c4.two_colouring(A) is not None else c4.theta_complement(V, A))
        res = dict(sizes={f"{k}{l}": int(S[(k, l)].sum()) for k, l in PAIRS}, thetas=[round(t, 4) for t in thetas])
        if min(thetas) > 3 - 1e-6:
            X = enc.full_graph(sets)
            res["complement_colourable"] = c4.colourable_complement(X, p)
            res["realize_residual"] = float(c4.realize(X, p, restarts=realize_restarts)[0])
        results.append(res)
    stats["theta_pass"] = sum(min(r["thetas"]) > 3 - 1e-6 for r in results)
    return stats, results


def _rotations(masks, rot, e):
    cur = list(masks)
    for _ in range(e - 1):
        cur = [rot(m) for m in cur]
        yield cur


if __name__ == "__main__":
    pmin, pmax = int(sys.argv[1]), int(sys.argv[2])
    emax = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    for p in range(pmin, pmax + 1):
        if p < 2 or any(p % q == 0 for q in range(2, int(p ** 0.5) + 1)):
            continue
        for e in range(2, emax + 1):
            if (p - 1) % e:
                continue
            t = time.time()
            stats, results = scan(p, e)
            print(f"p={p} e={e}: {stats} ({time.time() - t:.0f} s)", flush=True)
            for r in results:
                print("   ", r, flush=True)

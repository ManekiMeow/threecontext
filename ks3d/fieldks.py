"""KS sets among rays with coordinates a + b sqrt(d), |a| <= h, |b| <= hb.

usage: python3 fieldks.py d h hb T   -> decide whether a KS set of size <= T exists
"""
import random
import sys
from itertools import permutations, product

from cegar import decide, ray_perms
from minks import Instance, log, sym_actions
from qrays import field_rays, prune, structure


def field_signed_perms(K):
    out = []
    for p in permutations(range(3)):
        for s in product((1, -1), repeat=3):
            out.append(lambda v, p=p, s=s: K.canon(tuple((s[i] * v[p[i]][0], s[i] * v[p[i]][1])
                                                          for i in range(3))))
    return out


def show(K, r):
    def f(x):
        a, b = x
        if b == 0:
            return str(a)
        rb = f"{'' if abs(b) == 1 else abs(b)}r{K.d}"
        if a == 0:
            return ("-" if b < 0 else "") + rb
        return f"{a}{'+' if b > 0 else '-'}{rb}"
    return "(" + ", ".join(f(x) for x in r) + ")"


if __name__ == "__main__":
    d, h, hb, T = map(int, sys.argv[1:5])
    random.seed(0)
    K, R = field_rays(d, h, hb)
    R = prune(K, R)
    e, t = structure(K, R)
    inst = Instance(R, e, t)
    G = field_signed_perms(K)
    log(f"Q(sqrt{d}) h={h} hb={hb}: {inst.n} rays, {len(e)} pairs, {len(t)} bases; target <= {T}")
    S = decide(inst, T, sym=sym_actions(inst, G), perms=ray_perms(inst, G))
    if S:
        log("KS set:", " ".join(show(K, R[i]) for i in S))

"""Small finite groups as multiplication tables (used by ks/cayley4.py).

A Group has elements 0..n-1 with 0 the identity, mul[a][b] = ab and inv[a] = a^{-1}.
The list below is not a classification: it is every abelian group, the dihedral, dicyclic and
metacyclic groups Z_m ⋊ Z_k, and a few direct products, for each order in a range.
"""
import itertools
from math import gcd

import numpy as np


class Group:
    def __init__(self, name, elems, op):
        elems = list(elems)
        idx = {x: i for i, x in enumerate(elems)}
        n = len(elems)
        table = [[idx[op(a, b)] for b in elems] for a in elems]
        e = next(i for i in range(n) if all(table[i][j] == j for j in range(n)))
        order = [e] + [i for i in range(n) if i != e]                 # identity first
        pos = {old: new for new, old in enumerate(order)}
        self.name, self.n = name, n
        self.mul = np.array([[pos[table[a][b]] for b in order] for a in order], dtype=np.int32)
        self.inv = np.array([int(np.where(self.mul[a] == 0)[0][0]) for a in range(n)], dtype=np.int32)
        self.abelian = bool((self.mul == self.mul.T).all())

    def invariant(self):
        """Sorted element orders (cheap isomorphism invariant, used to drop duplicates)."""
        orders = []
        for a in range(self.n):
            k, x = 1, a
            while x != 0:
                x, k = self.mul[x][a], k + 1
            orders.append(k)
        return (self.n, self.abelian, tuple(sorted(orders)))


def cyclic(n):
    return Group(f"Z{n}", range(n), lambda a, b: (a + b) % n)


def invariant_factors(ns):
    """Prime-power cyclic orders -> invariant factors (Z2 x Z9 -> Z18, Z2 x Z2 x Z3 -> Z2 x Z6)."""
    by_p = {}
    for q in ns:
        p = next(p for p in range(2, q + 1) if q % p == 0)
        by_p.setdefault(p, []).append(q)
    cols = max(len(v) for v in by_p.values())
    out = [1] * cols
    for v in by_p.values():
        for i, q in enumerate(sorted(v, reverse=True)):
            out[cols - 1 - i] *= q
    return out


def abelian(ns):
    name = "x".join(f"Z{m}" for m in invariant_factors(ns))
    return Group(name, itertools.product(*[range(m) for m in ns]),
                 lambda a, b: tuple((x + y) % m for x, y, m in zip(a, b, ns)))


def gen_dihedral(ns):
    """A ⋊ Z2 with Z2 acting by inversion on the abelian group A = prod Z_m."""
    name = "Dih(" + "x".join(f"Z{m}" for m in invariant_factors(ns)) + ")"
    def op(a, b):
        x, s = a
        y, t = b
        y = y if s == 0 else tuple((-c) % m for c, m in zip(y, ns))
        return (tuple((u + v) % m for u, v, m in zip(x, y, ns)), s ^ t)
    return Group(name, itertools.product(itertools.product(*[range(m) for m in ns]), range(2)), op)


def dihedral(n):
    """Order 2n: (i, s) = r^i s^s."""
    return Group(f"D{2 * n}", itertools.product(range(n), range(2)),
                 lambda a, b: ((a[0] + (b[0] if a[1] == 0 else -b[0])) % n, a[1] ^ b[1]))


def dicyclic(n):
    """Order 4n: <a, x | a^{2n} = 1, x^2 = a^n, x a x^-1 = a^-1>; (i, j) = a^i x^j."""
    def op(p, q):
        (i, j), (k, l) = p, q
        if j == 0:
            return ((i + k) % (2 * n), l)
        if l == 0:
            return ((i - k) % (2 * n), 1)
        return ((i - k + n) % (2 * n), 0)
    return Group(f"Dic{4 * n}", itertools.product(range(2 * n), range(2)), op)


def metacyclic(m, k, r):
    """Z_m ⋊ Z_k, the generator of Z_k acting by multiplication by r (r^k = 1 mod m)."""
    return Group(f"Z{m}:Z{k}[{r}]", itertools.product(range(m), range(k)),
                 lambda a, b: ((a[0] + pow(r, a[1], m) * b[0]) % m, (a[1] + b[1]) % k))


def direct(G, H):
    els = list(itertools.product(range(G.n), range(H.n)))
    return Group(f"{G.name}x{H.name}", els,
                 lambda a, b: (int(G.mul[a[0]][b[0]]), int(H.mul[a[1]][b[1]])))


def symmetric4():
    return Group("S4", itertools.permutations(range(4)), lambda a, b: tuple(a[b[i]] for i in range(4)))


def alternating4():
    even = [p for p in itertools.permutations(range(4))
            if sum(p[i] > p[j] for i in range(4) for j in range(i + 1, 4)) % 2 == 0]
    return Group("A4", even, lambda a, b: tuple(a[b[i]] for i in range(4)))


def sl23():
    mats = [m for m in itertools.product(range(3), repeat=4) if (m[0] * m[3] - m[1] * m[2]) % 3 == 1]
    return Group("SL(2,3)", mats, lambda a, b: ((a[0] * b[0] + a[1] * b[2]) % 3, (a[0] * b[1] + a[1] * b[3]) % 3,
                                                (a[2] * b[0] + a[3] * b[2]) % 3, (a[2] * b[1] + a[3] * b[3]) % 3))


def abelian_types(n):
    """All abelian groups of order n, as lists of cyclic orders (prime-power factors)."""
    fac, m, p = {}, n, 2
    while m > 1:
        while m % p == 0:
            fac[p] = fac.get(p, 0) + 1
            m //= p
        p += 1

    def partitions(k, top=None):
        top = k if top is None else top
        if k == 0:
            yield []
        for i in range(min(k, top), 0, -1):
            for rest in partitions(k - i, i):
                yield [i] + rest
    per = [[[p ** e for e in part] for part in partitions(k)] for p, k in fac.items()]
    for choice in itertools.product(*per):
        yield sorted(x for part in choice for x in part)


def groups_of_order(n):
    out = []
    for t in abelian_types(n):
        out.append(cyclic(n) if len(invariant_factors(t)) == 1 else abelian(t))
    if n % 2 == 0 and n >= 6:
        out.append(dihedral(n // 2))
        for t in abelian_types(n // 2):
            if len(invariant_factors(t)) > 1:
                out.append(gen_dihedral(t))
    if n % 4 == 0 and n >= 8:
        out.append(dicyclic(n // 4))
    for k in range(2, n):
        if n % k:
            continue
        m = n // k
        for r in range(2, m):
            if gcd(r, m) == 1 and pow(r, k, m) == 1:
                out.append(metacyclic(m, k, r))
    small = {6: [dihedral(3)], 8: [dihedral(4), dicyclic(2)], 12: [alternating4(), dicyclic(3), dihedral(6)],
             24: [symmetric4(), sl23(), direct(alternating4(), cyclic(2))]}
    for a in (6, 8, 12):
        if n % a == 0 and n // a > 1:
            for H in small[a]:
                for t in abelian_types(n // a):
                    out.append(direct(H, cyclic(n // a) if len(invariant_factors(t)) == 1 else abelian(t)))
    out += small.get(n, []) if n == 24 else []
    seen, uniq = set(), []
    for G in out:
        key = G.invariant()
        if key not in seen:
            seen.add(key)
            uniq.append(G)
    return uniq


if __name__ == "__main__":
    for n in range(18, 41):
        gs = groups_of_order(n)
        print(n, len(gs), " ".join(G.name for G in gs))

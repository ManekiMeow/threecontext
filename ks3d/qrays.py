"""Rays in R^3 with coordinates in Q(sqrt d), exact arithmetic.

A field element a + b*sqrt(d) is a pair (a, b) of Fractions. A ray is normalised so that its
first nonzero coordinate equals 1.
"""
from fractions import Fraction as F
from itertools import product


class QF:
    def __init__(self, d):
        self.d = d

    def mul(self, x, y):
        return (x[0] * y[0] + self.d * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    def inv(self, x):
        n = x[0] * x[0] - self.d * x[1] * x[1]
        return (x[0] / n, -x[1] / n)

    def add(self, x, y):
        return (x[0] + y[0], x[1] + y[1])

    def sub(self, x, y):
        return (x[0] - y[0], x[1] - y[1])

    def canon(self, v):
        for x in v:
            if x != (0, 0):
                ix = self.inv(x)
                return tuple(self.mul(y, ix) for y in v)
        raise ValueError("zero")

    def dot(self, u, v):
        s = (F(0), F(0))
        for a, b in zip(u, v):
            s = self.add(s, self.mul(a, b))
        return s

    def cross(self, u, v):
        m, s = self.mul, self.sub
        return (s(m(u[1], v[2]), m(u[2], v[1])), s(m(u[2], v[0]), m(u[0], v[2])),
                s(m(u[0], v[1]), m(u[1], v[0])))


def field_rays(d, h, hb=None):
    """Rays spanned by vectors with coordinates a + b sqrt(d), |a| <= h, |b| <= hb (default h)."""
    hb = h if hb is None else hb
    K = QF(d)
    vals = [(F(a), F(b)) for a in range(-h, h + 1) for b in range(-hb, hb + 1)]
    out = set()
    for v in product(vals, repeat=3):
        if any(x != (0, 0) for x in v):
            out.add(K.canon(v))
    return K, sorted(out)


def structure(K, rays):
    idx = {r: i for i, r in enumerate(rays)}
    n = len(rays)
    zero = (0, 0)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if K.dot(rays[i], rays[j]) == zero]
    tris = set()
    for i, j in edges:
        k = idx.get(K.canon(K.cross(rays[i], rays[j])))
        if k is not None:
            tris.add(tuple(sorted((i, j, k))))
    return edges, sorted(tris)


def prune(K, rays):
    while True:
        _, tris = structure(K, rays)
        keep = sorted({i for t in tris for i in t})
        if len(keep) == len(rays):
            return rays
        rays = [rays[i] for i in keep]

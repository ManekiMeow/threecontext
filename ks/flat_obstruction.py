"""Counting obstruction to flat realizations (RESULTS.md §3.13, Proposition 27).

Flat: every nonzero overlap between different bases has modulus 1/sqrt(h). Gauge the rays of N(u) so that
<u|z> > 0. If every non-orthogonal pair has t common neighbours per basis, the t terms of <u|y> = sum_c <u|c><c|y>
force every edge inside N(u) to have phase within theta_max = arccos(sqrt(h) - t + 1). For orthogonal x, x' in N(u)
and a third basis b, the common neighbours of x and x' in b that lie in N(u) give terms with real part
>= cos(2 theta_max)/h. A flat realization therefore needs  inside * cos(2 theta_max) <= outside  for every such
(u, x, x', b). A negative slack proves that no flat realization exists.

usage: python3 flat_obstruction.py ks/cayley4_candidate_pP.json ...
"""
import itertools
import sys

import numpy as np

from realize_search import load_pattern


def t_range(X):
    """Smallest and largest number of common neighbours, in one other basis, of a non-orthogonal pair."""
    n = len(X)
    d = n // 4
    part = np.arange(n) // d
    ts = set()
    for k in range(4):
        u = k * d
        for y in np.where(X[u])[0]:
            for b in range(4):
                if b not in (k, part[y]):
                    ts.add(int((X[u] & X[y] & (part == b)).sum()))
    return min(ts), max(ts)


def obstruction(X):
    n = len(X)
    d = n // 4
    part = np.arange(n) // d
    h = int(X[0, d:2 * d].sum())
    tmin, tmax = t_range(X)
    if tmin < np.sqrt(h):                                # a single pair already violates t >= sqrt(h)
        return h, (tmin, tmax), 0.0, None, (-np.inf, "t_min < sqrt(h)"), 1
    t = tmin                                             # weakest phase bound over all pairs (conservative)
    c = np.clip(np.sqrt(h) - (t - 1), -1, 1)
    cos2 = np.cos(2 * np.arccos(c))
    worst, count = (np.inf, None), 0
    for k in range(4):
        u = k * d                                        # one ray per basis suffices (Z_p acts transitively)
        N = np.where(X[u])[0]
        inN = np.zeros(n, bool)
        inN[N] = True
        for x, xp in itertools.combinations(N, 2):
            if X[x, xp]:
                continue
            for b in range(4):
                if b in (k, part[x], part[xp]):
                    continue
                com = X[x] & X[xp] & (part == b)
                ins, out = int((com & inN).sum()), int((com & ~inN).sum())
                slack = out - ins * cos2
                count += slack < 0
                if slack < worst[0]:
                    worst = (slack, (k, int(x), int(xp), b, ins, out))
    return h, (tmin, tmax), np.degrees(np.arccos(c)), cos2, worst, count


if __name__ == "__main__":
    for f in sys.argv[1:]:
        h, ts, th, cos2, worst, count = obstruction(load_pattern(f))
        if cos2 is None:
            print(f"{f}: h = {h}, t in {ts}, sqrt(h) = {np.sqrt(h):.2f}: some pair has t < sqrt(h) "
                  f"-> NO FLAT REALIZATION", flush=True)
            continue
        print(f"{f}: h = {h}, t in {ts}, theta_max = {th:.1f} deg, cos(2 theta_max) = {cos2:.3f}; "
              f"min slack {worst[0]:.3f} at (basis of u, x, x', b, inside, outside) = {worst[1]}; "
              f"{count} contradictions -> {'NO FLAT REALIZATION' if count else 'no contradiction'}", flush=True)

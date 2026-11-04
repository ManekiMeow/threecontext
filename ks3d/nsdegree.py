"""Degree of algebraic refutations of 010-colourability (RESULTS.md §7.4).

010-colourings of a set with orthogonality graph G and bases T are the 0/1 solutions of
    x_i^2 = x_i,   x_i x_j = 0 (ij in E),   g_T := sum_{t in T} x_t - 1 = 0 (T in T).
Work in the quotient ring R = Q[x]/(x_i^2 - x_i, x_i x_j : ij in E), whose monomials are x_S for
independent sets S. A Nullstellensatz refutation of degree d is an identity
    1 = sum_{T, |S| <= d-1} c_{S,T} x_S g_T   in R,
i.e. 1 lies in the span of the reduced rows x_S g_T. Every quantum state gives a functional L with
L(1) = 1 killing all rows of degree <= 2 (L(x_i) = <Pi_i>, L(x_i x_j) = Re<Pi_i Pi_j>), so d >= 3
for every set of rays in C^3, and the maximally mixed state gives one killing all rows of degree
<= 3 (Theorem 7.3). This script checks that functional and finds the smallest refutation degree,
exactly over Q (Fractions) and over GF(2), GF(3) and two 31-bit primes.
"""
import itertools
import sys

import numpy as np

import spectral

PRIMES = (2_147_483_629, 2_147_483_587)


def independent_sets(n, A, dmax):
    """All independent sets of size <= dmax, as sorted tuples."""
    out = [()]
    level = [()]
    for _ in range(dmax):
        nxt = []
        for S in level:
            start = S[-1] + 1 if S else 0
            for v in range(start, n):
                if all(not A[v, u] for u in S):
                    nxt.append(S + (v,))
        out += nxt
        level = nxt
    return out


def rows(n, A, T, d):
    """Reduced rows x_S g_T for |S| <= d-1, as dicts {monomial: coeff}."""
    for S in independent_sets(n, A, d - 1):
        for t3 in T:
            row = {}
            for t in t3:
                if t in S:
                    U = S
                elif any(A[t, u] for u in S):
                    continue                      # x_S x_t = 0
                else:
                    U = tuple(sorted(S + (t,)))
                row[U] = row.get(U, 0) + 1
            row[S] = row.get(S, 0) - 1
            row = {k: c for k, c in row.items() if c}
            if row:
                yield row


def refutes(n, A, T, d, p):
    """Is 1 (the empty monomial) in the span of the degree-d rows, modulo p?
    Gaussian elimination with pivots on the largest monomial in (degree, lex) order."""
    key = lambda U: (len(U), U)
    piv = {}                                      # leading monomial -> reduced row
    for row in rows(n, A, T, d):
        row = {k: c % p for k, c in row.items() if c % p}
        while row:
            lead = max(row, key=key)
            if lead not in piv:
                inv = pow(row[lead], p - 2, p)
                piv[lead] = {k: c * inv % p for k, c in row.items()}
                break
            c = row[lead]
            for k, v in piv[lead].items():
                row[k] = (row.get(k, 0) - c * v) % p
                if row[k] == 0:
                    del row[k]
    # 1 is in the span iff the empty monomial is a pivot of a row with no other terms after
    # full reduction; since () is the smallest monomial, check whether () is a leading monomial.
    return () in piv


def mixed_state_functional(g):
    """L(x_S) = Re tr(prod_{s in S} Pi_s) / 3 for |S| <= 3 (the maximally mixed state).
    Real parts of traces of products of <= 3 projectors are symmetric in the factors, and
    sum_{t in T} Pi_t = I, so L kills every row of degree <= 3 (Theorem 7.3)."""
    P = [np.outer(v, v.conj()) for v in g.rays]

    def L(S):
        M = np.eye(3)
        for s in S:
            M = M @ P[s]
        return np.trace(M).real / 3
    return L


def check_functional(g, L, d):
    """Largest |L(row)| over the rows of degree <= d, and L(1)."""
    worst = 0.0
    for row in rows(g.n, g.A, g.T, d):
        worst = max(worst, abs(sum(c * L(S) for S, c in row.items())))
    return worst


def refutes_exact(n, A, T, d):
    """Same elimination as refutes() over Q (Fractions): an exact certificate check."""
    from fractions import Fraction
    key = lambda U: (len(U), U)
    piv = {}
    for row in rows(n, A, T, d):
        row = {k: Fraction(c) for k, c in row.items()}
        while row:
            lead = max(row, key=key)
            if lead not in piv:
                c0 = row[lead]
                piv[lead] = {k: c / c0 for k, c in row.items()}
                break
            c = row[lead]
            for k, v in piv[lead].items():
                row[k] = row.get(k, 0) - c * v
                if row[k] == 0:
                    del row[k]
    return () in piv


if __name__ == "__main__":
    sets = spectral.known_sets()
    names = sys.argv[1:] or ["KS 31 (height 2)", "KS 33 (Peres)", "Yu-Oh closure 25"]
    for name in names:
        g = spectral.OG(sets[name])
        L = mixed_state_functional(g)
        print(f"{name}: maximally mixed functional, max |L(row)| over rows of degree <= 3: "
              f"{check_functional(g, L, 3):.1e}, degree 4: {check_functional(g, L, 4):.3f}", flush=True)
        print(f"   over Q (exact): refutation of degree 3: {refutes_exact(g.n, g.A, g.T, 3)}, "
              f"degree 4: {refutes_exact(g.n, g.A, g.T, 4)}", flush=True)
        for p in (2, 3) + PRIMES:
            d = next((d for d in range(1, 5) if refutes(g.n, g.A, g.T, d, p)), None)
            print(f"   over GF({p}): smallest refutation degree <= 4: {d}", flush=True)

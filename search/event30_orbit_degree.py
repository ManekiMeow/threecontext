"""Exact degree-pair bounds for automorphism orbits in a 30-event support.

Let an automorphism group preserve each size-ten context and have an orbit of
size m in one context. Its vertices have a common pair of cross-degrees q,r.
Every event has c<1/9. Hence alpha>(q-1)/9 and beta>(r-1)/9, and Lemma14 gives

    c <= 1 / (1/alpha + 1/beta + 1/(1-alpha-beta)).

Maximize this expression over the relaxed closed lower bounds alpha>=a,
beta>=b. Equivalently minimize the strictly convex reciprocal sum. With a>=b,
the optimum is (1/3,1/3,1/3) if a<=1/3, (a,(1-a)/2,(1-a)/2) if b<=(1-a)/2,
and (a,b,1-a-b) otherwise. All numbers used here are exact fractions.

For m<10 the orbit mass is >(m-1)/9, from the strict probability bound on its
complement; for m10 it is exactly1. If m times the relaxed event upper bound
is at most (m-1)/9 for m<10, or less than1 for m10, the orbit is impossible.
No unitary realization of the automorphism or Gram averaging is assumed.

The resulting sorted degree pairs are componentwise bounded by:
 m1: q+r<=10;
 m2: (3,7), (4,6), or (5,5);
 m3,4,5: (3,6) or (4,5);
 m6,7,8,9,10: (4,5).
Passing this necessary condition never establishes a quantum realization.
"""
from fractions import Fraction as Q
import json


def relaxed_event_upper(q, r):
    """Supremum under degree-implied neighbor-mass lower bounds."""
    assert 1 <= q <= 9 and 1 <= r <= 9
    a, b = sorted((Q(q - 1, 9), Q(r - 1, 9)), reverse=True)
    if a + b >= 1:
        return Q(0)
    if a <= Q(1, 3):
        return Q(1, 9)
    if b <= (1 - a) / 2:
        return a * (1 - a) / (1 + 3 * a)
    return 1 / (1 / a + 1 / b + 1 / (1 - a - b))


def possible_degree_pair(m, q, r):
    assert 1 <= m <= 10
    bound = m * relaxed_event_upper(q, r)
    return bound >= 1 if m == 10 else bound > Q(m - 1, 9)


def maximal_degree_pairs(m):
    pairs = [(q, r) for q in range(1, 10) for r in range(q, 10)
             if possible_degree_pair(m, q, r)]
    return [(q, r) for q, r in pairs
            if not any(s >= q and t >= r and (s, t) != (q, r) for s, t in pairs)]


def main():
    expected = {
        1: [(1, 9), (2, 8), (3, 7), (4, 6), (5, 5)],
        2: [(3, 7), (4, 6), (5, 5)],
        3: [(3, 6), (4, 5)],
        4: [(3, 6), (4, 5)],
        5: [(3, 6), (4, 5)],
        **{m: [(4, 5)] for m in range(6, 11)},
    }
    results = {m: maximal_degree_pairs(m) for m in range(1, 11)}
    assert results == expected
    assert relaxed_event_upper(1, 8) == Q(7, 135)
    assert relaxed_event_upper(1, 7) == Q(2, 27)
    assert relaxed_event_upper(4, 6) == Q(5, 69)
    assert relaxed_event_upper(5, 5) == Q(2, 27)
    assert relaxed_event_upper(1, 6) == Q(5, 54)
    print(json.dumps({"maximal_sorted_degree_pairs": results,
                      "scope": "necessary conditions only; all arithmetic exact"}, indent=2))


if __name__ == "__main__":
    main()

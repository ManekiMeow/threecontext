# Exact arithmetic behind RESULTS.md 2, Corollary 5 (n >= 30, d >= 18).
# A realization in dimension d with n events has F-spectrum {3} + (d-1) eigenvalues f in (0,3) (3 is simple for
# d <= 33) with sum f = n - 3 and sum g(f) = -6, g(f) = f(f-1)(f-2).  The convex envelope of g on [0,3] is
# -f/4 on [0,3/2] and g on [3/2,3], so (d-1) * env((n-3)/(d-1)) <= -6 is necessary; print the largest such n.
from fractions import Fraction as Fr
env = lambda m: -m / 4 if m <= Fr(3, 2) else m * (m - 1) * (m - 2)
print(" d | n allowed by the spectrum   (Theorem 1: n >= 27; Corollary 5: n >= 30, and n >= 33 with a 9-event context)")
for d in range(17, 23):
    ok = [n for n in range(27, 80) if (d - 1) * env(Fr(n - 3, d - 1)) <= -6]
    assert ok == list(range(27, max(ok) + 1))
    print(f"{d:2d} | 27 <= n <= {max(ok)}")
# d = 17 needs n <= 29 < 30: so d >= 18.  A 9-event context forces d >= 20 (Theorem 4), so d = 18, 19 need all
# contexts of size >= 10 and n <= 32 resp. n <= 34.

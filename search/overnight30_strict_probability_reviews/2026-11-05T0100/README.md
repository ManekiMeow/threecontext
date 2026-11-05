# Strict probability LP trial on the 218 residual C5 quotients

This bounded trial adds strict orbit-probability bounds and the strict
weighted-degree inequalities of Theorem 56 to the existing sector and
disjoint-neighborhood moment constraints. It uses the new independent checker
`search/event30_strict_probability_lp.py`; the older certificate formats and
immutable catalogs are unchanged.

The checker accepts a zero right side only when the exact nonnegative Farkas
combination includes a positive multiplier on a mathematically strict
inequality. The total multiplier on strict rows is normalized to one. A
nonpositive right side then contradicts the strictly smaller left side.
Floating-point LP output is never accepted without this rational identity.

All **218 residual published quotient classes** were attempted in 5.33 seconds
of a 25-second budget, after the ordinary 01:00 review completed the fixed
published pool. **No new exact certificate was produced.** This is an
inconclusive search result, not a proof of feasibility of either the scalar
relaxation or a quantum realization. Subsequent phase enumeration can and does
exclude classes that these scalar constraints do not certify.

Two exact degree-two/degree-five zero-margin controls are accepted, while six
false certificates (including a zero identity with no strict multiplier) are
rejected. The helper explicitly refuses optimized Python, which disables
assertions. Run:

```sh
python search/event30_strict_probability_lp.py
python search/overnight30_strict_probability_reviews/2026-11-05T0100/verify.py
```

The trial metadata is bound to the two previous residual catalogs. There are
no quantum feasibility claims and no solver-derived exclusions in this bundle.

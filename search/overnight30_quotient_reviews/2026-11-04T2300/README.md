# Exact C5 quotient review

Published experiment commit `f37683f803842de5b3a8d7fcffcbdd781bcff492`
contains 25,675 scored C5-double records with 5,273 different quotient classes
under context permutations and within-context orbit swaps. The bounded review
selects the 150 classes with highest numerical representative scores.

The harmonic probability inequalities give exact rational Farkas certificates
for 126 selected classes. Of these, 108 passed the fast screen preceding
probability propagation. Every certificate excludes every support with its
actual automorphism-orbit quotient, regardless of cyclic edge phases. The
other classes are unresolved by this test; a feasible probability relaxation
does not establish a quantum realization.

Run using standard Python:

```
python search/overnight30_quotient_reviews/2026-11-04T2300/verify.py
```

This checks all 126 rational contradictions, the 150 representative graphs'
triangle-freeness and actual simultaneous C5 automorphisms, quotient counts,
and distinctness under all 48 permitted orbit relabelings. It does not trust
the numerical scores or an equitable partition alone. The certificate
framework received independent mathematical review, a feasible regular(3,3)
probability control, and fourteen rejected tamper cases.

`discover.py --data-repo /workspace/ghz30-overnight` reproduces collection and
bounded certificate discovery from the pinned published input using the
research Python environment. It overwrites this review, so verification is
the appropriate command for preserving the existing evidence. Discovery uses
SciPy and SymPy; proof verification has no third-party dependency.

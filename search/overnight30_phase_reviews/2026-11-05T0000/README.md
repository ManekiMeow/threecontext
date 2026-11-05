# Complete singleton-tree C5 phase exclusions

These five fixed quotient classes admit no thirty-event GHZ realization with
the simultaneous C5 action on six consecutive five-vertex orbits. The proof
covers every cyclic phase choice; it assumes no symmetry or uniformity of the
realizing vectors. Quotient indices refer to the immutable, hash-bound
`search/overnight30_quotient_reviews/2026-11-04T2300/summary.json`.

| Index | Normalized triangle-free assignments | Unit representatives | Exact uniform theta upper bound |
|---|---:|---:|---|
| 1 | 330 | 83 | 2991515946933/1000000000000 |
| 14 | 2448 | 612 | 579998681849/200000000000 |
| 75 | 408 | 102 | 177916739719/62500000000 |
| 107 | 96 | 24 | 2837622201683/1000000000000 |
| 124 | 816 | 204 | 115348376959/40000000000 |

All **1025** representatives have graph-bound exact integer factor certificates.
The largest bound is **2.991515946933 < 3**. Singleton tree edges fix the five
relative fiber translations; common nonzero multipliers identify further
phases. These are coverage representatives, not a graph-isomorphism census.

From the repository root, verification needs only standard Python:

```sh
for index in 001 014 075 107 124; do
  python3 search/event30_c5_phases.py verify search/overnight30_phase_reviews/2026-11-05T0000/index_$index
done
python3 search/overnight30_phase_reviews/2026-11-05T0000/independent_coverage.py
```

The independent audit uses fixed-order explicit residue sets, checks all 3125
relative translations for one representative of each class, and verifies every
exact dual against independently generated graph edges. Its saved report binds
the manifest hashes. The worker loads a class only after complete verification.

Together with the separate generalized-tree directory, these are the nine
finite phase exclusions of Theorem 54. They do not exhaust all C5 quotients or
all (10,10,10) supports.

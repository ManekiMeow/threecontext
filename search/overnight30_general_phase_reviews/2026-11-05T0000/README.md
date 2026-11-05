# Complete proper-subset-tree C5 phase exclusions

These four fixed quotient classes admit no thirty-event GHZ realization with
the simultaneous C5 action on six consecutive five-vertex orbits. All cyclic
phases are covered without assuming symmetry or uniformity of realizing vectors.
The indices refer to the immutable, hash-bound
`search/overnight30_quotient_reviews/2026-11-04T2300/summary.json`.

| Index | Normalized triangle-free assignments | Unit representatives | Exact uniform theta upper bound |
|---|---:|---:|---|
| 3 | 20 | 6 | 579826603527/200000000000 |
| 11 | 60 | 16 | 2877076947253/1000000000000 |
| 103 | 960 | 240 | 2911162685863/1000000000000 |
| 147 | 162 | 41 | 2833038240961/1000000000000 |

All **303** representatives have exact integer factor certificates; the largest
uniform bound is **2.911162685863 < 3**. Any nonempty proper subset of Z5 has
trivial translation stabilizer. Choosing its least translate along a spanning
tree therefore fixes the five relative fiber translations uniquely. Common
units act followed by renormalization. Counts are coverage representatives,
not a graph-isomorphism census.

From the repository root, verification needs only standard Python:

```sh
for index in 003 011 103 147; do
  python3 search/event30_c5_general_phases.py verify search/overnight30_general_phase_reviews/2026-11-05T0000/index_$index
done
python3 search/overnight30_general_phase_reviews/2026-11-05T0000/independent_coverage.py
```

The independent audit uses fixed-order explicit residue sets and verifies all
303 graph-bound dual proofs, 19,232 unit-composition identities, and 12,500
relative translation controls. Its saved report binds the manifest hashes.
This separate verifier preserves the earlier singleton-tree proof directory.

Together the two directories prove the nine fixed-quotient exclusions in
Theorem 54: **5300 normalized assignments and 1328 certified representatives**.
They do not exhaust the full C5 family or all (10,10,10) supports.

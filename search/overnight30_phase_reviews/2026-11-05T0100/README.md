# Five further complete C5 phase exclusions

These are the five highest-scoring quotient classes left by the midnight
analytic expansion. The table indices refer to records in that expansion's
immutable `results.jsonl.gz`, not to all C5 graphs. `targets.json` preserves
and hashes the original residual list; the phase manifests bind that file.

| Expansion index | Normalized triangle-free assignments | Unit representatives | Exact uniform theta upper bound |
|---|---:|---:|---|
| 1 | 2778 | 695 | 722377013729/250000000000 |
| 3 | 288 | 72 | 92252147651/31250000000 |
| 4 | 612 | 154 | 582133481299/200000000000 |
| 5 | 5364 | 1342 | 1483513397801/500000000000 |
| 10 | 20 | 6 | 2866089290867/1000000000000 |

Every cyclic phase of each displayed quotient is covered. There are **9062**
normalized assignments and **2269** representatives under common units. These
are coverage counts, not pairwise graph-isomorphism counts. All representatives
have exact integer-factor theta certificates; the largest bound is
**1483513397801/500000000000 = 2.967026795602 < 3**. Theta here concerns the
exclusivity complement of the saved allowed nonzero-overlap support.

The existing proper-subset spanning-tree gauge and verifier in
`event30_c5_general_phases.py` are reused unchanged. Numerical discovery used
one bounded worker at a time and finished in 82.87 seconds under a 90-second
aggregate budget. Some numerical solves reported an inaccurate status; every
exclusion nevertheless passed the independent exact integer certificate check.
No numerical failure or timeout is used as an exclusion.

Standard-library verification from the repository root:

```sh
python search/overnight30_phase_reviews/2026-11-05T0100/verify.py
python search/overnight30_phase_reviews/2026-11-05T0100/independent_coverage.py
```

The second command uses the preceding independently developed fixed-order
residue-set enumeration, rather than the MRV bitmask discovery algorithm. It
checks every dual against independently reconstructed graph edges, all unit
compositions, and all 3125 relative translations of one representative per
class. These proofs require an actual simultaneous C5 automorphism on the six
orbits, but assume no symmetry or uniformity of quantum vectors or weights.
They do not exclude the entire C5-double family.

# Eleven further complete C5 phase exclusions

The source list is the **213 residual quotient classes** from the completed
01:00 review of the fixed pool of **5,589 classes**. `targets.json` reconstructs
that list from the two immutable analytic result archives and removes the five
01:00 phase exclusions. These targets are ranked by the stored review
representatives; later published phases of the same quotient can score higher.

The bounded batch selected target indices 0–6 and 9–12: a leading group plus
small phase enumerations. Every phase of all eleven selected classes is now
exactly excluded. None of the previously certified classes was repeated.

| Target index | Normalized triangle-free assignments | Unit representatives | Exact uniform theta upper bound |
|---|---:|---:|---|
| 0 | 2784 | 697 | 2970481016859/1000000000000 |
| 1 | 926 | 232 | 2849042154913/1000000000000 |
| 2 | 1224 | 306 | 739168347673/250000000000 |
| 3 | 330 | 83 | 1430143582057/500000000000 |
| 4 | 408 | 102 | 2960137384731/1000000000000 |
| 5 | 3576 | 894 | 587830416027/200000000000 |
| 6 | 1056 | 264 | 2949100968581/1000000000000 |
| 9 | 240 | 60 | 1414018822703/500000000000 |
| 10 | 120 | 30 | 2824291745289/1000000000000 |
| 11 | 1056 | 264 | 2821806718979/1000000000000 |
| 12 | 10 | 3 | 564220428159/200000000000 |

**11,730 normalized assignments and 2,935 unit representatives** are covered.
The largest exact upper bound is **2.970481016859 < 3**. These are phase
coverage counts, not pairwise graph-isomorphism counts. Theta concerns the
exclusivity complement of the saved allowed nonzero-overlap supports.
This upper bound is loose: its graph's saved numerical optimum is about
2.64945317. The highest numerical score in this batch is about 2.85216582.
`upper_bound_comparison.json` also checks earlier phase bundles: their largest
upper bound 2.991515946933 belongs to a graph scoring about 2.73287287, not to
a near-GHZ numerical candidate. The strongest saved phase graph across the
bundles scores about 2.89988931 and passes the current exact fast screens.

Discovery took **104.60 seconds**, below the 120-second aggregate budget, with
one numerical worker at a time. Every numerical certificate was reconstructed
and checked with integer arithmetic. Failed or inaccurate numerical results
are never accepted as exclusions. The phase enumeration and certificate
algorithms remain the unchanged version-one implementations.

The combined catalogs exclude **5,387 of the fixed 5,589 classes**, leaving
**202 unresolved**. The subsequent pinned chunk 87 contributes another 34
quotient classes, so the updated published pool has **5,623 classes**, of which
**236 remain outside these catalogs**. This is not a whole-family classification.
The incremental census is recorded separately from the fixed-pool proof coverage.

Standard-library proof verification:

```sh
python search/overnight30_phase_reviews/2026-11-05T0200/verify.py
python search/overnight30_phase_reviews/2026-11-05T0200/independent_coverage.py
python search/overnight30_phase_reviews/2026-11-05T0200/verify_pool.py --data-repo /workspace/ghz30-overnight
```

The independent fixed-order residue-set audit checks all exact graph-bound duals,
187,680 unit compositions and 34,375 relative translations. The discovery
algorithm instead uses MRV bitmask enumeration. Actual simultaneous C5 graph
automorphisms on six five-event orbits are required; neither vector symmetry
nor uniform quantum probabilities is assumed.

An exact whole-orbit closure prescreen found no Corollary 57b exclusion among
the 213 targets. Catalog integration, compact cut replay and seed provenance
are recorded in `integration_checks.json`. Raw launchers, logs and diagnostic
candidate lists are in the experiment archive named in `raw_data.json`.

The catalog has 249,390 permuted signatures and 72 count auxiliaries. All 528
new-class relabeling/SAT/compact-replay controls and eleven altered-hash
controls pass. The old 2.8521658155 seed is now excluded across every phase.
The refreshed guide has numerical theta 2.8549236048 and an individual exact
upper bound 57298781491/20000000000 < 3, already published with chunk 71. It is
retained only for mutations; the whole quotient remains unresolved. Its source
is verified against experiment commit `3536304` and record 50212 of chunk 71.

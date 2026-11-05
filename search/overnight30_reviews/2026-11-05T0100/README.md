# Review at 01:00 UTC, 5 November 2026

This review fixes the data repository at published commit
`5a2b45478d066d5180d1364259e1c7b1d233fc33` and compares completed jobs with
the midnight review at `5484ce3`. The nineteen new completed jobs are 61–79.
Running jobs and auxiliary benchmark archives are excluded from the main
coverage counts. The review does not modify the live search or data repository.

## Historical and active leaders

The asymmetric worker improved its historical numerical maximum to
**2.9166384145141824**, graph ID `0711be2a3c6aff655ef8`. Theorem 53 now
excludes this support directly: context B partitions into vertices 10–14,
whose A-neighborhoods are pairwise disjoint, and vertices 15–19, whose
C-neighborhoods are pairwise disjoint. Their total probability is bounded by
`15/32 + 15/32 = 15/16 < 1`, contradicting context normalization. This is an
exact exclusion, independent of the numerical theta score.

Job 69 retired that leader, preserving its graph, certificate, and historical
score. Job 79 retains it in `retired_best_records`; the snapshot and source
checkpoint are bound by the verifier. The highest active construction guide
at the cutoff is **2.9073614218247865**, graph ID `1a41a2ad0d88e6f018bc`.
It passes the configured structural screens, but the supplemental exact
review gives theta upper bound **7271/2500 = 2.9084 < 3**. Thus neither graph
is a GHZ realization. The lower active score reflects a valid filter upgrade,
not loss of the historical maximum.

The historical leader has full automorphism group order two and
context-preserving group order one. The active guide has trivial full
and context-preserving automorphism groups. `check_automorphisms.py`
reproduces both exhaustive VF2 counts; these separate diagnostics require
NetworkX. Both graphs come from the asymmetric construction worker, which
does not itself guarantee a trivial automorphism group.

All theta values here concern the exclusivity complement of the stored
allowed nonzero-overlap support. Numerical optimizer ranks are not physical
realization dimensions.

## Coverage and exact review

The newly completed archives contain 422,441 records, including 46,427 scored
records. Their 422,146 distinct graph6 labels include **421,977 labels absent
from the preceding completed pools**. There are **46,053 newly scored labels**,
out of 46,195 distinct scored labels in the new chunks. These are labeled
support counts, not counts of graph isomorphism classes.

Selecting the highest 120 newly scored labels plus each family maximum gives
128 distinct selected labels. Exact uncolored graph-isomorphism comparisons
within the selection and against previously exactly reviewed representatives
yield 111 new representatives and 17 equivalences. All saved equivalences
include explicit vertex permutations. No whole-pool isomorphism census was
performed; a representative may match an older unreviewed support.

Twenty selected representatives received exact exclusions. The twelve leading
new asymmetric representatives all have disjoint-neighborhood partition
certificates: eleven use budget 15/16 and one uses 1411/1440. Other family
leaders are excluded by rational theta duals or one-sided regularity. The
C5-double family maximum reuses its published integer-factor dual, with its
original compressed source and hash preserved. Five earlier certificates
are referenced by path and content hash rather than copied.

Including the supplemental active-guide review and one numerical timeout,
this directory contains **23 exact certificates**: 22 global exclusions and
one additional dimension-19 support obstruction. The selected-set count of
20 excludes those two supplemental global checks. Numerical discovery took
about 2.30 seconds in the collector and 1.80 seconds for the active guide,
using one numerical thread. Every exclusion is established by exact replay.

## Solver outcomes and filter behavior

The new archives contain seven ordinary `deadline before numerical solve`
records and one `numerical solve timeout`. The timeout, from job 73 and graph
ID `3621d670f46eaadcc482`, is independently excluded by Theorem 56's one-sided
regularity condition. No failed numerical solve is treated as a proof.
The scored records report 46,311 `optimal` and 116 `optimal_inaccurate`
statuses; neither status alone certifies a quantum realization or exclusion.

The C2-twopairs and C3-onecycle jobs generated tens of thousands of supports
without SDP scores because the exact filters rejected their samples. Their
principal rejection counters are degree budgets, probability propagation,
and bridge/partition-related cuts. This does not establish that either whole
family is impossible. Job 74 reports the already excluded C5-mixed family
exhausted without a SAT proof log; the earlier analytic theorem supplies that
family exclusion.

`job_audit.jsonl` separates new archive counts from cumulative checkpoint
counters and explicitly named batch fields. `best_retirement_audit.jsonl`
checks active and retired records against the current configured filters.
The job 79 history snapshot demonstrates that full historical records remain
available after retirement.

## Verification

Run from the source checkout:

```sh
python search/overnight30_reviews/2026-11-05T0100/verify.py
python search/overnight30_reviews/2026-11-05T0100/verify.py --inputs
/workspace/threecontext-venv/bin/python search/overnight30_reviews/2026-11-05T0100/check_automorphisms.py
```

The first command uses only the Python standard library to verify every local
certificate, graph binding, saved isomorphism permutation, referenced older
proof, and artifact hash. The second additionally checks immutable input blobs,
published factor-dual provenance, and checkpoint-to-history snapshot equality.
The verifier refuses optimized Python because its exact checks use assertions.

The global event and dimension bounds remain unchanged.

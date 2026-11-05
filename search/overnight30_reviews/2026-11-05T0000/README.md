# Review at 00:00 UTC, 5 November 2026

This review fixes the data repository at published commit
`5484ce3069b6073ad5bbedf63ee368dafb062c53` and compares completed jobs with
the 23:00 review at `f37683f`. The nineteen newly completed jobs are
41 and 43–60. Running jobs are excluded. The collection and exact checks
used one numerical thread; the first pass took about nine seconds. The
review does not modify the live search or its data repository.

## New leader

The asymmetric search improved its numerical maximum from approximately
2.91343416 to **2.916220709092425**. The new leader has graph ID
`5e8f40b8d72c15cfd969`; its completed-job-60 snapshot is preserved in
`snapshot_asymmetric_job60_best.json`. It survives the current combinatorial,
probability, regular-pair, and quotient filters, including the necessary
dimension-19 pair diagnostics. Its exact rational theta certificate proves

\[
\vartheta(H)\leq\frac{29173}{10000}=2.9173<3.
\]

It therefore remains excluded as a 30-event GHZ support. The twelve reviewed
asymmetric leaders are mutually nonisomorphic as uncolored graphs and are
not isomorphic to the previously exactly reviewed representatives. Their
numerical values range from 2.91459956 to 2.91622071; all twelve have exact
rational theta exclusions. This is evidence of new supports beyond the old
C5 leader, without closing the event or dimension gaps.

An exhaustive colored VF2 diagnostic finds context-preserving automorphism
group order one for each of the four strongest representatives. Thus these
supports have no nonidentity automorphism fixing each specified context.
This does not address automorphisms that interchange contexts. The results
are in `automorphism_audit.json`; rerun `check_automorphisms.py` in the
project environment to reproduce this separate exact combinatorial check.

## Coverage and scope

The new completed archives contain 312,222 records, including 61,550 scored
records. They contain 311,757 distinct graph6 labels, of which **311,477**
were absent from the earlier completed archives. There are **60,905 newly
scored labels**, out of 61,146 distinct scored labels in these new chunks.
These are labeled-support counts, not isomorphism-class counts.

Selection took the highest 120 new scored labels and each family maximum,
giving 130 distinct selected labels. Exact uncolored graph-isomorphism
comparisons within this selection and against prior exactly reviewed
representatives yielded 92 new representatives and 38 equivalences.
`isomorphic_candidates.jsonl` records a vertex permutation for every claimed
equivalence. No whole-pool isomorphism census was performed; a representative
may still match an older support outside the previously reviewed selection.

Twenty-two selected representatives received exact exclusion proofs. The
C5-double family maximum reuses its published integer-factor dual, copied
locally and bound to the immutable published source. Other family maxima
were excluded by exact theta or probability-propagation certificates. The
C2-onepair maximum also has a separate dimension-19 obstruction. Proofs for
equivalent older representatives are referenced by path and SHA-256 rather
than regenerated.

## Errors and filter behavior

There were eight routine deadline-before-solve records and one numerical
solve timeout. The timeout support, graph ID `f5e74eb4ad78f36b635f`, was
rechecked and has exact theta upper bound `13991/5000 = 2.7982 < 3`.
The scored records report 61,197 `optimal` and 353 `optimal_inaccurate`
statuses; those labels alone are not exact certificates.

Current filters exclude six historical active-best snapshots: four by
probability propagation and two copies of the former C5 leader by its
quotient obstruction. The later asymmetric best snapshots (jobs 50, 53,
56, and 60) all survive these filters. Job 55 has retired the old C5-double
leader and records a new active maximum of approximately 2.89072359, itself
excluded by the reused factor dual. The current retirement audit verifies
the exclusion of all 29 saved retired-best entries it encountered.

The newer C5-fixed and C2-onepair jobs 57 and 59 saved 81,575 and 44,507
models respectively without numerical scores, consistent with the stronger
filters eliminating those samples before SDP. Job 58 reports the already
excluded C5-mixed family exhausted; its status contains no proof log and is
not used here as an independent SAT proof.

The asymmetric checkpoint preserves the historical excluded-best value
2.91343416, but jobs 56 and 60 write a lower current-chunk graph into their
`best_excluded` file. This is a record-selection defect, not evidence that
the historical score was lost; it was reported to the worker maintainer.

## Verification

Run from the source checkout:

```sh
python search/overnight30_reviews/2026-11-05T0000/verify.py
python search/overnight30_reviews/2026-11-05T0000/verify.py --inputs
```

The first command uses only the Python standard library to verify local
certificates, graph bindings, saved isomorphism permutations, referenced
older proofs, and artifact hashes. The second also checks the published
input blobs and the provenance of the reused factor dual against the data
repository. Numerical discovery requires the existing project environment;
certificate replay does not. `job_audit.jsonl` explicitly separates new
archive counts from cumulative checkpoint counters and `batch_*` fields.

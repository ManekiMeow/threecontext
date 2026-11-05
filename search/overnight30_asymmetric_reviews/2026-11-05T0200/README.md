# Natural repair performance and paired-block control

This review is pinned to experiment commit
`3536304a360437b6c54413602377c31438c719f0`. It reads completed asymmetric
chunks 82, 85, 88, 92, 95 and 98, comparing exact labels with preceding completed
asymmetric main-worker pools. These are labeled-support counts, not graph
isomorphism classes or global novelty counts across all symmetry and auxiliary
pools. No new SDP or SAT search was run for this review.

The active numerical guide rose from 2.9073614218 before chunk 88 to
2.9113269650, 2.9141879080, 2.9147908752 and 2.9156526869 in the next four
asymmetric chunks. The historical excluded leader remained 2.9166384145.
The four new active leaders were directly logged as an edge kick and three
population crossovers. Complete ancestry was not recorded, so these improvements
cannot be causally attributed to the repair move or injected guide.

## Natural move observations

For the four completed chunks using partition-witness repairs:

| Move | Proposals | Accepted | Scored records | Labels unseen in prior asymmetric scores |
|---|---:|---:|---:|---:|
| Edge kick | 16,413 | 3,358 | 5,987 | 5,842 |
| Population crossover | 5,742 | 2,499 | 2,801 | 2,622 |
| Cluster partial rewire | 10,952 | 383 | 2,128 | 2,124 |
| Vertex rewire | 9,458 | 207 | 2,123 | 2,109 |
| Conflict-clearing switch | 5,910 | 92 | 1,832 | 1,830 |
| Partition-witness repair | 258 | 67 | 106 | 71 |

All 258 repair proposals were valid; none needed a fallback. Ninety-six of
their scored records survived the configured exact screens. Their highest
direct score was the already certified, excluded 2.9092320745 mutation guide.
These natural observations do not justify increasing repair frequency or
retuning the other moves. The table omits one saved-seed injection, which was
scored and accepted once.

Across all six reviewed chunks, excursions evaluated 5,177 candidates, scored
4,444 new records, and accepted 2,134 moves including 955 downhill transitions.
Their continuing exploratory behavior gives no reason to change the existing
temperature, length or cadence. Per-move solver times were not logged; elapsed
worker times cannot isolate a move's cost. Screens and seeds also changed
during this period, preventing a causal comparison of chunk throughput.

## A remaining exact obstruction

The natural repair logs contain 154 distinct labels. Exact replay excludes
80 by a disjoint-neighborhood partition of sizes **5+4+1**, with budget

`15/32 + 2/5 + 1/9 = 1411/1440 < 1`.

One additional label fails probability propagation; 73 survive these checks.
The 81 graph-bound certificates are in `natural_certificates.jsonl.gz`.
Screen survival is not evidence of a quantum realization.

A single new shared-neighbor edge in one five-event block leaves four events
with disjoint neighborhoods, while the other five-event block remains intact.
Deleting triangle conflicts cannot destroy that disjointness. Further edges
added during maximal completion can escape this residual certificate, but
the observed 80 exclusions motivate repairing two blocks before completion.

## Matched combinatorial control and implementation

Both procedures used the same frozen historical parent, original source commit
`fd05e9b`, and 128 matched deterministic seeds:

| Procedure | Valid proposals | Survived fast screens | Distinct labels | Exact rejections |
|---|---:|---:|---:|---|
| Repair one block | 128 | 55 | 91 | 73 disjoint-partition certificates |
| Repair two blocks | 128 | 119 | 127 | 9 probability-propagation certificates |

The complete control took 8.73 seconds including screening and deterministic
replays. No theta scoring was performed. This supports improved admission and
labeled diversity for this fixed parent, not improved fitness or a general
performance guarantee.

The worker retains its one-attempt-per-128-ordinary-proposals cadence and
alternates single- and paired-block attempts. Paired attempts preserve both
new sharing pairs, clear one edge per triangle conflict, try at most eight
target pairs, and fall back to the existing single-block move when needed.
Existing screening, temperatures and acceptance rules are unchanged. New
nested counters default safely for older checkpoints. The ordinary single
repair helper is unchanged.

Production controls cover 64 paired proposals and JSON RNG replays on the
original 5+5 parent and a real 5+4+1 witness, exact agreement with the approved
prototype, unchanged single moves, mismatched-certificate rejection, a real
one-block fallback, and the actual run-loop selector across old JSON checkpoints.
Independent review additionally checked protected edges, nontrigger and
excursion paths. No numerical solves were needed.

`raw_data.json` binds the immutable experiment archive in carrier chunk 104;
it contains all 256 controlled proposal records and their exact certificates.
These auxiliary records must not be counted as main-worker search coverage.

Verification from the repository root:

```sh
python search/overnight30_asymmetric_reviews/2026-11-05T0200/verify.py
/workspace/threecontext-venv/bin/python search/overnight30_asymmetric_reviews/2026-11-05T0200/worker_controls.py
```

The first command uses standard Python and pinned Git inputs. It replays
natural counters, all 81 exact natural exclusions, and all 256 controlled
proposals, checking every protected sharing edge and minimal conflict removal.
The second checks the current worker implementation without numerical solves.

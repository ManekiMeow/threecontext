# 22:00 UTC review

The strongest new support is the 95-edge `c5_double` graph with numerical
theta `2.9134341612056254`. It passes the current combinatorial screens and
the necessary dimension-19 pair-connectivity and matching tests. Its exact
rational dual certificate nevertheless gives

`theta(complement(H)) <= 1821/625 = 2.9136 < 3`.

The review read 69,112 scored records from nineteen completed jobs, plus an
atomic snapshot of the current job 23 best. These contain 68,486 distinct
**labeled** supports; resumed archives and best files can repeat records.
The selection consists of the top eighty labeled supports and each family's
strongest saved support. After overlap, these are 89 labeled graphs. Exact
uncolored VF2 isomorphism checks inside Weisfeiler-Lehman hash buckets reduce
this selection to eighteen representatives. This is not an isomorphism
census of the full pool.

Seventeen representatives received new exact exclusions: sixteen rational
theta duals and one same-context contraction witness. One additional
representative is the already-excluded initial asymmetric seed, with its
previous `(9,10,11)` coloring certificate referenced rather than duplicated.
The new certificate files are self-contained. The historical `/tmp` reference
in that one reused record is provenance, not required to verify any new result.

Nine of the newly screened representatives also fail a necessary condition
for dimension nineteen. The seven leading new asymmetric representatives,
with theta about `2.8873` to `2.8877`, all have exactly eighteen edges between
their second and third contexts. Their bipartite graph on twenty vertices
cannot be connected. The `c2_onepair` best similarly has only seventeen such
edges. This suggests the asymmetric search has concentrated on a graph
structure that is unsuitable for dimension nineteen; it is not a general
obstruction to thirty-event realizations in higher dimensions. The new
`c5_double` leader has pair edge counts `(30,35,30)` and avoids that restriction.

`strongest_candidates.jsonl` stores the selected representatives and results;
`isomorphic_candidates.jsonl` records the seventy-one exact isomorphism
duplicates within the selection. Family maxima in `summary.json` include
previously excluded mutation seeds, so they should not be read as viable
quantum candidates. The current job's complete atomic best record is preserved
as `snapshot_00023_c5_double_best.json`, with its SHA-256 hash in the summary.

The single-thread collection and certification pass took about 84 seconds.
The bounds on event number and dimension remain unchanged. Verify all new
certificates and the snapshot hash with standard Python:

```
python search/overnight30_reviews/2026-11-04T2200/verify.py
```

Adding `--isomorphism` also rechecks all seventy-one isomorphism claims from
their saved graph6 supports; that optional check requires NetworkX.

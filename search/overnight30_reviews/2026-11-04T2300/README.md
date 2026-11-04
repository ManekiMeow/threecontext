# 23:00 UTC published-pool review

This audit uses completed jobs in published commit
`f37683f803842de5b3a8d7fcffcbdd781bcff492`. Input bytes were checked against
their Git blob IDs; no data-worktree files were changed.

The strongest new selected support is asymmetric, with numerical theta
`2.9130510891771215`. It passes the current combinatorial and necessary
dimension-19 screens, but its exact rational dual gives

`theta(complement(H)) <= 29141/10000 = 2.9141 < 3`.

The overall numerical leader remains `2.9134341612`: the new labeled copies
are isomorphic to the C5-double support already excluded in the previous
review and by Proposition 46. The new asymmetric support is not isomorphic
to that leader or any of the other previously selected representatives.
The event and dimension bounds are unchanged.

## Coverage and selected representatives

Twenty newly completed chunks contain 108,351 archived records and 94,192
numerically scored records. After exact labeled deduplication against the
earlier completed chunks and selected snapshots, they add **107,777 labeled
supports**, of which **93,632 have numerical scores**. These are labeled
counts, not counts of graph isomorphism classes.

The top 160 newly scored labels plus each family's highest new scored label
give 169 distinct selected labels. Exact uncolored graph-isomorphism checks
identify 44 equivalences, including two to the prior C5-double leader, and
125 representatives distinct from one another and from the prior selected
representatives. No claim is made about isomorphism classes across the full
pool or novelty relative to previously unselected classes.

All 125 representatives and all 44 explicit vertex mappings are saved.
Thirty-three leading or family representatives received exact exclusions:
24 rational theta duals, eight degree-mass certificates, and one bridge
shortage. The other selected representatives were not certified in this
bounded audit. The proof files contain no numerical feasibility claims.

Two improved family supports, C2-fixed at `2.829880019` and C2-two-pairs at
`2.795056564`, were sent to the asymmetric worker as explicitly excluded
mutation seeds. Their exact degree budgets are `98/99` and `88/91`.

## Numerical failures

The archives contain eleven end-of-budget messages, one numerical timeout
in C2-two-pairs job 29, and one unreliable CLARABEL result in asymmetric job
37. The latter two supports were independently rechecked:

* The timeout support is excluded by a degree budget of `121/126`.
* The unreliable-result support has exact theta upper bound
  `5637/2000 = 2.8185`.

Thus neither unresolved numerical result conceals a feasible support among
these two records. The eleven routine deadline records were not re-solved.
There are also 1,581 C5-double scores marked `optimal_inaccurate`, alongside
92,611 scores marked `optimal` across all new chunks. These remain numerical
search scores; the saved exact exclusions do not rely on either status.
Archive counts describe the new chunks; checkpoint counters are cumulative
unless their names explicitly start with `batch_`.

## Files and verification

`summary.json` records counts and scope; `job_audit.jsonl` records per-chunk
archive activity and checkpoint counters. `strongest_candidates.jsonl`
contains the 33 certified representatives, `selected_representatives.jsonl`
contains all 125 selected representatives, and `isomorphic_candidates.jsonl`
contains explicit permutations proving the 44 equivalences.

`inputs.jsonl` binds the source files to the published commit. Existing proof
references are hash-checked, and rerunning collection reuses local exact
certificates. The first collection/certification pass took about ten seconds
on one CPU thread; the subsequent check reused those certificates.

Verify all 36 exact certificates, explicit isomorphism mappings, referenced
proofs, and artifact hashes using standard Python:

```
python search/overnight30_reviews/2026-11-04T2300/verify.py
```

Add `--inputs` to also recheck the published input bytes from the data
worktree or its Git object database. `collect_review.py` reproduces the
bounded selection using NetworkX and the search environment; it writes only
this review directory and should be run with one BLAS thread.

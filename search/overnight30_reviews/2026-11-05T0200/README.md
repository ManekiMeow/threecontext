# Review at 02:00 UTC, 5 November 2026

The immutable data cutoff is `3536304a360437b6c54413602377c31438c719f0`,
compared with the 01:00 review at `5a2b454`. Nineteen newly completed chunks,
80–98, are included. Running chunks and auxiliary research archives are not
counted as main-worker coverage. This review does not modify the live runner
or experiment Git index.

## Best candidates and exact findings

The highest newly scored support and active asymmetric guide has numerical
**theta 2.9156526868614114**, graph ID `b9a140d8fb027c569420`. Its exact rational
dual proves **theta <= 29167/10000 = 2.9167 < 3**. Full and context-preserving
automorphism groups are both trivial. It passes the configured structural and
dimension-19 support diagnostics, so it remains useful as a mutation guide;
it is not a GHZ realization.

The historical numerical maximum remains **2.9166384145141824**, already
excluded by the exact disjoint-neighborhood partition budget **15/16**. Its
full automorphism group has order two, while its context-preserving group is
trivial. The job-98 history snapshot retains that graph, score and proof.

All theta values concern the exclusivity complement of the saved allowed
nonzero-overlap support. Numerical optimizer ranks are not realization
dimensions. No event-number or dimension bound changed.

## Coverage and selected review

The new archives contain **271,321 records and 46,832 scored records**.
They contain 270,804 distinct labeled supports and 46,472 distinct scored
labels. Compared with all earlier completed main-worker pools, **270,430
labels and 46,181 scored labels are new**. These are graph6 label counts,
not graph-isomorphism counts.

Selecting the highest 120 newly scored labels and each family maximum gives
**127 selected labels**. Exact uncolored isomorphism comparisons within this
selection and against earlier exactly reviewed representatives give **124
new representatives and three explicit isomorphism mappings**. The result is
restricted to the selected sets, not a whole-pool classification.

Nineteen selected representatives have exact exclusions, including the twelve
leading new representatives and all newly selected family leaders. There are
**21 local exact certificates**: nineteen selected-graph exclusions plus two
supplemental solver-failure exclusions. All are global, with no additional
conditional dimension-19 certificate in this bundle. The C5-double leader
reuses its published exact integer-factor dual and retains its source path,
content hash and bound. Two earlier proof files are referenced by content
hash, including the historical leader's direct partition proof.

Numerical certificate discovery used **5.31 seconds** and eighteen attempts
with one numerical thread. Every exclusion is established by exact replay;
no numerical failure or bounded search is used as a nonexistence proof.

## Natural closure-filter audit

The new `degree_two_closure` filter produces **281 exact compact dominance
cuts** in these completed symmetry chunks: 50 in C2-twopairs, 225 in
C3-onecycle and six in C2-pairs. All 281 saved graph witnesses and generated
SAT clauses replay independently. Forty use the whole thirty-vertex set;
**241 use proper closed subsets of 23–29 vertices**. This demonstrates use
of the stronger closed-subset result, not only the whole-graph special case.
These are cut counts, not independent graph-isomorphism counts.

`closure_cuts.jsonl.gz` contains the exact compact witnesses, source archive
paths and zero-based source-line indices. `audit_closure_cuts.py --inputs`
checks the copied cuts against the pinned published blobs. The ordinary
verifier also replays all of these certificates.

The new archives report **two substantive numerical errors**: one solve
timeout in chunk 83, and one unreliable CLARABEL residual/status outcome in
chunk 85. The corresponding supports are independently excluded by an exact
theta dual and a probability-propagation certificate, respectively. Six
`deadline before numerical solve` records are normal budget endings. Scored
records report 46,784 `optimal` and 48 `optimal_inaccurate` statuses; neither
status alone proves realization or exclusion. Zero-score bounded chunks do
not exclude their whole symmetry families.

## Verification

From the source checkout:

```sh
python search/overnight30_reviews/2026-11-05T0200/verify.py
python search/overnight30_reviews/2026-11-05T0200/verify.py --inputs
/workspace/threecontext-venv/bin/python search/overnight30_reviews/2026-11-05T0200/check_selection.py
/workspace/threecontext-venv/bin/python search/overnight30_reviews/2026-11-05T0200/check_automorphisms.py
```

The first command uses standard Python and checks local certificates, graph
bindings, explicit isomorphism maps, referenced proofs, history, all 281 cuts
and artifact hashes. `--inputs` additionally verifies pinned blobs, reconstructs
the full old/new label and score counts, and binds every selected graph,
iteration/index and numerical score to its original archive record. The two
NetworkX diagnostics reproduce bounded-set nonisomorphism and the two leaders'
automorphism counts. Exact verifiers refuse optimized Python.

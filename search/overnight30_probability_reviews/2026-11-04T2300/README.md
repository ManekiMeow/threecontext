# Exact probability propagation review

This bounded review revisits the 125 selected isomorphism representatives in
`search/overnight30_reviews/2026-11-04T2300/selected_representatives.jsonl`.
It is not a census of the complete search pool.

The iterated strict-probability screen excludes 16 supports that passed the
previous combinatorial and degree-budget screens. Every new exclusion occurs
within three synchronous rounds, using upper bounds rounded upward to a
denominator of one million. The strongest has numerical theta about
2.9077694209, but its exclusion uses no numerical theta: the second round
forces the B-context probability sum below `993319/1000000`.

The 33 representatives already selected for exact SDP/combinatorial review
contain nine previous fast exclusions and 24 that still pass this probability
test. All 16 new certificates concern additional representatives. The updated
screen does not improve the global bounds of 30–33 events or 19–20 dimensions.

Verification uses only standard Python:

```sh
python search/overnight30_probability_reviews/2026-11-04T2300/verify.py
python search/event30_probability_propagation.py self-test
```

The proof and rectangle maximization are documented in
`search/event30_probability_propagation.py`. Actual upper bounds, every update,
and the terminal contradiction are replayed exactly. A passing graph is only
not excluded by this necessary test.

The eight-round discovery review took about 2.83 seconds for 125 supports on
the running search host. A separate deterministic reservoir sample of 32
records each from completed chunks 36, 39 and 40 found 46 further exclusions
beyond the old fast screen among 96 sampled supports. Three rounds took about
0.92 seconds for those 96 supports. These timing and sampling observations
justify a three-round production cap; they are not nonexistence proofs for
any search family. Exact sample counts, source identifiers and timings are
recorded in `summary.json`.

# Supplemental exact review of the 2.9202387 leader

The live chunk-108 best record was snapshotted at 02:25 UTC. Its graph,
`097316a813a1e997f5f6`, has numerical theta **2.920238739653219**, improving
the earlier historical maximum. A single numerical certification process,
limited to five seconds, completed its screening and exact dual check in
3.52 seconds. The resulting rational certificate proves
**theta <= 29213/10000 = 2.9213 < 3**. The graph is excluded as a GHZ
realization, while remaining useful as a mutation guide.

The graph passes the current structural, coloring and dimension-19 support
diagnostics. Its full and context-preserving automorphism groups are both
trivial, independently counted by exhaustive VF2 enumeration.

This was a carried best record, not a new chunk-108 discovery. Its saved
scoring record is in **chunk 101**, iteration **287180**, with direct move
**`population_crossover`**. The best timestamp is 02:03:35 UTC. This does not
establish ancestry of the crossover parents and does not attribute the graph
to the later paired-repair move.

The live snapshot is byte-identical to the already published chunk-101 best
record in experiment commit `d35f5efae769ed72ff6a575a807fd9a4a37e032c`.
`published_provenance.json` binds the scoring archive, SHA256 and zero-based
source line 6874. No new raw-data archive or data-worktree change is needed.
The frozen 02:00 audit is unchanged, and this graph is not added to its counts.

```sh
python search/overnight30_reviews/2026-11-05T0225_leader/verify.py --inputs
/workspace/threecontext-venv/bin/python search/overnight30_reviews/2026-11-05T0225_leader/check_automorphisms.py
```

Exact certificate replay uses standard Python and no numerical solver.
Theta concerns the exclusivity complement of the saved allowed support.
No global event-number or dimension bound changes.

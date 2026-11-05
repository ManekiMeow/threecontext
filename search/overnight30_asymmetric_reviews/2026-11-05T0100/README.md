# 01:00 witness-guided asymmetric construction

A bounded 19.634-second batch broke the historical leader's certified
`disjoint_partition` witness through shared-neighbor additions and minimal
triangle-conflict clearing. It processed 223 labels absent from the pinned
scored pool: 118 exact pre-score exclusions and 105 numerical scores. The latter
contain 55 nonisomorphic representatives within the batch. This does not count
isomorphism classes relative to the full previous pool.

The strongest new guide has numerical theta 2.9092320745 and trivial full
automorphism group. `best_exact_dual.json` proves the individual upper bound
29103/10000 < 3. The guide is excluded and is used only for subsequent mutations.
It improves the configured-screen surviving guide's score at the published
cutoff, but no global event-number or dimension bound changes.

`summary.json` binds the pinned input pool, raw archive hash and location,
construction coverage, exact controls, and the new sparse worker proposal.
`verification.json` records exact replay of all 118 pre-score certificates,
50 isomorphism maps, 55 pairwise nonisomorphic representatives, the best rational
dual, and 32 seeded worker repairs plus their deterministic replays. The complete
replay scripts and raw inputs are in the experiment archive named in the summary.
The normal worker weights, screening and acceptance rules are unchanged; one
certified historical-partition repair is attempted every 128 ordinary proposals.

Independent review added explicit singleton-block handling: verified singleton
blocks have no target context and are skipped when constructing repairs. The
original distinct-pair loop already skipped them, so this changes no behavior. A real
5+4+1 certificate of budget 1411/1440 passed 32 further deterministic proposal
controls. The additional control archive is referenced in `summary.json`;
the original raw archive and published review remain unchanged.

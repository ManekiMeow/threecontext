# Exact transfer of probability-certificate templates

This follow-up uses no numerical solver. It transfers the 48 distinct
coefficient/descriptor templates from the frozen 93-certificate probability
review to a separately frozen sample of nearby saved supports.

For a target graph, the old sector and subset rows are checked against its
actual allowed support. Their old weighted sum may now have residual
coefficients r. In each ten-event context, set lambda=max(r_i), add
`(lambda-r_i)` times `c_i < 1/9`, and subtract lambda times its normalization
equation. All added inequality multipliers are nonnegative, coefficients
cancel, and the right side changes by `(lambda-sum(r_i))/9`. A negative
right side, or zero with positive strict multiplier, gives a contradiction.
Multipliers are renormalized and the ordinary exact vertex-probability
certificate verifier checks the result independently.

The benchmark finished in **1.794 seconds**, below its 15-second budget:

| Sample | Processed | Already structurally excluded | New transferred certificates | Inconclusive |
|---|---:|---:|---:|---:|
| Previously LP-inconclusive controls | 4 | 0 | 0 | 4 |
| Previously unreviewed chunk98 labels | 120 | 6 | 52 | 62 |

Thus transfer certifies **52 of the 114 new-sample labels that pass the existing
fast screens**. Forty-four retain the old exact coefficient cancellation;
eight use nonzero residual correction, adding between one and ten upper-bound
rows. The highest-scoring transferred label has numerical theta **2.9091815205**.
These are additional individually excluded labeled supports, not new
isomorphism classes. The unsuccessful transfers make no feasibility claim.

Selection is the 120 highest-scoring distinct graph6 labels in pinned
experiment chunk98, after excluding every graph6 label in the saved previous
and current review selections. The 93 original certificate sources are also
excluded. Input metadata records every exclusion-source hash, the immutable
Git blob hash, source line numbers, tie ordering and the full sample. This
is a fixed local sample, not a full-pool or family classification.

All 52 new certificates replay exactly. The 48 templates also reproduce their
original certificates without corrections. To verify proofs, provenance,
benchmark counts and the frozen selection using standard Python:

```sh
python search/overnight30_probability_transfer/2026-11-05T0200/verify.py --inputs
```

An independent replay separates the measured costs: 0.319 seconds for cold
catalog hash/exact validation, 1.402 seconds for the 124 existing fast screens,
0.021 seconds for target geometry, 0.122 seconds for transfer arithmetic and
0.187 seconds for independent verification of the 52 emitted certificates.
Thus the incremental transfer cost, including full proof replay, is about
0.330 seconds over 118 fast-screen survivors. The prepared source templates
and verified moment bounds can be cached safely per process; emitted target
proofs still require exact verification. `timing.json` records this separate
replay. Six preexisting exclusion certificates are also saved and checked,
without modifying the frozen 52-proof discovery outputs.

The main probability review and its helper are unchanged. The benchmark above
preceded integration and changes no global event or dimension bound.

The asymmetric worker now calls `event30_probability_transfer.find_transfer`
after its existing exact screens. The 48 verified templates are cached once
per process; no target cache or numerical solver is added. Emitted proofs use
the self-contained global wrapper `vertex_probability`, with an exact
inner/outer graph-binding check. Symmetry-worker discovery is unchanged.

`integration_checks.py` verifies all 52 wrappers, exact H33/H36 scope controls,
hash and proof tampering, catalog-independent replay and historical record
preservation. A separate partition-repair parent is retained when a higher
excluded record uses a different method, including when an external historical
partition record is reclassified during recovery. The original single/paired
repair helpers and 128-proposal cadence are unchanged. The previously frozen
paired `worker_controls.py` passes all 64 proposal/replay/prototype controls,
eight single-repair controls and the actual selector/fallback checks without
editing its harness. Scoreless valid partition-parent fixtures remain supported.
The separate reports and hashes preserve the original discovery evidence.

The current numerical guide at theta 2.9202387397 remains inconclusive for
template transfer; its separate exact theta dual is retained. No numerical
solve was run for integration or these checks.

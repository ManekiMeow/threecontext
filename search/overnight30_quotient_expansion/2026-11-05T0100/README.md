# Completion of the published C5 quotient review at 01:00 UTC

The immutable experiment checkpoint is
`5a2b45478d066d5180d1364259e1c7b1d233fc33`. Its five completed C5-double pools
contain **34,721 scored records and 5,589 quotient classes**, identifying
context permutations and the two orbit swaps in each context. These quotient
counts are not graph-isomorphism counts.

The review skipped the original 150 selected classes and the 3,560 classes
reviewed at 00:00. It examined all **1,879 remaining published classes** in
31.61 seconds of a 90-second analytic budget, obtaining **1,843 exact
conditional quotient exclusions**:

- 534 sector-probability Farkas certificates;
- 502 affine and Cauchy-equality saturation certificates;
- 807 Farkas certificates using the disjoint-neighborhood moment bounds.

Thirty-six newly reviewed classes remain inconclusive. Together with the
earlier 182 inconclusions, there are **218 residual classes**. Before any
new finite-phase work, the combined analytic and earlier phase catalogs thus
exclude **5,371 of all 5,589 published classes**. Every published class has
now been reviewed. This is a complete review of the fixed published pool,
not an exhaustive classification of possible C5-double supports.

Each conditional proof requires an actual context-preserving group with
the specified five-event orbits; equitability alone is insufficient. Such
an exclusion covers every cyclic shift assignment with that quotient.
Floating-point optimization is used only to discover rational multipliers.
An absent certificate is inconclusive.

The five input archives are pinned by commit, path and SHA-256. A compressed
representative pool records every published quotient class, permitting an
independent exact check that the old and new selections partition the entire
fixed pool. Run with standard Python:

```sh
python search/overnight30_quotient_expansion/2026-11-05T0100/verify.py
python search/overnight30_quotient_expansion/2026-11-05T0100/verify.py --inputs
```

The second command independently reads the pinned Git blobs, reconstructs the
complete quotient census from all 34,721 scored records, and checks every saved
representative's source binding. The first replays all 1,843 new proofs and
verifies the actual C5 actions of the saved representatives.

The separate strict-probability trial uses Theorem 56 and exact zero-margin
Farkas certificates. It reviewed all 218 residual classes and found no new
exclusion; its failure to certify does not establish feasibility.

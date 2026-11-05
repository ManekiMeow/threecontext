# Exact bounds from disjoint cross-neighborhoods

For any subset of events in one context whose neighborhoods in another context
are pairwise disjoint, the saved rational moment certificates bound its total
state probability. They apply to arbitrary complex realizations, without a
symmetry, dimension, or equal-probability assumption. The proof uses the Gram
matrix of state-weighted event vectors and their projections onto the other
two contexts; see `search/event30_disjoint_probability.py`.

| Subset size | Certified total upper bound |
| --- | --- |
| 2 | 7/32 |
| 3 | 5/16 |
| 4 | 2/5 |
| 5 | 15/32 |
| 6 | 21/40 |
| 7 | 23/40 |
| 8 | 49/80 |
| 9 | 13/20 |
| 10 | 27/40 |

Each certificate consists of a rational 3-by-3 and 4-by-4 positive semidefinite
dual matrix and an exact affine trace identity. All saved 4-by-4 dual matrices
are positive definite, so the displayed bounds are actually strict. The
numbers are convenient certified bounds, not claims of optimality.

On thirty events, a six-event subset has probability greater than 5/9 by the
strict 1/9 bound on the other four events. Since 21/40 is below 5/9, six events
in one context cannot have pairwise disjoint neighborhoods in another.
Context partitions into such subsets of sizes 5+5, 5+4+1, or 5+3+2 also give
contradictions. Different blocks of a partition may use different other
contexts. The last partition uses strictness because its upper bounds sum to
exactly one.

For a true five-vertex automorphism orbit, weighted-Gram averaging makes its
event probabilities equal. Quotient column counts at most one toward another
context guarantee disjoint neighborhoods, so the orbit's per-event probability
is at most 3/32. Adding these valid upper bounds to the rational sector LP
excludes ten of the nineteen classes left by the earlier LP and saturation
catalogs: source indices 18,20,25,94,97,104,106,112,116,135. These certificates
cover every support with the stated actual group-orbit quotient, including
all cyclic phases. An equitable partition alone does not establish the
necessary symmetry.

Together the three probability catalogs exclude **141 of the 150 previously
selected quotient classes**. Nine remain unresolved by these probability
methods. Separate exhaustive phase certificates may provide further
exclusions. None of these selected-class counts is a global GHZ nonexistence
claim.

Verify using standard Python:

```sh
python search/overnight30_disjoint_reviews/2026-11-05T0000/verify.py
python search/overnight30_disjoint_reviews/2026-11-05T0000/independent_audit.py
python search/event30_disjoint_probability.py search/overnight30_disjoint_reviews/2026-11-05T0000/disjoint_05_mean_bound.json
```

The verifier checks all moment matrices, affine identities, quotient geometry,
Farkas multipliers, saved representative graphs, and a known H33 moment
control. Numerical optimization is used only to discover rational
certificates; numerical solver status does not establish an exclusion.
The independent audit uses a separate determinant expansion for all 198
principal minors and independently reconstructs all 54 affine coefficients.
The new certificate helper refuses optimized Python (`-O`), which disables
the assertions used for exact verification.

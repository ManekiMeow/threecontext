# Exact saturation exclusions for orbit quotients

This review applies rational affine equations and Cauchy equality conditions
to the 150 selected six-orbit quotients in the earlier probability-LP review.
It preserves that earlier review unchanged.

The exact saturation rules exclude five of the 24 quotients left unresolved
by the LP certificates. One is the already excluded Proposition 46 quotient;
the other four are new quotient-class exclusions. Together, the LP and
saturation certificates exclude **131 of these 150 selected quotient classes**.
The other 19 remain unresolved by these tests. This is neither a census of all
C5-double quotients nor a nonexistence proof for (10,10,10) GHZ graphs.

Each certificate applies to every support admitting the stated quotient under
an actual context-preserving group with two five-vertex orbits in each context.
An equitable partition alone is insufficient. All cyclic edge phases are
covered. Averaging the weighted Gram matrix may increase its rank; these are
exclusions in every dimension, not assertions about symmetry of a minimum-rank
realization.

All five final contradictions involve an edge whose saturated endpoints have
probability 1/10 but neighbor masses 1/5 and 2/5 toward each other. Equality
would require that edge's normalized overlap to be both 1/2 and 1/4. The exact
transcripts derive the necessary constant probabilities and masses first.

For quotient indices 5 and 6, the A-neighbor mass at B0 and C-neighbor mass at
B1 are both 1/5. Thus both B-orbit probabilities are at most 1/10, and context
normalization forces equality. Equality at B0 gives `x_C0+3*x_C1=2/5`, so both
C-orbit probabilities are also 1/10. C0 has A-neighbor mass 2/5 and B-neighbor
mass 1/5, while B0 has A-neighbor mass 1/5 and C-neighbor mass 2/5. The B0–C0
edge has the contradictory overlaps. No value of either A-orbit probability
is required for this contradiction.

Run the standard-library checks with:

```sh
python search/event30_quotient_saturation.py self-test
python search/overnight30_saturation_reviews/2026-11-04T2300/verify.py
```

The verifier replays all rational derivations, checks artifact hashes, and
checks the saved representative graphs' actual simultaneous C5 automorphisms
and quotient counts. Numerical theta values identify the source candidates;
they are not used in any proof.

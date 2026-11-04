# References

* `Yu_Simnacher_Nguyen_Guhne_PRXQuantum_3_010340_2022.pdf`: X.-D. Yu, T. Simnacher, H. C. Nguyen and O. Gühne,
  *Quantum-inspired hierarchy for rank-constrained optimization*, PRX Quantum **3**, 010340 (2022),
  doi:10.1103/PRXQuantum.3.010340. Open access, CC BY 4.0. Recommended by the user's mentor (30 Oct 2026).
* The paper this project follows up, arXiv:2208.07794 (Liu et al., Sci. Adv. 2025), is in `Perkel/2208.07794v3.pdf`.

## Why the rank-hierarchy paper matters here

**It is our problem.** Whether a four-context skeleton (RESULTS.md §3.12) is realizable is exactly the question of
§III B of the paper: does the non-orthogonality graph H on 4d rays have an orthonormal representation in dimension d?
That means a Gram matrix Γ ⪰ 0 with Γ_ii = 1, Γ_ij = 0 on the non-edges (which include all same-basis pairs), and
rank Γ ≤ d. With four complete bases, rank ≤ d already forces Γ² = 4Γ.

**The Lovász bound is the first level, and our candidates pass it.** The bound dim ≥ ϑ(H) is level N = 1 of the
paper's hierarchy. For the p = 97, …, 193 candidates, ϑ(H) = d exactly, so it certifies nothing.

**What the paper adds.** It writes the rank constraint as optimisation over separable states and adds symmetric
extensions and the PPT condition. That gives a complete hierarchy of SDPs, and infeasibility at any level is a
rigorous certificate of non-realizability, which is the direction we cannot reach numerically. Their Fig. 3 is an
11-vertex graph with ϑ = 4 whose N = 2 PPT relaxation already excludes dimension 4. Two practical points:
* the U(k) symmetry of the auxiliary space removes the rank k from the variable size (Eq. 25);
* k may be treated as continuous, which guards against numerical error (Appendix C).

**Cost.** Level N = 2 works with n² × n² matrices, n = 4p = 388–772 here, which is far beyond direct solution. The
authors' own N = 2 max-cut SDPs for n = 64 and 72 took 37 and 100 hours with SCS on one cluster node (their
ref. 52). Our skeletons, however, are Cayley graphs with large symmetry:
* the translations of Z_p split C^n ⊗ C^n into p blocks of size 16p, by total frequency;
* the multipliers H_e identify all but e + 1 of these blocks.

For p = 97 that leaves five distinct blocks of size 1552 (complex), about 10⁷ variables in all. That is a job for a
first-order solver with GPU eigendecompositions on a cluster node, not for this container; see the 2026-10-31
briefing for the estimate.

**Directly usable at small size.** The minimal realization dimension of the GHZ graphs H33, H36, H37 and H39
(RESULTS.md §1.6, n = 33–39) is a rank-constrained SDP of the same form, small enough for this machine once the
graph symmetry is used.

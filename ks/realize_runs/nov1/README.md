# Realizability re-examined, 1 Nov 2026 (RESULTS.md §3.13 rank-hierarchy audit, §3.14)

Run from `ks/` (the scripts add `ks/` to the path themselves).

* `../../unrealizable_examples.py`: the KS-uncolourable skeletons without a realization. It covers the Shrikhande
  pattern, `../../unrealizable_p41.json` and `../../unrealizable_p137.json`, and checks K4-freeness, the SAT colouring,
  (DEG)/(C2)/(C4)/(TRI), and local ϑ with exact rational dual certificates (about 4 minutes).
* `nearmiss_sat.py P E SECONDS OUT.json`: runs the CEGAR search of `cayley4.py` on F(p, e), but computes ϑ at all
  four bases for every SAT model and keeps the model with the largest smallest ϑ.
  * `nearmiss_sat_p41_e4.log`: 113 models in 184 s, then UNSAT. The best one, ϑ = 2.7924 at every basis, is
    `unrealizable_p41.json`.
* `nearscan.py E P...`: for every K4-free single-class pattern of F(p, e) passing `cyclo4.check`, the local ϑ at the
  four bases; it lists the near misses (ϑ < 3 somewhere).
  * `nearscan_e4.log`: p = 89–173.
* `local_face.py PATTERN [u...]`: the optimal face of the local ϑ SDP. For p = 97: rank 70, dual slack rank 2, real
  face dimension 288, and its projection to the moduli has dimension 69, i.e. it is unconstrained beyond the
  basis sums.
* `closed_relations.py PATTERN u`: the completeness relations that close inside {u} ∪ N(u). At p = 97 there are 12,
  the ones Proposition 27 uses.
* `smalld_cayley.py LIMIT ORDER[:indices]...` and `smalld_cayley_results.log`: Cayley skeletons with d ≤ 17. None has
  a model of the combinatorial encoding.
* `rankhier_logs/`: SCS logs of `../../rankhier.py` (level 2, PPT, complex).
  * The paper's Fig. 3 graph: k = 4, 4.005 (stall); 4.01, 4.02, 4.05, 4.1, 4.3, 5 (feasible).
  * The Shrikhande pattern: k = 4 (feasible), 3.5 (margin −0.0127).
  * `real_fig3_*.txt`: the real variant (`feasible_real`) on the Fig. 3 graph: k = 4.01, 4.05 feasible; k = 4 margin −8e-6 (inaccurate).

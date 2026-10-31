# Smaller three-context GHZ graphs, lower bounds, and the four-context Kochen–Specker problem

Follow-up to *Exploring the boundary of quantum correlations with a time-domain optical processor*
(arXiv:2208.07794; Sci. Adv. 2025). Everything numerical is reproducible from `search/`; every
claim marked **(exact)** is checked in rational arithmetic by `search/verify_all.py`,
`search/minimal37.py` or `search/dualcert.py`.

Notation (following Theorem 1 of the paper). `H` is the graph whose vertices are the events and
whose edges join **non-orthogonal** pairs. In the paper's figure this is the Perkel graph. The
exclusivity graph is `G = complement(H)`. A three-context GHZ paradox is the same thing as a graph
`H` with

* `H` triangle-free (`alpha(G) = 2`),
* `H` 3-colourable (`chi(H) = 3`; the colour classes are the contexts),
* `theta(complement(H)) = 3`.

Because `theta(Hbar) <= chi(H) = 3`, only the inequality `theta(Hbar) >= 3` needs a certificate.

**Scope of this file.** §1–2: three-context GHZ graphs (smaller examples than Perkel, lower bounds). §3: KS sets
covered by four cliques (four contexts). The separate question of the **smallest KS set in dimension 3** (bounds
[24, 31], exact searches over coordinate universes) is in [`ks3d/RESULTS.md`](ks3d/RESULTS.md), with code in `ks3d/`.
Its §7 explains why the spectral methods of §2 below cannot work there: every spectral bound is blind to the 3D KS
property, and the contradiction first shows at algebraic degree 4.

## Summary

| question | result |
|---|---|
| (1) smaller graph than Perkel (57 vertices) | **33 vertices, dimension 20** (exact, new; Section 1.7): H33, built from the affine plane AG(2,3), with the 9 points as one context and two copies of the 12 lines as the other two. Its realization is unique and has dimension 20, the smallest known. It was found from the equality case of Theorem 5, not by search. It is the smallest member of a family built from affine planes and nets (Section 1.8), and the operator lift of the qutrit Hesse configuration (Section 1.9). Previous record: **36 vertices** (exact), an 8-regular, vertex-transitive tri-Cayley graph over Z2×Z6 with integer spectrum {8, 2^16, 0^9, −4^10}. It is certified by the **unweighted** Hoffman bound and has a real realization in dimension 26 in which every overlap is ±1/4 and every event has probability 1/12. Also found: a 37-vertex graph that has realizations of rank 32 (and approximately 24), and a 39-vertex graph over Z13. See Section 1.6 on dimensions. |
| (2) minimality | **Theorem: every three-context GHZ graph has at least 30 vertices, and every realization needs dimension at least 18** (exact proofs, Section 2, Corollary 6; earlier bounds 28 and 17). Every event has probability ≤ 1/9 (Lemma 3). A context with exactly 9 events forces at least 12 in each other context and dimension ≥ 20 (Theorem 5); H33 is the only equality case. Two events of probability exactly 1/9 in one context already force ≥ 11 events in each other context, three force ≥ 12 (Theorem 21); so a smaller example has at most one such event per context, even approximately (Corollary 22). So **n_min ∈ [30, 33] and d_min ∈ [18, 20]**. A smaller example would need contexts of sizes 10–12 only. The 33-, 36- and 37-vertex graphs are vertex-critical (exact). The symmetric families searched (tri-Cayley graphs over all groups of order ≤ 12, all Cayley graphs of order 30 and 33) contain nothing below 36; H33 is not vertex-transitive and its contexts have unequal sizes. |
| (3) four-context KS set | Not found. Sufficient criterion: four orthonormal bases whose non-orthogonality graph is K4-free. Consequences in that setting: every ray is non-orthogonal to at least 30 rays in the other bases, and the dimension must be at least 18 (Corollary 6; the GHZ bound can never give more than 20, the dimension of H33). The K4-free case is the only one in which this bound applies. Numerical searches found nothing. **Reductions (Section 3.7):** rays outside every basis are redundant, so the problem is exactly "four disjoint orthonormal bases with a KS union". Only internal bases meeting ≥ 3 of the four can block a K4. A symmetry group transitive on all four bases forces the K4-free case. Four Pauli eigenbases are never KS. In the K4-free case every cross overlap is at most 1/3. **Computer-assisted: no example exists in d ≤ 5, so d ≥ 6** (Section 3.8, Theorem 15). **The orthogonality graph is never strongly regular, and a vertex-transitive one forces the K4-free case** (Section 3.10). **H33 cannot be used in the K4-free route** (two of its contexts have overlap 1/2 > 1/3; Theorem 19). Attempts to build a KS set on H33 in d = 20 all fail, and the structure behind this (four 5-dimensional class blocks, a duality W, a 10000-ray orbit of GHZ states) is in Section 3.11. **Symmetric families (Section 3.12, computer-assisted):** Cayley skeletons over Z_p that are invariant under the e-th powers are empty for 56 pairs (p, e) with p ≤ 151, including every e = 2, 3. **Multiplier families (Sections 3.12–3.13):**
* e = 2 and e = 3 give no example for any p (Theorem 24, Corollary 26).
* e = 4 gives exactly ten skeletons passing every necessary condition we have (Theorem 28; p = 97 to 457, 388 to 1828 rays).
* e = 5 gives 19 and e = 6 at least 111.

These conditions (K4-free, ϑ = 3 at every ray, ϑ(G) = 4, the flat spectral signature) are therefore far from sufficient. No realization was found numerically, and flat realizations are excluded for p = 97 by a counting proof (Proposition 27). The path forward is in Section 3.8. |

---

## 1. Smaller examples

### 1.1 Tri-Cayley graphs

For a group Γ and subsets S0, S1, S2 ⊂ Γ, let `T(Γ; S0, S1, S2)` be the graph with vertex set `Z3 × Γ`, where
`(a, g) ~ (a+1, g·s)` for every `s ∈ S_a`.

* It is 3-partite, and the classes `a = 0, 1, 2` are the three contexts.
* It is triangle-free iff `e ∉ S0·S1·S2`.
* Γ acts by automorphisms. For abelian Γ the Lovász SDP splits into |Γ| Hermitian 3×3 blocks (`search/ftheta.py`).
* The Perkel graph is `T(Z19; S, 7S, 49S)` with `|S| = 3`.

**Uniform weights.** If each `|S_a| = 3`, the optimal weights are forced to be uniform. A full scan (`search/tricayley333.py`)
shows that **q = 19 is the smallest cyclic group that works**. So the Perkel graph is the smallest 6-regular member of
this family, and smaller examples need larger connection sets.

**Complete scans** (all connection-set sizes, all S up to symmetry, exact SDP):

* Γ = Z_q with q ≤ 12: **nothing** (`search/fsall.log`).
* Γ = Z3², D5, Z2×Z6, D6, A4, Dic3 (all non-cyclic groups of order 9–12): exactly **one** graph up to isomorphism (`search/groups.log`). It appears both over
  Z2×Z6 and over Dic3, and it is H36 below.
* Γ = Z13: two graphs (39 vertices). Γ = Z14: several (42 vertices).

**Vertex-transitive graphs below 36.** If H is vertex-transitive, then `θ(H̄) ≤ χ_f(H) = n/α(H)`. Also `α ≥ n/3` because H is
3-colourable. So `α = n/3`, which forces `3 | n` and `deg ≤ α = n/3` (neighbourhoods are independent). Moreover
`θ(H) = n/θ(H̄) = n/3` forces `λ_min(A) ≤ −k/2` (Hoffman). Combined with Theorem 2 (n ≥ 28), the only VT orders below 36
are 30 and 33. An exhaustive scan of **all** triangle-free Cayley graphs on Z30, D15, Z3×D5, Z5×S3 and Z33
(`search/cayley_vt.py`, `search/cayvt.log`; best θ = 2.839) finds nothing. So **H36 is the smallest Cayley-graph
example**. Non-Cayley vertex-transitive graphs of order 30 and 33 were not checked.

### 1.2 The 36-vertex graph (previous record)

Let `H36 = T(Z2×Z6; S0, S1, S2)`, where elements of Z2×Z6 are written as pairs (a, z), and

    S0 = S1 = {(0,0), (0,3), (0,4), (1,4)},   S2 = {(0,1), (1,0), (1,1), (1,3)}.

* 36 vertices with contexts of size 12, 12, 12. 144 edges, 8-regular (4 neighbours in each other context).
  Triangle-free, girth 4, diameter 3. Every vertex has a unique vertex at distance 3.
* Adjacency spectrum: `8^1, 2^16, 0^9, −4^10`. Since `λ_min = −4 = −k/2`, Hoffman's bound gives
  `θ(H̄) ≥ 1 − 8/(−4) = 3` **(exact; the certificate is the 0/1 adjacency matrix)**.
* **Quantum realization.** Take the Gram matrix `I + A/4` (PSD of rank 36 − 10 = **26**). It is the only one (Section 1.6). The state ψ is the normalized
  sum of the vectors of any one context. Then every event has probability exactly 1/12, every non-orthogonal pair has
  overlap 1/4, and each context sums to 1.
* **(exact)** It is vertex-critical: by vertex-transitivity one deletion suffices, and a rational dual certificate gives
  `θ(complement(H36 − v)) ≤ 2.959`.
* Quantum–classical ratio 3/2, as for every three-context GHZ graph.

### 1.3 The 37-vertex graph (smallest dimension found: 24)

Let `H37 = T(Z13; {0,1,3,9}, {0,1,4}, {1,2,5,7})` with the vertices `(0,5)` and `(2,11)` removed.

* Contexts of size 12, 13, 12. 129 edges, degrees 6/7/8, triangle-free.
* **(exact)** `θ = 3`, via the rational certificate `search/shr_q13_1_0139_014_1257_cert.json`. This gives weights
  `c_i > 0` and `Y_ij` on edges with `B = Y + diag(c)` PSD, `B_ij = 0` on non-edges, and `sum(B)/tr(B) = 3`.
  The event probabilities are `c_i`.
* **(exact)** It is vertex-critical: every one-vertex deletion has `θ ≤ 2.9857` (`search/minimal37.py`).
* Unlike Perkel and H36, its optimal Gram matrix is **not unique** (Section 1.6). Walking inside the optimal face
  reaches an extreme point of rank **32** (non-edge residual 1e-10). A gradient search gives an approximate rank-**24**
  realization (orthogonality residual 3e-6, context sums 0.99997), which is not certified.

### 1.4 The 39-vertex graph over Z13

`T(Z13; {0,1,3,9}, {0,1,10}, {1,6,8})` has 130 edges and degrees 6/7.

* Certificate **(exact)**: `W = 1/4 A_{01} + 1/3 A_{12} + 1/3 A_{20}`. Every vertex has W-weight exactly 1 into each
  other class, so `W1 = 2·1`.
* `W + I ⪰ 0`, with spectrum of W equal to `2, 0.714^12, 0.244^12, −0.958^12, −1^2`. Hence `θ(H̄) ≥ 1 + 2/1 = 3`.
* All 39 events have probability 1/13. Deleting any edge orbit destroys the paradox.

### 1.5 Local search

* Randomized vertex deletion from all tri-Cayley hits never went below 36 or 37.
* Annealing over edge-maximal tripartite triangle-free graphs on 36 vertices (starting near H37) reached θ = 2.989 at
  best.

### 1.6 Is the realization dimension determined by the graph?

Any realization (vectors g_i, state ψ) gives an optimal point `B_ij = <ψ|g_i><g_i|g_j><g_j|ψ>/3` of the Lovász SDP.
Here rank B = dimension. The graph fixes only the **zero pattern** of B (zeros on non-edges). The edge entries and the
weights are free, subject to B ⪰ 0, tr B = 1, sum B = 3. So the dimension is fixed by the graph only when this optimal
face is a single point. In general the realizations form a face whose points have different ranks. The Cholesky rank
of the solver output is the **maximum** rank on that face, not the minimum.

`search/face2.py` computes the face around the interior-point optimum B₀, i.e. the affine space of `U M Uᵀ` with
`range U = range B₀` satisfying the constraints. It then walks to extreme points (numerical, tolerance 1e-7):

| graph | max rank (solver/Cholesky) | face dimension | lowest rank found |
|---|---|---|---|
| Perkel (57) | 37 | **0**, so the Gram matrix is unique | 37 |
| H36 | 26 | **0**, so it is unique: `I + A/4` | 26 |
| H37 | 35 | 92 | 32 at an extreme point; about 24 approximately |
| H33 | 20 | **0** (Theorem 5, step 6), so it is unique | 20 |

So for Perkel the answer to the paper's open question ("we do not know if a lower-dimensional realization exists")
is no among real realizations. **37 is forced**, because the optimal Gram matrix is unique, as expected for an
edge-transitive graph where symmetrization loses nothing. The same holds for H36 with 26. Complex realizations were
not analysed: a complex Gram matrix could add antisymmetric imaginary parts on edges. For H33 the uniqueness proof
(Theorem 5) covers complex realizations as well: 20 is forced.

### 1.7 The 33-vertex graph H33 from the affine plane AG(2,3) (new record)

This graph came out of the equality case of Theorem 5 (Section 2), not out of a search.

* **Contexts.** A = the 9 points of AG(2,3) = Z3². B and C = two copies of its 12 lines.
* **Edges (non-orthogonal pairs).** A point and a line (in B or in C) are adjacent iff the point lies on the line. A
  line b ∈ B and a line c ∈ C are adjacent iff they are distinct and parallel. That gives 96 edges, with degree 8 for
  points and 5 for lines.
* **Triangle-free.** A triangle would need a point on two parallel lines.
* **Realization (exact).** The Gram matrix has 1/√12 on point–line edges, 1/2 on line–line edges and 0 on non-edges.
  * Its spectrum is `3^1, 2^3, (3/2)^16, 0^13`, so it is PSD of rank **20**. The check uses the rational matrix
    `D⁻¹GD` and the identity `K(K − 3/2)(K − 2)(K − 3) = 0`.
  * The state is ψ = (1/3)Σ_A g_a = (1/√12)Σ_B g_b = (1/√12)Σ_C g_c.
  * Every point has probability 1/9 and every line 1/12, so each context sums to 1 and `θ(H̄33) = 3`.
* **Unique realization.** By Theorem 5 the realization is unique up to unitaries, so dimension 20 is forced. For
  comparison: H36 needs 26, Perkel 37, and H37's approximate 24 is not certified.
* **Vertex-critical and edge-critical (exact).**
  * Deleting a vertex leaves a proper 3-colouring with a class of at most 9 vertices, out of only 32. By Lemma 3 and
    Theorem 5 that is impossible.
  * Deleting a point–line edge violates step 3 of Theorem 5. Deleting a line–line edge violates step 6.
* **3-colourings.** Every proper 3-colouring has the 9 points as one class; the 16 colourings with a fixed colour for one
  point all have this shape. Every vertex sees both other colours in each colouring, as the rigidity argument of
  Section 2 requires.
* **Why the earlier searches missed it.** The contexts have unequal sizes (9, 12, 12) and the graph is not regular,
  so it lies outside the tri-Cayley and Cayley families of Section 1.1.
* Files: `search/h33.py` (exact verification, including the uniqueness step), `search/h33.log`, `search/h33.g6`.

### 1.8 A family: GHZ graphs from affine planes and nets (new)

H33 is the first member of an infinite family (`search/affine_ghz.py`, log `search/affine_ghz.log`).

* **Construction.** Take a net of order q with r parallel classes: q² points and rq lines, r ≤ q + 1, and r = q + 1 is
  the affine plane. The graph H(q, r) has three contexts:
  * A = the points; B and C = two copies of the lines;
  * a point and a line are adjacent iff incident; a B-line and a C-line are adjacent iff they are distinct and parallel.
* **The paradox is the parallel postulate.** A noncontextual assignment would be a point on two distinct parallel
  lines. So H(q, r) is triangle-free.
* **Realization.** Expanding ψ in each context forces the probabilities 1/q² (points) and 1/(rq) (lines), and the
  overlaps 1/√(rq) (incidences) and 1/(q − 1) (parallel pairs).
  * This Gram matrix is PSD, i.e. a GHZ realization, **iff r(q − 2) ≥ 2(q − 1)**.
  * The dimension is r + (q − 1)(2r + q + 1), minus r(q − 1) in the equality case.
* **Members:**

  | (q, r) | n | contexts | dimension |
  |---|---|---|---|
  | (3, 4) | **33** | 9, 12, 12 | **20** (H33) |
  | (4, 3) | 40 | 16, 12, 12 | 27 |
  | (4, 4) | 48 | 16, 16, 16 | 43 |
  | (5, 3) | 55 | 25, 15, 15 | 51 |
  | (4, 5) | 56 | 16, 20, 20 | 50 |
  | (5, 6) | 85 | 25, 30, 30 | 78 |

* **Failures, consistent with §2.** (2, 3) fails, because its points would have probability 1/4 > 1/9 (Lemma 3). (3, 3)
  fails too; it has 27 vertices, and Theorem 2 excludes 27.
* **Minimality within the family.** Only the two equality cases, (3, 4) and (4, 3), drop in dimension. H33 is the
  smallest member, since q = 3 needs all four classes and q ≥ 4 gives n ≥ 40.

### 1.9 Is there a link with the 3D KS problem? (new)

The two minimum problems have similar numbers: KS sets in ℝ³ need between 24 and 31 rays (`ks3d/RESULTS.md`), and
three-context GHZ graphs between 30 and 33 events. The numbers measure different things, and there is no reduction
either way.

* **Dimension separates them.**
  * A GHZ realization needs dimension ≥ 18 and contexts of ≥ 9 events (§2).
  * So no ray set in ℂ³ contains a three-context GHZ paradox, and no GHZ graph is realizable in ℂ³.
  * In particular, the Xu–Chen–Gühne reduction of a 3D KS set at a ray (§3.2) always yields all-versus-nothing proofs
    with many contexts of 2–3 rays, never three contexts.
* **"33 = 33" is a coincidence.** Peres's 33 rays in ℝ³ also split 9 + 12 + 12: the rays with entries in {0, ±1}, and
  two families of 12 with √2. But the 9 rays are orthogonal to 0 or 4 rays of one family and 0 or 2 of the other. An
  affine plane would need 4 and 4 for every point.

**The genuine link is with qutrit geometry (complex dimension 3), exact; `qutrit_picture` in `search/affine_ghz.py`.**
* **H33 inside qutrit operator space.** The 9 points of H33 span S₁, which is the operator space M₃(ℂ). The dictionary
  is:
  * **points:** the qutrit phase-point operators A_a/√3, where A_a = I − 2|f_a⟩⟨f_a| and the f_a form the Hesse SIC;
  * **lines:** the S₁-parts of the lines are the 12 MUB projectors P_ℓ/2;
  * **state:** ψ = I/√3, the maximally mixed qutrit state.
* **Orthogonality flips.** In ℂ³ the Hesse configuration has ⟨f_a|m_ℓ⟩ = 0 iff a ∈ ℓ. This is the orthogonality
  pattern of the 21-ray state-independent contextuality set of Bengtsson, Blanchfield and Cabello (arXiv:1109.6514).
  In operator space it becomes tr(A_a P_ℓ) = 1 iff a ∈ ℓ: incidence turns into *non*-orthogonality.
* **From inequality to paradox.** Those 21 rays are not a KS set: the SIC rays lie in no basis, so they can be coloured
  0. Lifting them, plus 11 extra dimensions that make each copy of the lines orthonormal, turns the same incidence
  structure into the smallest known all-versus-nothing paradox with three contexts.
* **Other dimensions.** The same dictionary holds for H(q, q + 1) with the phase space of a q-level system (checked for
  q = 3, 5, 7).
* **Why the 3 of "three contexts" selects the qutrit.**
  * Lemma 3's bound c ≤ 1/9 is a three-context bound. It excludes the qubit plane (points of probability 1/4).
  * Theorem 5 forces a context of exactly 9 events to be AG(2,3), the qutrit phase space.
  * So among phase-space constructions the qutrit is the smallest one with three contexts.

Conclusion: no structural link between the two minimum problems was found, and none is forced. The real 3D KS sets
(Conway–Kochen 31, Peres and Schütte 33) have cubic and octahedral structure, unrelated to the affine plane. The shared
"3" is the qutrit's: it connects H33 with qutrit state-independent contextuality (Bengtsson–Blanchfield–Cabello,
Yu–Oh type), not with the minimum KS set in ℝ³.

The methods do not transfer either (`ks3d/RESULTS.md` §7). Here the classical contradiction is free (H is
triangle-free) and the hard part is quantum realizability, which is what spectra measure. In 3D KS the roles are
reversed: every spectral or degree-≤3 algebraic bound is satisfied by the rays themselves.

---

## 2. Lower bounds (exact proofs)

**Setup.** Let `g_i` (`i ∈ V(H)`) be unit vectors with `g_i ⊥ g_j` for non-edges, let `ψ` be the state, and let
`I_1, I_2, I_3` be the colour classes. Each `I_k` is orthonormal and `Σ_{i∈I_k} |<ψ|g_i>|² = 1`, so `ψ ∈ S_k := span(I_k)`.

Let `Π_k` be the projector onto `S_k` and set `F = Π_1 + Π_2 + Π_3` on `span{g_i}` (dimension d). Then:

* `F = Σ_i g_i g_i†`, so `tr F = n`;
* `Fψ = 3ψ`;
* the spectrum of F lies in [0, 3].

**Lemma (triangle-free ⇒ cubic trace identity).** We have
`tr(Π_1Π_2Π_3) = Σ_{i∈I1, j∈I2, l∈I3} <g_i|g_j><g_j|g_l><g_l|g_i> = 0`. Every term vanishes, because a non-zero term would need a triangle.

Expanding `tr F³` over the 27 words in `Π_1, Π_2, Π_3` gives `tr F³ = 3 tr F² − 2 tr F`, that is,

    Σ_f f (f − 1)(f − 2) = 0   (sum over the eigenvalues f of F).

**Theorem 1 (n ≥ 27, and d ≥ 17).** Use the identity `f(f − 3/2)² = f(f−1)(f−2) + f/4`. Summing over the spectrum gives
`Σ_f f (f − 3/2)² = n/4`. The eigenvalue 3 (from ψ) contributes `27/4`, so

    n − 27 = 4 Σ_{f ≠ ψ} f (f − 3/2)²  ≥ 0.

For the dimension, the eigenvalue 3 contributes +6 to `Σ f(f−1)(f−2) = 0`. Every other eigenvalue contributes at
least `−2/(3√3)`, so at least `9√3 > 15.5` eigenvalues besides 3 are non-zero. Hence `d ≥ 17`.

**Theorem 2 (n = 27 is impossible, so n ≥ 28).** If n = 27, then every eigenvalue other than ψ's is 0 or 3/2, so
`F = 3/2·(1 + ψψ†)` on `C^d` (with d = 17).

* For `i ∈ I_1`, `F g_i = g_i + Π_2 g_i + Π_3 g_i`. So `Π_2 g_i + Π_3 g_i = g_i/2 + (3/2) a_i ψ`, where `a_i = <ψ|g_i>`.
* Pair this with `g_m` for an edge `i ~ m ∈ I_3`. The term `<g_m|Π_2 g_i>` vanishes (no common neighbour,
  by triangle-freeness). This gives `<g_m|g_i> = 3 a_i ā_m` for **every edge**.
* Expand `a_i` in the basis `I_3`: this gives `c(N_3(i)) = 1/3`, where `c_j = |a_j|²`. Any vertex with `a_i = 0` could be deleted,
  contradicting n ≥ 27, so all `a_i ≠ 0`.
* Then `<g_i|F g_i>` forces `c_i = 1/9` for all i. After rephasing so that every `a_i = 1/3`, the Gram matrix is
  exactly `I + A/3`, with A the 0/1 adjacency matrix of H.
* Its spectrum `{3, 3/2, 0}` means A has the eigenvalue 3/2. A rational non-integer cannot be an eigenvalue of an
  integer matrix. Contradiction. ∎

(Equivalently, H would have to be strongly regular with parameters (27, 6, 0, 3/2). This is consistent with the
paper's Theorem 3.)

(The numbering continues from Theorem 2 and skips 4, because §3 uses "Theorem 4" for the paper's Theorem 4.)

**Lemma 3 (one identity per event, so every event has probability ≤ 1/9).** Let `E = F − 3/2 − 3/2·ψψ†` and
`c_i = |<ψ|g_i>|²`. Then for every event i

    ||E g_i||² = (1 − 9 c_i) / 4.

*Proof.*
1. Let i ∈ I_1 and put p = Π_2 g_i and q = Π_3 g_i, so that F g_i = g_i + p + q.
2. The neighbours of i are pairwise non-adjacent (H is triangle-free), so their vectors are orthonormal and p ⊥ q.
3. We have `<g_i|p> = ||p||²` and `<g_i|q> = ||q||²`. Because ψ ∈ S_2 ∩ S_3, `<ψ|p> = <ψ|q> = <ψ|g_i>`.
4. Expand `||p + q − g_i/2 − (3/2)<ψ|g_i>ψ||²`. The terms in ||p||² and ||q||² cancel, and 1/4 − (9/4)c_i is left. ∎

Summing over i gives back `n − 27 = 4 Σ ||E g_i||²` (Theorem 1). Three consequences:

* `c_i ≤ 1/9`, so every context has at least 9 events of positive probability. This is the GHZ case of Lemma 14
  (§3.7), with a shorter proof.
* `c_i = 1/9` holds iff `F g_i = (3/2) g_i + (3/2)<ψ|g_i> ψ`.
* The Cauchy–Schwarz form of Lemma 14 gives more: `c_i (1/α + 1/β + 1/ρ) ≤ 1`, where α and β are the probabilities
  of the neighbourhood of i in the other two contexts and ρ = 1 − α − β.

`search/vertex_identity.py` checks the identity and the inequality on the realizations of H33, H36, H37 and H39. The
identity holds to within 2e-15. The maximum event probabilities are 1/9, 1/12, 0.094 and 1/13.

**Theorem 5 (a context with exactly 9 events).** Suppose a context A contains exactly 9 events of positive probability.
Then each of the other two contexts B and C contains at least 12 events of positive probability, and the realization
has dimension at least 20. If both B and C have exactly 12, the graph is H33 (§1.7) and its realization is unique up
to unitaries.

*Proof.* Delete the events of probability 0. This keeps a realization and can only lower the dimension. By Lemma 3
every a ∈ A has `c_a = 1/9`, so `F g_a = (3/2) g_a + (3/2)<ψ|g_a> ψ`. Rephase so that `<ψ|g_a> = 1/3`. Then
ψ = (1/3)Σ_A g_a.

1. *T = Π_2 + Π_3 on S_1.* `T g_a = F g_a − g_a = g_a/2 + ψ/2`. So S_1 is T-invariant, T = 1/2 on the 8-dimensional
   space S_1 ⊖ ψ, and Tψ = 2ψ.
2. *Two projections.* If Tv = λv with λ ≠ 1, then `<v|Π_2 v> = <v|Π_3 v> = λ||v||²/2`. (Expand
   `||Π_3 v||² = ||λv − Π_2 v||²`.) With λ = 1/2 and polarization, `<g_a|Π_2 g_a'> = δ_aa'/4 + 1/12`.
3. *Neighbours in A.* Take an edge a ~ b with b ∈ B.
   * Triangle-freeness gives `<g_b|Π_3 g_a> = 0`, so `<g_b|F g_a> = 2<g_b|g_a>`.
   * Comparing with `F g_a = (3/2) g_a + ψ/2` gives `<g_b|g_a> = <g_b|ψ>`.
   * Then `<g_b|ψ> = (1/3) Σ_A <g_b|g_a> = (|N_A(b)|/3) <g_b|ψ>`, so every b ∈ B has exactly 3 neighbours in A.
4. *Counting.* For a ≠ a', steps 2 and 3 give `Σ_{b ∈ N_B(a) ∩ N_B(a')} c_b = <g_a|Π_2 g_a'> = 1/12 > 0`.
   * So each of the 36 pairs of A has a common neighbour in B.
   * Each b ∈ B covers exactly 3 pairs, so |B| ≥ 12. The same holds for C.
5. *Dimension.*
   * Π_2 − Π_3 anticommutes with T − 1. So it maps the 1/2-eigenspace S_1 ⊖ ψ injectively into the 3/2-eigenspace of
     T.
   * That eigenspace is orthogonal to S_1. So Π_1 vanishes on it and F = 3/2 there, which gives 8 more dimensions.
   * The 17-dimensional sum is invariant under Π_2, and Π_2 has rank 1 + 8·(1/4) + 8·(3/4) = 9 on it.
   * The remaining |B| − 9 ≥ 3 dimensions of S_2 lie in its orthogonal complement. Hence d ≥ 17 + 3 = 20.
6. *Equality.* Suppose |B| = |C| = 12.
   * Each pair of A then has exactly one common neighbour in B. So the sets N_A(b) form a Steiner triple system on 9
     points, i.e. the 12 lines of AG(2,3), and step 4 gives `c_b = 1/12`. The same holds for C.
   * An edge b ~ c needs N_A(b) ∩ N_A(c) = ∅.
   * Rephase so that all `<g_b|ψ>, <g_c|ψ>` equal 1/√12. For a ∉ N_A(c), the non-edge (a, c) gives
     `Σ_{b ∈ N_B(a), N_A(b) ∩ N_A(c) = ∅} <g_c|g_b> = 1/2`.
   * These equations are solvable only if N_A(c) is a line of B's plane. The solution is then 1/2 on its two parallels
     (exact check over all 84 triples in `search/h33.py`).
   * So C's lines are B's lines, the B–C edges join distinct parallel lines with `<g_c|g_b> = 1/2`, and the Gram
     matrix is that of §1.7. ∎

**Corollary 6 (n ≥ 30 and d ≥ 18).** Every three-context GHZ graph has at least 30 vertices, and every realization
has dimension at least 18.

*Proof.* By Lemma 3 every context has ≥ 9 events of positive probability. If one context has exactly 9, Theorem 5 gives
n ≥ 33 and d ≥ 20. Otherwise every context has ≥ 10, so n ≥ 30.

For the dimension, the d − 1 eigenvalues f ≠ 3 satisfy `Σ f = n − 3` and `Σ f(f−1)(f−2) = −6`. (The eigenvalue 3 is
simple for d ≤ 33.) The convex envelope of `f(f−1)(f−2)` on [0, 3] is −f/4 on [0, 3/2] and `f(f−1)(f−2)` on
[3/2, 3]. So `(d−1)·env((n−3)/(d−1)) ≤ −6`. This allows the following n (`search/bounds30.py`, exact):

| d | 17 | 18 | 19 | 20 | 21 | 22 |
|---|---|---|---|---|---|---|
| largest n | 29 | 32 | 34 | 36 | 38 | 40 |

So d = 17 would need n ≤ 29, contradicting n ≥ 30. ∎

**Theorem 21 (rigid events).** Call an event *rigid* if its probability is exactly 1/9, i.e. E g = 0 (Lemma 3).
If a context A contains two rigid events, each other context contains at least 11 events. If A contains three, each
other context contains at least 12. A need not consist of rigid events; the counting in Theorem 5 is the case of nine.

*Proof.* Rephase every event x so that ⟨g_x|ψ⟩ = √c_x > 0, and let B be another context.
1. *Edges at a rigid event.* Step 3 of Theorem 5 uses only c_a = 1/9. So for every neighbour x of a rigid event a,
   ⟨g_x|g_a⟩ = ⟨g_x|ψ⟩ = √c_x. Expanding ⟨ψ|g_a⟩ in B then gives Σ_{b∼a} c_b = 1/3.
2. *Pairs of rigid events.* For rigid a, the vector g_a − ψ/3 is an eigenvector of T = Π₂ + Π₃ with eigenvalue 1/2.
   * The two-projection step of Theorem 5 needs only such eigenvectors, not a whole rigid context. So for rigid
     a ≠ a′ in A, ⟨g_a′|Π₂|g_a⟩ = 1/12.
   * By step 1 this reads Σ_{b∼a, a′} c_b = 1/12.
   * In particular, an event of B adjacent to two rigid events has probability at most 1/12.
3. *Counting.* Sort the events of B by the rigid events they are adjacent to. Every event has probability at most
   1/9 (Lemma 3).
   * Two rigid events a, a′: the classes "both", "only a", "only a′" and "neither" carry probability 1/12, 1/4, 1/4
     and 5/12. They need at least 1, 3, 3 and 4 events, so |B| ≥ 11.
   * Three rigid events: let t be the probability carried by the events adjacent to all three. The classes carry t,
     1/12 − t (each of the three pairs), 1/6 + t (each single) and 1/4 − t. Every t ∈ [0, 1/12] needs at least 12
     events; for example t = 0 gives 3 + 3·2 + 3. More rigid events only add constraints. ∎

Checks (`search/rigid_events.py`, log `search/rigid_events.log`):
* An exact LP and a MILP over the neighbourhood sets reproduce the minima: 11 for two rigid events, 12 for three to
  nine.
* The bound for three is sharp. In H33 every triple of points meets exactly 12 B-lines: 3 + 6 + 3 if the points are
  not collinear, 1 + 9 + 2 if they are.

**Corollary 22 (rigid and nearly rigid events below 33).** Suppose a realization has n ≤ 32 events of positive
probability, so every context has at least 10 events.
* No context contains three rigid events. A context with two rigid events forces n = 32 with context sizes 10, 11, 11.
  For n ≤ 31, every context contains at most one rigid event.
* *Stability.* There is an ε > 0 such that the same holds with "rigid" replaced by "probability above (1 − ε)/9".
  * For a fixed graph the realizations form a compact set, and there are finitely many graphs with at most 32
    vertices.
  * A sequence of counterexamples would therefore converge to a realization of the same graph in which those events
    are exactly rigid.
  * Theorem 21 applied to the limit gives the other contexts at least 11 (or 12) events of positive probability, and
    these stay positive along the sequence.
  * The ε is not explicit.

So the concentrated case, a context whose deficit sits on one or two events, is excluded. What is left is the
*spread* case: in every context, all events but at most one (two when n = 32) have probability at most (1 − ε)/9.
An example would be ten events of probability 1/10 each.

**Remark (towards the spread case).** For any events a, a′ of A, rigid or not, put w_a = g_a − ⟨ψ|g_a⟩ψ and
J = Π₂ − Π₃. Then exactly

    ⟨g_a′|Π₂|g_a⟩ = δ_aa′/4 + (3/4)⟨g_a′|ψ⟩⟨ψ|g_a⟩ + (1/2)⟨g_a′|E|g_a⟩ + (1/2)⟨w_a′|J|w_a⟩.

This follows from Π₂ = (T + J)/2 and from Ew = (T − 1/2)w on S₁ ⊖ ψ.
* |⟨g_a′|E|g_a⟩| ≤ (1/2)·min(√δ_a, √δ_a′), where δ = 1 − 9c is the deficit (Lemma 3).
* J maps the T-eigenvalue λ to 2 − λ, so |⟨w_a′|J|w_a⟩| ≤ √δ_a + √δ_a′ + √(δ_a δ_a′).
* Both corrections vanish for rigid pairs, which is step 2.
* The identity holds to 3·10⁻¹⁶ on H36 and H37, whose events are not rigid.

In principle this identity gives a quantitative Corollary 22, but the class probabilities in step 3 have slack of only about
1/36, so ε stays tiny. Excluding n = 30 and 31 this way would need ε > 1/9, because every 10-event context carries
total deficit 1. Excluding n = 32 would need ε > 1/8. The spread case needs a non-perturbative idea.

**Where minimality stands.** **n_min ∈ [30, 33] and d_min ∈ [18, 20]**, and H33 attains both upper bounds.

* By Theorem 5, a graph with fewer than 33 vertices has every context of size 10, 11 or 12, with n ∈ {30, 31, 32}.
* A realization in dimension 18 or 19 has no 9-event context (that needs d ≥ 20), so every context has ≥ 10 events,
  and n ≤ 32 or n ≤ 34 respectively.
* By Lemma 3 such contexts are "near-rigid": `Σ_{i ∈ I_k} ||E g_i||² = (|I_k| − 9)/4 ≤ 3/4`.
* By Corollary 22 their deficits are spread: at most one (nearly) rigid event per context, two in the 10-event
  context when n = 32.

The natural next step is the spread case: contexts of 10–12 events in which (almost) every event has
probability bounded away from 1/9. A plain
exhaustive search is not an option. The number of triangle-free graphs is already about 3·10^13 at 17 vertices (OEIS A006785),
and `search/rigid.c` handles only n ≤ 13.

**Rigidity and minimum degree (used for pruning).** In a vertex-minimal example every vertex has positive weight. By
the θ-argument (θ(H̄) = 3), every proper 3-colouring of H is a valid choice of contexts. Hence in **every** proper
3-colouring, every vertex has neighbours of both other colours. Otherwise moving it to another colour class would
force its weight to 0.

* In particular the minimum degree is at least 2, so vertices of degree 0 or 1 can be discarded.
* A vertex of degree 2 must have neighbours that receive different colours in every proper 3-colouring.
* By Lemma 3 such a vertex has small probability. With one neighbour in each other context, `c_i ≤ 7/135 ≈ 0.052`.
  In a context of 10 events the shortfalls 1/9 − c_i add up to only 1/9.
* No triangle-free graph on ≤ 11 vertices is rigid in this sense, and the rigid ones on 12–13 vertices have
  theta ≤ 2.5.

---

## 3. Four-context Kochen–Specker sets

**Update on the bounds.** Section 2 now gives n ≥ 30 and d ≥ 18 for three-context GHZ graphs (Corollary 6). Some
subsections below were written earlier and still say "≥ 28 rays" or "d ≥ 17" for the K4-free route (Route A). The
improved values ≥ 30 and d ≥ 18 apply there. The condition "d ≤ 16" (then every ray lies in a K4) becomes d ≤ 17.
These bounds come from a three-context GHZ graph at every ray, so they cannot exceed d_min(GHZ) ≤ 20 (H33). They say
nothing about Route B, where the best bound is d ≥ 6.

**Reduction.** Take four orthonormal bases `B_1, …, B_4` of `C^D` such that no 4-tuple (one vector from each) is pairwise
non-orthogonal. That is, their 4-partite non-orthogonality graph is K4-free. Then `B_1 ∪ … ∪ B_4` is a KS set.

*Proof:* a KS colouring must pick one vector per basis, and the picked vectors must be pairwise non-orthogonal.

The converse holds when the union contains no further bases. In general, extra bases inside the union could help,
and this route stays open (see below).

**Necessary local structure (K4-free case).** Fix `w ∈ B_4`. Its neighbours `N_k(w) ⊂ B_k` (k = 1, 2, 3) are
orthonormal, `w ∈ span N_k(w)` (since B_k is complete), and the non-orthogonality graph on `N(w)` is triangle-free.
So `(w, N(w))` is a **three-context GHZ paradox with state w**. By Section 2:

* every ray is non-orthogonal to **at least 30** rays of the other three bases (Corollary 6; 28 by Theorem 2);
* `D ≥ 18` (Corollary 6; 17 by Theorem 1).

So a K4-free four-context KS set would contain a three-context GHZ graph at every vertex. This
also recovers the paper's Theorem 4: shared vectors between bases are impossible in this setting.

**Searches (negative so far).**

* `search/k4opt.py` minimizes `Σ_{4-tuples} Π_{pairs} |<v_a|v_b>|²` over four real orthonormal bases, for
  D = 8, 12, 17, 20, 24 with 2 random starts each. The minimum stays at 37–50% of the Haar-random value and never
  approaches 0. This is weak evidence only: the landscape is highly non-convex, and any solution is non-generic.
* In the group-covariant ansatz `B_k = {U_g v_k}`, with an abelian G acting regularly, edges come from the Fourier supports of the
  unimodular functions `ū_k u_l`. Unimodular functions with small Fourier support are essentially functions on
  quotient groups. This makes that ansatz tensor-product-like, i.e. Bell-type, which the paper argues cannot work
  with fewer than four contexts. So a genuinely new structure is needed.

**Input from Lisoněk–Badziąg–Portillo–Cabello, PRA 89, 042101 (2014) and Cabello, arXiv:2011.13790.**

1. *The 21-ray, 7-basis set (d = 6) needs exactly 5 contexts* (`ks/ks21.py`, exact arithmetic in Z[ω]).
   * Its orthogonality graph has 105 edges and is exactly the line graph of K7, with no hidden orthogonalities. So its clique cover number is 5.
   * Any 4 orthogonal groups cover at most 18 of the 21 rays.
   * A 5-cover: B1, B2, B3 plus the two triangles {45,46,56} and {47,57,67} of the K4 on {4,5,6,7}, each completed to a basis with 3 new rays (`ks/ks21_cover5.py`).
   * Result: a 27-ray KS set covered by **5 contexts in d = 6**. This ties the Mermin star's 5 contexts in a lower dimension (6 instead of 8).

2. *Theorem 4 of the paper and the K4-free criterion assume that the four covering bases are the only bases.*
   * Both arguments put ψ = v₀ and use only the covering contexts.
   * A general KS set V ⊂ B1 ∪ … ∪ B4 can also contain **internal bases** that mix rays from several B_k. A KS colouring must hit each of these exactly once too.
   * A pairwise non-orthogonal transversal (v1, …, v4) is then not a colouring if some internal basis avoids all four rays. So a K4 is allowed, as long as every K4 is "blocked" by an internal basis.
   * This is exactly the mechanism of Cabello's true-implies-false sets (01-gadgets). It is also how the 7-context set works: its 7 bases are internal to any 5-cover.
   * So the general four-context problem is strictly larger than the K4-free one. Theorem 4's "disjoint" conclusion, and "no parity proof", are proven only in the K4-free setting.

3. *Dimension constraints in the general setting.*
   * A 4-context KS set has |V| ≤ 4d rays.
   * The smallest KS sets have 18 rays in d = 4 (proven minimal) and at least 24 rays in d = 3 (computer-assisted lower bound; see `ks3d/RESULTS.md`). Both exceed 4d, so **d ≥ 5**.
   * In d = 6 the only known KS set with ≤ 24 rays is the 21-ray set above, which needs 5 contexts.
   * Lisoněk et al. exhaustively list vertex-transitive "parity" candidates up to 31 vertices: 18 (d=4), 21 (d=6), 26 (d=4), 27 (d=6), 30 (d=4). Every one other than the 21-ray set violates |V| ≤ 4d.

### 3.1 Gadgets inside the four covering bases: what is possible (new)

**Gadget lemma.** Take a three-context GHZ paradox (ψ; I_1, I_2, I_3) and complete each I_k to an orthonormal
basis B_k. In any KS colouring of any set containing B_1, B_2, B_3, f(ψ) = 0. Proof: if f(ψ) = 1, every ray
orthogonal to ψ is 0, so the 1 of B_k lies in I_k. Those three rays are pairwise non-orthogonal, which is a triangle in H.
So each new GHZ graph (H36, H37, H39, Perkel) is a "**0-gadget for ψ that uses only 3 bases**". A four-context KS set is
then exactly a fourth basis B_4 all of whose rays are forced to 0, using only the internal bases of the union.

**Linking lemma (internal bases inside 4 bases).** Let B' ⊂ B_1 ∪ … ∪ B_4 be an orthonormal basis different from the B_k.
Then its pieces A_k = B' ∩ B_k span mutually orthogonal subspaces. In the two-basis case, W = span(A_k) = span(B_l \ A_l)
is spanned by subsets of **both** B_k and B_l. Each such common coordinate subspace W yields exactly two new bases,
and the classical constraint is only "B_k has its 1 inside W iff B_l does".

**Measured on the actual constructions** (`ks/ghzgeom.py`, `ks/ghz_gadget.py`, `ks/ghz_gadget2.py`):
* For H36 (d = 26): S_k ∩ S_l = span(ψ) exactly and dim (S_k + S_l)^⊥ = 3.
* For Perkel (d = 37): S_k ∩ S_l = span(ψ) and (S_k + S_l)^⊥ = 0.
* So the three GHZ contexts share **no** common coordinate subspace other than through ψ, and the only available links
  are B_4 ↔ B_k via S_k or S_k^⊥.
* Completing H36 with every link we could create gives 90–104 rays and 4–6 internal bases. The union is
  **KS-colourable** (SAT). A surviving colouring puts its B_4-one on a ray orthogonal to ψ; nothing forces those rays to 0.
* Conclusion: one GHZ gadget plus links forces only ψ. Each of the other d − 1 rays of B_4 needs its own forcing
  mechanism, and that mechanism must reuse the same three bases.

**Symmetric lifting = unimodular phase problem (lemma).** Suppose an abelian group G with |G| = D permutes all four
bases regularly, so the Gram matrix is G-invariant. Then each Fourier block of W is a 4×4 matrix with zero diagonal and
W² = 2W + 3I. This forces it to be 4|v⟩⟨v| − I with |v_a| = 1/2. So every basis is a G-orbit of a Fourier-flat vector with
phase function u_a, and rays (a, g), (b, h) are non-orthogonal iff h − g ∈ T_ab = supp FT(ū_a u_b).
* **Quadratic phases (stabilizer bases) never give a K4-free set.** The supports are cosets β_b − β_a + R_ab, so
  x_a = β_a is always a K4. This is why Pauli/Mermin-type constructions need a fifth context.
* **Cubic phases** over Z₂ⁿ, n = 5…8 (`ks/covariant4.py`, hill climbing on the exact K4 count): best K4 counts per vertex
  were 24, 96, 204, 214 (one search per n). None reached 0. A K4-free solution must also satisfy Σ_a |T_a4| ≥ 28 and
  D ≥ 17 (Section 2), which rules out n ≤ 4 and most of n = 5.

**Alphabet search, d = 6** (`ks/alphabet6.py`, `ks/ks4search.py`). All 1365 rays with entries in {0, 1, ω, ω²} and
their 1609 bases. Local search over 4-tuples of bases, counting all internal bases and KS colourings exactly by SAT:
the best unions contain 16 internal bases but still have ≥ 162 KS colourings.

**Two routes, two different dimension regimes.**
* **Route A: four bases are the only contexts (K4-free).** Every ray's neighbourhood is a three-context GHZ paradox.
  So d ≥ 17 and each ray is non-orthogonal to ≥ 28 others, **also over ℂ**: the trace identity uses only Hermitian
  projectors, and every term of tr(Π₁Π₂Π₃) vanishes. Every construction built from the GHZ graphs lies here or uses
  them as gadgets, which also needs d ≥ 17.
* **Route B: internal bases block the K4s.** KS demands *exactly* one 1 in every basis of the set. This is stronger
  than exclusivity, so neighbourhoods need not be GHZ paradoxes, and the bounds above do not apply. The only bound
  known here is d ≥ 5 (|V| ≤ 4d against known minimal KS sizes). Small-dimensional attempts (d = 6–8 alphabets) belong
  to this route and need small Cabello-type gadgets, not GHZ graphs.

**Refinement: when does a KS set contain a GHZ graph?** Let V be covered by B_1, …, B_4 and take w ∈ B_4.
* The rays non-orthogonal to w lie in B_1, B_2, B_3, and each restricted context has probability 1 in state w.
  So θ(G_w) = 3.
* If those rays contain **no** pairwise non-orthogonal triple, then α(G_w) = 2. G_w is then a three-context GHZ
  graph (Theorem 1) and **d ≥ 17**.
* Hence a four-context KS set with d < 17 must have, for **every** ray w, such a triple among the rays non-orthogonal
  to w. Every ray lies in a pairwise non-orthogonal 4-tuple, and all of these must be killed by internal bases
  ("Route B everywhere").
* Illustration (`ks/ks_neighbourhood.py`): in the 21-ray d = 6 KS set, the paradox derived from any ray w has
  10 events, clique cover 3, α = 2 but θ = 2.5. Its contradiction is a parity argument over 5 contexts, which the
  exclusivity graph does not capture. This is why KS sets exist far below the GHZ dimension bound.

### 3.2 Reversing the Xu–Chen–Gühne reduction (arXiv:2001.07656)

XCG turn a KS set into a GHZ-type proof. Take a ray ψ as the state and delete ψ and all rays orthogonal to it. The
remnants of the bases become *equality* contexts. Their GHZ-type proofs are more general than Theorem 1: their
minimal 10-event proof (from the 18-ray CEG set) has clique cover 4, α = 3, θ = 3.5 (`ks/xcg10.py`).

**Reverse-XCG criterion (proved).** S = B_1 ∪ … ∪ B_4 is KS **iff** for every u ∈ B_4 the XCG reduction of S at u is
classically infeasible. Proof: every KS colouring has exactly one 1 in B_4, and its restriction solves the reduction at
that ray. Conversely, a solution at u extends by f(u) = 1 and f = 0 on rays orthogonal to u.
* A GHZ graph (H36, H37, H39, Perkel) supplies the certificate at u = ψ.
* Reversing the reduction means supplying d − 1 further certificates, one for each u ∈ B_4 \ ψ, inside the same three
  bases. Route 1 gets them as symmetry images; Route 2 gets them from equality-type certificates built from
  internal bases.

**Abelian Route 1 with our graphs is blocked by cosets.** In a G-covariant construction, ψ's neighbours in B_k sit at
positions T_k ⊂ G with |⟨ψ|·⟩|² = 1/m each (uniform for H36, H39, Perkel). This requires a unimodular phase function
whose Fourier transform has constant modulus on T_k ("plateaued").
* If T_k = c_k + H is a coset of a subgroup, the neighbourhood Gram matrix depends only on H. It then has rank ≤ m
  (12, 13, 19), below the GHZ bound 17 and the rigid ranks 26 / 37 / 37. So this case is **impossible (proved)**.
* Non-coset plateaued supports: in Z_D (D = 26–48), random 12-element supports never admit such a function (residual
  ≥ 1; `ks/flatsupport.py`); only subgroup supports do. Non-coset plateaued functions do exist over Z₂ⁿ (Boolean
  plateaued functions), but there the support sizes are powers of 4, not 12, 13 or 19.
* The non-uniform example H37 (rank ≥ 32) is not covered by this argument.

### 3.3 Route 1 with non-abelian symmetry, and with H37

**General covariant setting.** Let a group Γ (|Γ| = D) act regularly on all four bases. Then B_4 = {δ_g} and
B_k = {λ_g b_k} with b_k a unitary element of ℂ[Γ]. The union is KS iff the neighbourhood of δ_e is
triangle-free. ψ = δ_e's class-k neighbours are λ_g b_k with b_k(g⁻¹) ≠ 0.

**Commuting-span obstruction (proved).** In any three-context GHZ paradox the context spans S_k cannot pairwise
commute. Otherwise F = Π₁+Π₂+Π₃ is diagonal with eigenvalues in {0,1,2,3}, and Σ f(f−1)(f−2) = 6·#{f = 3} ≥ 6 > 0,
contradicting the trace identity.
* **Consequence for Route 1, any group (abelian or not):** if the neighbour positions of each class are cosets of
  subgroups (b_k supported on a coset), the spans S_k are coordinate subspaces ℂ[·] of ℂ[Γ]. These commute, so the
  construction is **impossible**. This strengthens the abelian result of §3.2.
* **Our four graphs:**
  * H36, H39 and Perkel have S_k ∩ S_l = span ψ (`ks/intersections.py`).
  * H37 has context sizes 12, 13, 12, which forces trivial subgroup intersections.
  * Their natural placements (tri-Cayley over Z₂×Z₆, Z₁₃, Z₁₉, i.e. each class on a coset or a coset minus a point)
    collapse every class into one ℂ[H] of dimension ≤ 19, so they are excluded.
* Route 1 therefore needs **non-coset supports** for the b_k. Unitarity of b_k with support on a prescribed
  12–13-element set is then a "perfect sequence" condition, Σ_x b̄(x) b(xg) = 0 for all g ≠ e. Random supports never
  satisfy it (cf. `ks/flatsupport.py`), so targeted embeddings of H37 need structured non-coset supports. This is
  open.

**Generic non-abelian search** (`ks/nonabelian4.py`). Exact unitarity is enforced via a Newton–Schulz polar factor of
the convolution operator, and the smooth K4 count of the whole covariant 4-partite structure is minimized. Groups:
SL(2,3), S₄, D₁₂ (24), D₁₃ (26), Heis₂₇, S₃×S₃, A₄×Z₃ (36), Z₁₃⋊Z₃ (39, the symmetry of the Z₁₃ graphs), Z₇⋊Z₃×Z₂ (42),
S₄×Z₂ (48), Z₁₉⋊Z₃ (57, Perkel's symmetry) and A₅ (60).
* In every run the K4 loss drops only to 21–42 % of its start, and all supports stay full
  (`ks/na_*.log`). There is no sign of a sparse K4-free configuration.

### 3.4 Toward impossibility: Theorem 4 holds in general (new proof)

The paper's proof of Theorem 4 reasons with α and θ of the exclusivity graph left after removing v₀. That does
not account for internal bases, whose "exactly one" constraints the exclusivity graph cannot see (§3.1–3.2).
Here is a proof that includes them.

**Lemma (two-context colouring).** Let S be a set of rays and v₀ ∈ S. Suppose every ray non-orthogonal to v₀ lies in
C ∪ C′, the remnants (after deleting rays ⊥ v₀) of two complete bases contained in S. Then S has a KS colouring with
f(v₀) = 1.

*Proof.*
1. Let p(x) = |⟨x|v₀⟩|² on C and p′(y) = |⟨y|v₀⟩|² on C′. Each sums to 1, since v₀ lies in the span of each remnant.
2. For X ⊆ C, let Y ⊆ C′ be the rays orthogonal to all of X. Then X ∪ Y is orthonormal, so Bessel gives
   p(X) + p′(Y) ≤ 1, i.e. p(X) ≤ p′(N(X)).
3. By Hall's / the transportation theorem there is a coupling q(x, y) with marginals p, p′ supported on
   non-orthogonal pairs.
4. Any other basis B′ ∌ v₀ of S has remnant A ∪ A′ (A ⊂ C, A′ ⊂ C′, A ⊥ A′) of total probability 1 in state v₀. So
   P(x ∈ A) + P(y ∈ A′) = 1, and the two events are disjoint under q. Hence every (x, y) in the support of q hits B′
   exactly once.
5. f = 1 on {v₀, x, y} (pairwise non-orthogonal) and 0 elsewhere is a KS colouring. ∎

(Shared rays. If the two bases share rays, so that C ∩ C′ ≠ ∅, put mass p(z) on (z, z) for each shared z and
apply Hall to the rest. X ∪ Y ∪ (C ∩ C′) is still orthonormal, so step 2 goes through. In step 4,
E|{x, y} ∩ R| = p(R ∩ C) + p(R ∩ C′) − p(R ∩ C ∩ C′) = p(R) = 1 for the remnant R of B′. So the lemma, and
Theorem 4 below, need no disjointness assumption.)

**Consequences (internal bases allowed).**
* **Theorem 4 holds in general.** If a KS set is covered by four bases, they are pairwise disjoint (a ray in two of
  them has all its neighbours in the other two). So the union has 4D distinct rays.
* In a KS set covered by k bases, no ray lies in k − 2 or more of them.
* No KS set is covered by three bases. This extends the Xu–Chen–Gühne lemma, which assumed each vector lies in exactly
  one context.

**Necessary condition per ray.** Let w ∈ B₄ with neighbours in B₁, B₂, B₃ and probabilities p(w)_r = |⟨r|w⟩|².
Suppose p(w) is a convex combination of pairwise non-orthogonal transversal triples, i.e. p(w) ∈ STAB(G_w). Then every
internal basis is hit exactly once almost surely, so a KS colouring exists. Hence **in a four-context KS set every ray w
gives a state-dependent contextual point p(w) ∉ STAB(G_w)** for the other three bases.
* Route A (K4-free) satisfies this trivially, since there are no transversal triples.
* In Route B it is a genuine constraint.

**What an impossibility proof still needs.** A single ray cannot give the contradiction: the GHZ graphs are
maximally contextual for one state. The argument has to combine all rays, e.g. show p(w) ∈ STAB(G_w) for at least one
w, using Σ_{w∈B₄} p(w)_r = 1 for every r.

### 3.5 Orthogonal GHZ states for completed GHZ contexts (exact MILP)

Route A reformulated: three complete bases B₁, B₂, B₃ of ℂᵈ admit a fourth basis making a KS set iff there are
d mutually orthogonal **GHZ states**, i.e. states whose support in B₁ ∪ B₂ ∪ B₃ is triangle-free.
* A GHZ state is a point of ∩_k span(support in B_k); in all our examples this intersection is one-dimensional.
* So GHZ states are isolated, and descent or greedy searches miss them. Both heuristics failed the positive
  control (two block-diagonal copies of H36 have ≥ 2 orthogonal GHZ states), so their "1 found" results are void.

**Exact test** (`ks/ghzmilp.py`): a mixed-integer program over ψ ∈ ℂᵈ and binaries z_r (z_r = 1 ⇒ ⟨r|ψ⟩ = 0).
* Constraints: every triangle has a zero ray; ψ ⊥ ψ₀; ⟨r*|ψ⟩ = 1, enumerated over r* ∈ B₁.
* A unit GHZ state has some B₁-overlap ≥ 1/√d, so the big-M bound √d + 1 is safe and infeasibility for every r* is a proof.
* It passes the positive control (finds the second GHZ state of 2 × H36).
* Status codes: 0 = feasible, 2 = infeasible, 1 = time limit.

Results for the realization dimension d, completing each I_k by an ONB of S_k^⊥. Completions are random
(Haar), or symmetry-adapted (eigenvectors of the translation symmetry inside S_k^⊥):

| graph | d | completion | triangles | GHZ state ⊥ ψ₀? |
|---|---|---|---|---|
| H36 | 26 | random | 11816 | **none** (26/26 infeasible) |
| H36 | 26 | symmetric | 3410 | **none** (26/26 infeasible) |
| H39 | 37 | random | 39408 | **none** (37/37 infeasible) |
| Perkel | 37 | random | 27378 | none found (36 infeasible, 1 time limit) |
| H39 | 37 | symmetric | 5088 | undecided: r* ∈ completion rays (18) infeasible; r* ∈ I₁ undecided (30-min MILP limit) |
| Perkel | 37 | symmetric | 4122 | undecided: r* ∈ completion rays (18) infeasible; r* ∈ I₁ undecided (30-min MILP limit) |

So for the decided cases the largest orthogonal GHZ family is exactly 1 (ψ₀ alone), far from the d needed. Normalising on a ray of I₁ makes the big-M MILP much harder; those cases need a stronger formulation (e.g. symmetry breaking or exact support enumeration). This is
evidence for, not a proof of, an upper bound on orthogonal GHZ families. Remaining freedom: other completions,
and embedding in a larger space D > d.

### 3.6 Graph invariants alone cannot rule out four contexts

Route A in graph language: the orthogonality graph G on |V| = 4d rays is four disjoint d-cliques with α(G) = 3.
This forces ϑ(G) = 4 (cover by 4 cliques, and every state sums to 4) and ϑ(Ḡ) = d (the rays are an orthonormal
representation of Ḡ in dimension d, and ϑ(G)ϑ(Ḡ) ≥ n). It also gives χ(G), χ_f(G) ≥ 4d/3 > d automatically.

**Counterexample to a purely graph-theoretic proof.** The GHZ exclusivity graph (complement of the Shrikhande
graph) has n = 16, four disjoint 4-cliques, α = 3, ϑ(G) = ϑ(Ḡ) = 4 and χ_f ≥ 16/3 > 4. So it meets every condition
for d = 4. Yet no four-context KS set exists in d = 4, since every KS set has ≥ 18 rays (Xu–Chen–Gühne).
The obstruction is the rank: the Gram matrix P must satisfy P² = 4P with identity blocks. ϑ(Ḡ) = d is only its SDP
relaxation. Rank-sensitive arguments (the trace identity behind d ≥ 17 in Route A, the 18-ray bound in Route B)
are needed.

**Most promising routes now.**
1. Several GHZ states in the same three bases. One could look for a non-abelian symmetry (e.g. the multiplier
   group of Z13 acting on the H39 realization), mapping ψ to ψ' ⊥ ψ while permuting the completed contexts.
2. Cubic or higher-degree phase functions over larger or non-abelian groups (the covariant formulation above),
   searched with the K4 count as the objective.
3. (Route B) Larger alphabets ({0, ±1, ±i, ω^k}) in d = 6–8, using the exact SAT pipeline in `ks/ks4search.py`.

(Superseded by §3.8, which uses the reductions of §3.7.)

### 3.7 Rays outside the bases, incomplete cliques, and which internal bases matter (new)

Conventions: V is a finite set of rays in ℂ^d (d ≥ 3). A *clique* is a set of mutually orthogonal rays, and a
*basis* is a clique of size d (automatically a maximal clique). A KS colouring puts exactly one 1 on every basis
contained in V and never puts 1 on two orthogonal rays.

**Question: can a KS set contain a ray that lies in no orthonormal basis of the set?** Yes, but such a ray is never
needed. Also, the four cliques of a four-clique cover can always be taken to be complete bases.

**Lemma 7 (deletion).** Let V be KS and r ∈ V. If r lies in no basis contained in V, then V \ {r} is KS. The same
holds if r lies in exactly one basis B ⊆ V and every ray orthogonal to r lies in B.

*Proof.* Let f colour V \ {r}.
* First case: put f(r) = 0. The bases of V are those of V \ {r}, and a 0 never breaks exclusivity.
* Second case: B \ {r} is still a clique, so it carries at most one 1. Put f(r) = 1 − f(B \ {r}). Then B has exactly
  one 1, and if f(r) = 1, every ray orthogonal to r (all in B \ {r}) is 0.

Either way f colours V. ∎

* So in a critical KS set every ray lies in a basis and has at least d orthogonal neighbours. A ray in two bases has
  at least d + 1, because two distinct bases share at most d − 2 rays. For d = 3 this is the "minimum degree 3"
  property of Uijlen–Westerbaan (arXiv:1412.8544, Prop. 10).
* Check (`ks/reductions_check.py`): Mermin's star plus 7 rays in no basis (5 generic, 2 orthogonal to some star
  rays) is KS, and it stays KS after they are deleted. Peres' 33-ray set has no ray outside its 16 triads.

**Lemma 8 (monotonicity).** Any finite superset of a KS set is KS. A colouring of the superset restricts to a
colouring of the subset, because orthogonality is intrinsic and every basis of the subset is a basis of the superset.

**Corollary 9 (four cliques = four disjoint bases).** Suppose the orthogonality graph of a KS set is covered by four
cliques C₁, …, C₄ (the paper's "four context-cover"); the cliques need not be bases. Complete each C_i to an
orthonormal basis B_i. By Lemma 8, B₁ ∪ … ∪ B₄ is KS. By Theorem 4 (§3.4) the B_i are pairwise disjoint.
* So **a four-context KS set exists iff some ℂ^d has four pairwise disjoint orthonormal bases whose union (4d rays) is
  KS.** Rays outside the four bases never help, and assuming complete contexts, as the paper's Theorem 4 does, loses
  nothing.
* The other reading, "the orthogonality graph has only four maximal cliques", is impossible:
  1. Delete the rays that lie in no basis (Lemma 7). What remains is covered by those of the four cliques that are bases.
  2. Three or fewer bases are never KS (§3.4).
  3. Four bases are disjoint, so every orthogonal pair lies inside one B_i, and one ray from each B_i is a colouring.
* Critical cores can still have incomplete covering cliques. The critical 21-ray set needs five cliques, two of them
  triangles in d = 6 (§3, item 1); only its completion is covered by bases.

**Lemma 10 (common coordinate subspaces split the problem).** Let W be a proper nonzero subspace spanned by subsets
of all four bases. Then V is KS iff V ∩ W and V ∩ W^⊥ are both KS. Each of these is again a union of four disjoint
bases, of W and of W^⊥ respectively.

*Proof.* Every ray of V lies in W or in W^⊥. So every basis B′ ⊆ V splits into a basis of W and a basis of W^⊥, both
inside V, and every such pair is a basis. Rays of W are orthogonal to rays of W^⊥, so a colouring has all its 1s on
one side, where it restricts to a colouring. Conversely, a colouring of one side, extended by 0, colours V. ∎

So in the smallest dimension d_min with an example, the four bases share no proper coordinate subspace. An example
that has one satisfies d ≥ 2·d_min ≥ 10.

**Lemma 11 (two-piece internal bases never block a K4).** Let T = {a₁, …, a₄} with a_k ∈ B_k be pairwise
non-orthogonal, and let B′ ⊆ B_k ∪ B_l be an internal basis. Then T meets B′.

*Proof.* Put A_k = B′ ∩ B_k, A_l = B′ ∩ B_l and W = span A_k, so that span A_l = W^⊥. If T misses B′, then a_k is
orthogonal to A_k, so a_k ∈ W^⊥; and a_l is orthogonal to A_l, so a_l ∈ W. Then a_k ⊥ a_l, a contradiction. ∎

* **Consequence.** V is KS iff every K4 misses some internal basis that meets **at least three** of the B_k. If no
  internal basis meets three covering bases, V is KS iff it is K4-free (Route A), and then d ≥ 17. In particular the
  two-basis links of §3.1 can never help.
* Structure of a blocking basis B′ with pieces A_k (k ∈ S, |S| ≥ 3): the spaces W_k = span A_k are mutually
  orthogonal with ⊕ W_k = ℂ^d, and B_k \ A_k spans ⊕_{l≠k} W_l. B′ blocks exactly the K4s with a_k ∈ B_k \ A_k for all
  k ∈ S, and never a K4 through one of its own rays.
* The lemma is sharp. In d = 3, let e_k ∈ B_k with B_k \ {e_k} spanning e_k^⊥ (k = 1, 2, 3). Then {e₁, e₂, e₃} is a
  3-piece internal basis, and the K4 ((e₂+e₃)/√2, (e₁+e₃)/√2, (e₁+e₂)/√2, a₄) misses it (`ks/reductions_check.py`).

**Lemma 12 (averaging).** Let μ be a probability distribution on K4s and ρ a state with μ(a_k = r) = ⟨r|ρ|r⟩ for every
k and every r ∈ B_k. Then μ-almost every K4 is a KS colouring. In particular V is not KS.

*Proof.* For any basis B′ ⊆ V, |T ∩ B′| ≤ 1 because T is pairwise non-orthogonal, and
E_μ|T ∩ B′| = Σ_{r∈B′} ⟨r|ρ|r⟩ = tr ρ = 1. ∎

* ρ = |w⟩⟨w| with w ∈ B₄ gives the per-ray condition of §3.4. The two-context lemma is the case where Hall's theorem
  supplies μ.
* ρ = I/d gives χ_f(G) > d, the Ramanathan–Horodecki condition for state-independent contextuality (PRL 112, 040404
  (2014)). For four disjoint bases χ_f(G) ≥ 4d/α(G) = d always holds. So the condition says that the 4d rays cannot be
  fractionally partitioned into K4s.
* **Corollary 12a (symmetry kills Route B).** Suppose a group of unitaries or antiunitaries maps each B_k onto itself
  and acts transitively on each B_k. Then **every K4 is a KS colouring**, so V is KS iff it is K4-free, and then
  d ≥ 17. To see this, average one K4 over the group. The marginals are uniform, i.e. ρ = I/d, so every K4 in the
  orbit, T included, meets every basis exactly once.
  * Hence every covariant construction (§3.2–3.3: one group acting transitively on all four bases) is Route A, and
    minimizing the K4 count there was the right objective.
  * A Route-B example, if one exists, has no symmetry group that is transitive on all four bases.
* **Corollary 12b (local form).** Suppose the unitaries that fix w ∈ B₄ and preserve B₁, B₂, B₃ act transitively on
  each N_k(w). Then every triangle of N(w), together with w, is a colouring.

**Theorem 13 (Pauli bases).** Four eigenbases of maximal commuting sets of n-qubit Pauli operators (or n-qudit, p
prime) never form a KS set. Neither do any four bases that are invariant under the full Pauli group.

*Proof.*
1. The Pauli group preserves the eigenbasis of each Lagrangian L_k and acts transitively on it: a Pauli P shifts the
   eigenvalue characters by ω(P, ·). By Corollary 12a, the union is KS only if it is K4-free.
2. Take w in the eigenbasis of L₄. The span S_k of N_k(w) is a joint eigenspace of the abelian group L₄ ∩ L_k ⊆ L₄.
   So the three projectors commute, and tr(Π₁Π₂Π₃) = dim(S₁ ∩ S₂ ∩ S₃) ≥ 1, because w lies in all three.
3. Triangle-freeness gives tr(Π₁Π₂Π₃) = 0 (the commuting-span obstruction, §3.3). Contradiction.
4. Second claim: if the Pauli group preserves B_k, the Paulis that fix a ray of B_k form a subgroup of index ≤ d. They
   share an eigenvector, so they commute and form a Lagrangian, and B_k is its eigenbasis. ∎

* Check (`ks/pauli4.py`):
  * The Mermin star (5 bases, 40 rays, 25 bases inside) is KS, and each of its five 4-subsets is colourable.
  * 4000 random 4-tuples of the 135 three-qubit Lagrangians (half of them chains of overlapping Lagrangians) are all
    colourable, and their K4s always have uniform marginals.
  * Internal bases meeting ≥ 3 covering bases occur in about 30% of these 4-tuples but never block a K4.
* So Mermin-type constructions provably need a fifth context. This extends the abelian covariant argument of §3.2.
* **Consequence for covariant Route A.** Suppose the stabilizer H of a ray w acts on ℂ^d without multiplicities. Then
  the commutant of H is commutative, so all H-invariant subspaces have commuting projectors. The S_k(w) are
  H-invariant, so the contradiction of step 3 applies. A covariant Route-A example therefore needs ray stabilizers
  with multiplicities (as in regular actions, H = 1) and a symmetry group that does not contain the Pauli group.

**Lemma 14 (pairs outside every K4 have overlap ≤ 1/3).** Let a ∈ B_k and w ∈ B_l (k ≠ l) be non-orthogonal rays that
lie in no common K4. Then |⟨a|w⟩| ≤ 1/3.

*Proof.* Let B_m, B_n be the other two bases, X = N_m(a) ∩ N_m(w) and Y = N_n(a) ∩ N_n(w). Since (a, w) lies in no K4,
every ray of X is orthogonal to every ray of Y. So P = P_X, Q = P_Y and R = I − P − Q are orthogonal projectors.
Expanding in B_m gives ⟨a|w⟩ = ⟨a|P|w⟩, and expanding in B_n gives ⟨a|w⟩ = ⟨a|Q|w⟩. Hence ⟨a|R|w⟩ = −⟨a|w⟩. By
Cauchy–Schwarz, |⟨a|w⟩|² is at most each of ‖Pa‖²‖Pw‖², ‖Qa‖²‖Qw‖² and ‖Ra‖²‖Rw‖². The three ‖·a‖² sum to 1, and so do
the three ‖·w‖². By AM–GM the smallest of the three products is at most 1/9. ∎

* **Route A:** every non-orthogonal pair from different bases has |⟨a|w⟩| ≤ 1/3. So every ray is non-orthogonal to
  at least **9 rays of each other basis** (the per-basis version of the bound of 30). By Theorem 5, if it meets
  exactly 9 rays of one basis, it meets at least 12 of each of the other two.
* **Route B:** every cross pair with overlap above 1/3 lies in a K4.
* The same argument, with ψ in place of the fourth basis, applies to any three-context GHZ paradox: every event has
  probability ≤ 1/9, so every context has at least 9 events. The known examples comply: 1/12 (H36), 1/13 (H39), 1/19
  (Perkel), and a maximum of 0.094 for H37. The 9 points of H33 attain 1/9 (§1.7). Lemma 3 of §2 gives a one-line proof
  of the bound.

**Updated picture.**

| | Route A (K4-free) | Route B (K4s exist, all blocked) |
|---|---|---|
| internal bases | irrelevant | blockers must meet ≥ 3 covering bases (Lemma 11) |
| symmetry | allowed; covariant ⇒ KS iff one ray's neighbourhood is triangle-free | no group transitive on all four bases (Cor. 12a) |
| each ray | neighbourhood is a 3-context GHZ graph with ≥ 30 rays, ≥ 9 in each other basis (≥ 12 in the other two if one has exactly 9); all cross overlaps ≤ 1/3 (Lemma 14) | lies in a K4 if d ≤ 17; p(w) ∉ STAB(G_w); pairs with overlap > 1/3 lie in a K4 |
| dimension | d ≥ 18 (Corollary 6; this argument cannot give more than 20) | d ≥ 6 (Theorem 15, computer-assisted) |
| Pauli eigenbases (Lagrangians) | impossible (Thm 13) | impossible (Thm 13) |
| orthogonality graph | may be regular or vertex-transitive; never strongly regular (Thm 16) | never strongly regular, vertex-transitive (Thm 17) or doubly Hoffman-tight regular (Prop. 18) |

Parity proofs are not excluded once internal bases are counted. With disjoint covering bases, each ray lies in
1 + m(r) bases, where m(r) counts internal bases. A parity proof therefore needs every m(r) odd and an odd number M of
internal bases. Double counting, Σ_r (1 + m(r)) = d(4 + M), then forces d even. Every ray lies in an internal basis, so
M ≥ 4, i.e. M ≥ 5.

### 3.8 Path forward

The two routes now live in different regimes. Route B is asymmetric and may sit in small dimension. Route A sits in
d ≥ 17 and allows symmetry, under which KS reduces to a single-ray test. Recommended order:

1. **Route B in small d: exhaustive computer-assisted classification (d = 5 done, d = 6 next).** After §3.7 the problem in fixed d is
   finite combinatorics plus realizability. There are four disjoint d-cliques plus cross orthogonalities, and every K4
   must miss a d-clique that meets at least three of them.
   * `ks/smalld_sat.py` encodes this together with necessary geometric conditions:
     * each ray is orthogonal to ≤ d − 2 rays of another basis;
     * expanding ⟨x|y⟩ in a third basis forces the supports to meet when x ̸⊥ y, and not in exactly one ray when x ⊥ y;
     * no two rays are orthogonal to a common (d − 1)-clique;
     * the span equalities forced by internal bases hold (B_k \ A_k and B′ \ A_k have the same orthogonal rays);
     * no common coordinate subspace exists (Lemma 10);
     * every ray lies in a K4 (d ≤ 16).

     UNSAT would prove that no four-context KS set exists in ℂ^d.
   * **Validation, d = 3 and d = 4:** the same encoding is UNSAT for d = 3 (instantly) and for d = 4 (17 s). This
     reproduces "no four-context KS set in d ≤ 4" without the 18- and 24-ray theorems.
     * Ablation at d = 4: dropping the span, unique-completion or splitting clauses still gives UNSAT.
     * Dropping the support clauses, or "every ray lies in a K4", gives SAT within 1 s.
     * So at d = 4 the trace identity (through "every ray lies in a K4") together with Lemma 11 and the
       support-intersection combinatorics already rules everything out.
   * **d = 5 (Theorem 15 below).**
     * The instance has about 19k variables and 2.06M clauses.
     * Without the unique-completion and span clauses it is satisfiable in 0.4 s, by a model that completes the same
       4-clique twice, which is geometrically impossible. So the geometry carries the difficulty.
     * With all clauses, fix one K4 (WLOG the last ray of each basis) and an internal basis blocking it, made of the
       first rays of the bases it meets. Its piece sizes are 3+1+1, 2+2+1 or 2+1+1+1, the only ≥ 3-piece shapes in
       d = 5 (`--case 311 | 221 | 2111`, with lex-leader breaking of the remaining symmetry).
       * **All three cases are UNSAT** (CaDiCaL 1.5.3: 493 s, 1287 s and 2256 s; `ks/smalld_d5.log`).
     * Checks on the case split: with the geometric clauses removed, each of the three cases is SAT, and brute force
       confirms that the models satisfy "every ray lies in a K4" and "every K4 is blocked by a ≥ 3-piece basis".
     * Independent second encoding (`ks/smalld_cegar.py`):
       * the blockers are encoded with selection variables instead of one variable per d-subset;
       * unique completion and span equalities are added lazily (CEGAR);
       * it is about 10× smaller (228k static clauses at d = 5) and settles d = 4 in 3 s (both shapes UNSAT);
       * it is the intended tool for d = 6, where the first encoding does not fit (the five shapes there are 411,
         321, 222, 3111 and 2211).
     * Soundness: `ks/smalld_soundness.py` fixes the variables to the patterns of genuine configurations, and every
       geometric clause is satisfied. The configurations tested:
       * disjoint bases from the {0, ±1} alphabet, hill-climbed to have 3–5 internal bases meeting ≥ 3 covering bases
         (`--three`), with 47–50 cross orthogonalities;
       * the {0, 1, ω, ω²} alphabet;
       * designed configurations containing a blocking basis of each d = 5 shape.
   * Next steps:
     * CEGAR loop: test each model with `ks/realize.py`. Only the prescribed zeros need fitting, by Lemma 8. Block
       unrealizable patterns by minimal unrealizable sub-patterns, proved exactly with Gröbner bases, as in
       Li–Bright–Ganesh's minimum-KS work.
     * Add stronger geometric constraints: known zero patterns of d × d unitaries, and the support consistency of
       U₁₃ = U₁₂U₂₃.
     * Use the per-ray LP of Lemma 12 as a fast certificate of colourability for candidates.
   * **Theorem 15 (computer-assisted).** No four-context KS set exists in ℂ⁵. Hence every four-context KS set has
     d ≥ 6 (previously d ≥ 5).
     * The proof rests on three things:
       * the validity of each clause group (§3.7 and the list above);
       * the WLOG case split with lex-leader symmetry breaking;
       * the correctness of `ks/smalld_sat.py` and of CaDiCaL.
     * Supporting evidence: the soundness tests, the relaxed-instance checks, and d = 3, 4 reproduced from scratch.
     * Not yet done: an independently checked proof certificate (DRAT/LRAT), and an independent cross-check.
       * The second encoding `ks/smalld_cegar.py` hit its 7000 s limit in all three cases without a verdict. Its last
         logged rounds were 450, 850 and 850, each round slower than the last.
       * It never found a model satisfying all clauses, but it did not reach UNSAT either. So it neither confirms nor
         contradicts Theorem 15.
       * At d = 5, adding completion and span clauses lazily is slower than adding them all at once. So for d = 6 a
         hybrid is the better design: completion clauses static, span clauses lazy.
   * Next: d = 6 (24 rays, the dimension of the 21-ray set), using `ks/smalld_cegar.py` and the five shapes 411, 321,
     222, 3111 and 2211. Outcome: either a 24-ray four-context KS set, or d ≥ 7.
2. **Route A with rigid symmetry.** By Corollary 12a, in a covariant configuration KS holds iff the neighbourhood of
   **one** ray is triangle-free. Non-regular monomial constructions are discrete and can be enumerated.
   * Data: a finite group G, a representation ρ of degree d ≥ 17, and four inequivalent monomial structures
     ρ ≅ Ind_{H_k}^G χ_k (H_k of index d, χ_k linear). b_k spans the χ_k-eigenline of H_k, and B_k = ρ(G)·b_k.
   * Overlaps come from Mackey theory: ⟨b_k|ρ(g)b_l⟩ ≠ 0 only if χ_k = χ_l^g on H_k ∩ gH_lg⁻¹.
   * Filters: ρ restricted to H₄ has multiplicities; G does not contain the Pauli group; each ray has ≥ 28
     non-orthogonal rays.
   * Tool: GAP (SmallGroups; irreducible characters of degree 17–64 and their monomial structures).
   * The regular case H_k = 1 was searched in §3.3. The non-regular case is the one the Mermin star uses
     (|G/Z| = 64, d = 8, |H_k/Z| = 8). The Pauli group itself is excluded by Theorem 13, but no other group with
     non-regular monomial structures has been searched.
3. **Theory.**
   * Route B: show that internal bases blocking every K4 force a state ρ and a K4 mixture as in Lemma 12. A first
     test is to run the LP of Lemma 12 (`mixture_with_marginals` in `ks/fourbasis.py`) on random and structured
     configurations that have ≥ 3-piece internal bases, for ρ = I/d and ρ = |w⟩⟨w|.
   * Route A: find a global form of the trace identity. Σ_{a∈B₁, b∈B₂, c∈B₃} ⟨a|b⟩⟨b|c⟩⟨c|a⟩ = tr I = d, while
     K4-freeness forces Σ_{w∈B₄} |⟨w|a⟩⟨w|b⟩⟨w|c⟩|² = 0 for every pairwise non-orthogonal (a, b, c). Combining such
     identities over all four bases might raise d ≥ 17 or exclude Route A altogether.

### 3.9 Combinatorial skeleton of Route A (merged from branch `claude/ghz-contextuality-graphs-1ur6ur`)

This section was written in parallel with §3.7–3.8. Overlaps: Lemma A below is Lemma 11 of §3.7. Lemma B (≥ 3 events
per GHZ context) is superseded by Lemma 14 (≥ 9 events per context, every event probability ≤ 1/9). Open direction 2
(Route B in d = 5, 6) is done for d = 5 by Theorem 15 (§3.8).

**Lemma A (two-piece internal bases are redundant).** Let B′ = A₁ ∪ A₂ (A₁ ⊂ B₁, A₂ ⊂ B₂) be an internal basis.
For a transversal with ⟨a|b⟩ ≠ 0 (a ∈ B₁, b ∈ B₂), "exactly one of a ∈ A₁, b ∈ A₂" holds automatically: both is
impossible since A₁ ⊥ A₂; neither gives a ∈ span(B₁ \ A₁) = span(A₂) and b ∈ span(A₁), so a ⊥ b. Hence Route B needs
internal bases with **≥ 3 pieces**: a coordinate subspace of one basis spanned jointly by rays of two others. This is
why the linked completions of H36 in §3.1 stayed colourable.

**Lemma B (no GHZ class of size ≤ 2).** In a three-context GHZ paradox every context has ≥ 3 events. If N₃ = {c, c′},
rigidity (every event has neighbours in both other classes, which holds whenever the contexts are remnants of complete
bases) splits N₁ = A ⊔ A′ with A = N(c), A′ = N(c′), and N₂ = B ⊔ B′ likewise; triangle-freeness forbids A–B and
A′–B′ edges. Then {c, A′, B′} ∪ {c′, A, B} is a proper 2-colouring of H, so θ(H̄) ≤ χ(H) = 2 < 3. Size 3 admits odd
cycles (c₁–a₁–b₂₃–c₃–a₁₃), so the bound is exactly 3.

**Combinatorial skeleton of Route A** (`ks/skeleton_sat.py`, checker `ks/skeleton_check.py`). Necessary conditions
on the 4-partite non-orthogonality graph of four complete bases, all proved from completeness alone:
* (K4) no transversal 4-clique;
* (C2) a ~ b ⇒ a, b have a common neighbour in each other basis (⟨a|b⟩ = Σ_c ⟨a|c⟩⟨c|b⟩);
* (C4) a ⊥ b ⇒ a, b have 0 or ≥ 2 common neighbours in each other basis (the zero-pattern condition of unitaries);
* (DEG) ≥ 3 neighbours in each other basis (Lemma B) and ≥ 28 in total (§2);
* (TRI) each triple of bases has a triangle (Re tr U₁₂U₂₃U₃₁ = d).
With sound symmetry breaking (double-lex on the B₁–B₂ block, column-lex on B₂–B₃ and B₃–B₄):
* d = 4, 5: **UNSAT** (instant). d = 6: **SAT**, a perfectly regular skeleton (4 neighbours per other basis).
* With the full requirement n ≥ 28: d = 10 **UNSAT**; d = 17 **undecided** (CaDiCaL, 14000 s, no verdict). Deciding
  it needs a stronger encoding (e.g. lazy θ = 3 neighbourhood constraints, direction 1 below) rather than more time.
  *Later:* d = 17 is now excluded analytically (Corollary 6 gives d ≥ 18), and (DEG) can be raised to ≥ 30 in total,
  with ≥ 12 in the other two bases whenever one basis contributes exactly 9 (Theorem 5).
* So the combinatorial conditions alone cannot give impossibility. Since d ≥ 18 analytically, the content is in
  θ(H̄_w) = 3 for every neighbourhood. UNSAT at a given d ≥ 18 would still be a theorem for that d.

**Open directions.**
1. *CEGAR skeleton search (finding and per-d impossibility).* SAT skeleton → SDP check θ(H̄_w) = 3 for every w → block
   failing neighbourhoods with the sound clause "some non-edge inside N(w) becomes an edge or N(w) gains a vertex"
   (θ(H̄) is monotone under adding edges to H) → repeat; realizability of the final zero pattern by unitaries last.
2. *Route B in small d by case analysis.* Disjointness (§3.4), Lemma A (≥ 3-piece internal bases), XCG's bound (every
   ray needs ≥ 10 neighbours since its reduction is a GHZ-type proof), and |V| = 4d make d = 5, 6 finite case studies.
3. *Block-generic models for Route B.* Enumerate subspace lattices generated by 4 bases (which coordinate subspaces of
   one basis are sums of coordinate subspaces of others); orthogonality is then decided by the lattice, generic entries
   realize it automatically, and KS-colourability is a SAT question. Lemma A says only ≥ 3-piece lattice elements matter.
4. *Global analytic coupling.* The per-w identity tr((I−Q₁)(I−Q₂)(I−Q₃)) = 0 with Q_k the projector onto the rays
   of B_k orthogonal to w, summed over w ∈ B₄, ties the orthogonality counts z_k(w) to the T-masses between
   orthogonal-to-w sets; combined with the uniform marginals Σ_w |⟨r|w⟩|² = 1 it is the only identity we have that
   sees K4-freeness.
5. *Sparse-unitary theory.* Route A is a zero-pattern problem for U₁₂, U₂₃, U₃₄ with U₁₂U₂₃ = U₁₃ etc.; results on
   admissible zero patterns of unitary matrices may give conditions beyond (C4).

### 3.10 Can the orthogonality graph be strongly regular? Regular? (new)

Notation: G is the orthogonality graph of V = B₁ ⊔ … ⊔ B₄ (n = 4d), and H = Ḡ is the non-orthogonality graph
(4-partite, with the bases as parts). For every such V, ϑ(G) = 4 and ϑ(H) = d (§3.6). The bases are maximum cocliques
of H (α(H) = d), and ω(H) ≤ 4, with ω(H) = 4 exactly in Route B. For comparison:
* the paper's Theorem 3 says no strongly regular graph (SRG) hosts a three-context GHZ paradox;
* the four-context GHZ (Mermin) paradox *is* hosted by an SRG: G is the complement of the Shrikhande graph, (16, 9, 4, 6).

For KS sets the answer is negative.

**Theorem 16 (no strongly regular graph).** For four pairwise disjoint orthonormal bases whose orthogonality graph is
strongly regular, the union is KS-colourable. So the orthogonality graph of a four-context KS set is never an SRG
(equivalently, neither is H).

*Proof.* G is primitive: if G = 4K_d there are no cross orthogonalities, and if G is complete multipartite then H is a
union of cliques of size ≤ 4, whereas every ray has ≥ 2 non-orthogonal rays in each other basis. Let k > r > s be
its eigenvalues.
1. *Every K4 is a colouring.*
   * For an SRG, ϑ(G) equals the Hoffman bound −ns/(k − s), and ϑ(H) = n/ϑ(G) = 1 − k/s, the Delsarte clique bound
     of G. (The feasible point aI + bĀ and the dual point J + tA meet; `ks/srg_check.py` has the parameter checks.)
   * So −ns/(k − s) = 4 and 1 − k/s = d.
   * Every 4-coclique T of G (a K4 of H) attains Hoffman's bound, so 1_T − 𝟙/d is an s-eigenvector.
   * Every basis C ⊆ V, covering or internal, is a d-clique attaining Delsarte's bound, so 1_C − 𝟙/4 is an
     r-eigenvector.
   * These eigenspaces are orthogonal, so |C ∩ T| = 1. A K4 that meets every basis once is a KS colouring, so a KS
     set would have to be K4-free (Route A).
2. *Parameters.* A Delsarte clique has the property that every vertex outside it has the same number of neighbours
   in it. So every ray outside B_j is orthogonal to the same number e of its rays, with 1 ≤ e ≤ d − 2 by
   disjointness.
   * Hence k = d − 1 + 3e, r = d − 1 − e and s = −k/(d − 1) = −1 − 3e/(d − 1).
   * Since n = 4d ≢ 1 (mod 4), G is not a conference graph, so s is an integer: s ∈ {−2, −3}.
3. *s = −2.* Then d = 3e + 1. Integral multiplicities force e ∈ {1, 2, 5}, i.e. d ∈ {4, 7, 16}, but Route A needs
   d ≥ 17.
4. *s = −3.* Then H has the parameters (12t + 4, 3t + 3, 2, t + 1) of a GQ(3, t) collinearity graph, with d = 3t + 1.
   * Since λ_H = 2, the graph H[N(w)] is 2-regular for every ray w, i.e. a disjoint union of cycles.
   * The three remnants sum to 3 in state w, so ϑ(G_w) ≥ 3 for the orthogonality graph G_w on N(w).
   * ϑ of the complement of a union of cycles is the largest ϑ(C̄_m) over its components. That is 3 for m = 3, 2 for
     even m ≥ 4, and 1 + 1/cos(π/m) < 3 for odd m ≥ 5.
   * So H[N(w)] contains a triangle, i.e. w lies in a K4, contradicting step 1. ∎

* The Mermin GHZ graph is the t = 1 member of the s = −3 family: H is the Shrikhande graph, whose local graphs are
  6-cycles (ϑ = 2 < 3). This is why it cannot be completed to four bases of ℂ⁴ even though it passes all the global
  invariants of §3.6.
* The other (16, 6, 2, 2) graph is the rook graph L₂(4) with a Latin-square partition. There each of its 8 K4s meets
  each of its 24 candidate bases exactly once, as step 1 predicts (`ks/srg_check.log`).

**Theorem 17 (a vertex-transitive graph forces Route A).** If the orthogonality graph G is vertex-transitive, every K4
is a KS colouring. Here any graph automorphisms count, not only those induced by unitaries. So a four-context KS set
with vertex-transitive G is K4-free, with d ≥ 17 and all of Route A.

*Proof.* Average |gT ∩ C| over g ∈ Aut(G), for a K4 T and a basis C. The average is |T||C|/n = 1, and each term is
≤ 1 because gT is a coclique and C a clique of G. So every term is 1, including g = 1. ∎

* This is the combinatorial form of Corollary 12a. It covers Cayley graphs and all symmetric constructions.
* A Route-B example therefore needs an orthogonality graph that is not vertex-transitive and not doubly Hoffman-tight
  (Proposition 18, part 3).

**Proposition 18 (regular graphs are allowed, with spectral conditions).** Regularity is not excluded: covariant
constructions give regular graphs, and so do the four product bases of the Mermin star. If H is D-regular, i.e. every
ray is non-orthogonal to D rays:
1. λ_min(H) ≤ −D/3, from ϑ(H) = d ≤ −n λ_min/(D − λ_min). Equality holds iff the bases are Hoffman cocliques of H, and
   then every ray meets exactly D/3 rays of each other basis.
2. λ_min(G) ≤ −k_G/(d − 1), equivalently λ₂(H) ≥ (3d − D)/(d − 1), from ϑ(G) = 4 ≤ Hoffman's bound.
3. If both are equalities, every K4 meets every basis exactly once. The argument is step 1 of Theorem 16, with the
   eigenspaces E_{λmin(G)}(G) and E_{−1−λmin(H)}(G). So Route B is impossible there.
4. *Flat* configurations are those in which every nonzero overlap between different bases has modulus 1/√m, i.e. the
   transition matrices are scaled weighing matrices.
   * They are D-regular with D = 3m and λ_min(H) = −m, because the Gram matrix of the projectors, I + A_H/m, is PSD
     and singular.
   * In Route A they need m ≥ 10 (Lemma 14 and D ≥ 28).
   * Every non-orthogonal pair then has between 2√m and m − √m common non-orthogonal rays, by the three-way split in
     the proof of Lemma 14.

* Check (`ks/srg_check.py`): the four Mermin product bases (d = 8) are flat and 12-regular. They have
  λ_min(H) = −4 = −D/3 (equality in 1) and λ₂(H) = 4 > 1.71 (strict in 2), and they are colourable (Theorem 13).
* In Route A every local graph H[N(w)] is a three-context GHZ graph.
  * By the paper's Theorem 3 it is not strongly regular.
  * Since ϑ(complement) ≤ χ_f, it cannot be triangle-free with maximum degree ≤ 3, because such graphs have
    χ_f ≤ 14/5 (Dvořák–Sereni–Volec).
  * So in Route A every ray w has a non-orthogonal neighbour a with at least 4 common non-orthogonal rays, and H
    cannot be edge-regular with λ_H ≤ 3.
  * With Theorem 17, a vertex-transitive edge-regular H with λ_H ≤ 3 is impossible outright.
* **Consequence for searches.**
  * Route A searches may stay inside regular or vertex-transitive graphs. These must satisfy the inequalities of
    Proposition 18, d ≥ 17, D ≥ 28 with ≥ 9 per basis, and GHZ local graphs.
  * Route B searches must leave the vertex-transitive and strongly regular worlds altogether.

### 3.11 Building a four-context KS set from H33 (attempt; new)

Template (§3.1): put ψ in B₄. Complete the contexts of H33 to bases: B₁ = points + A′, B₂ = B-lines + B′,
B₃ = C-lines + C′. Then ψ is forced to 0, and every other ray of B₄ must also be forced to 0. Everything below is
reproduced by `ks/h33_attempt.py` (log `ks/h33_attempt.log`, about 30 s).

**Theorem 19 (H33 cannot be used in Route A).** Let B₁, …, B₄ be orthonormal bases of ℂ^D (any D). Suppose some ray w
has a neighbourhood whose non-orthogonality graph is H33. Then the union contains K4s, so the K4-free route is closed.

*Proof.*
1. By Theorem 5 the realization of H33 is unique. In particular, each B-line and each distinct C-line parallel to it
   have overlap exactly 1/2.
2. 1/2 > 1/3, so by Lemma 14 each of these 24 pairs lies in a K4. ∎

So a KS set built on H33 must use Route B: internal bases meeting ≥ 3 covering bases have to block those K4s. Among
the known GHZ graphs, H33 is the only one whose forced realization breaks the Route-A bound:

| graph | cross overlaps of its known realization |
|---|---|
| H36 | 1/4 |
| H39 | at most 1/3 |
| Perkel | exactly 1/3 (Gram matrix I + A/3, rank 37) |
| H37 | 0.51 in the certified realization, but H37's realizations are not unique, so no obstruction follows |

**Lemma 20 (state-averaging test).** Let V be the union of four orthonormal bases, with any internal bases. If the Born
vector `(⟨r|ρ|r⟩)_r` of some state ρ is a convex combination of K4 indicators, then V is not KS.

*Proof.* Each basis of V has Born weight 1. Each K4 meets each basis at most once. So every K4 in the support of the
combination meets every basis exactly once, which makes it a KS colouring. ∎

Moreover, no state passes this test iff there are weights λ_r with `Σ_{r∈T} λ_r ≥ 0` for every K4 T and
`Σ_r λ_r |r⟩⟨r| ≺ 0`. This follows by separating the two convex sets, and shifting λ by a constant on one basis.

* ρ = I/d is Lemma 12.
* ρ = |y⟩⟨y| is the per-ray condition of §3.4. In the form used below, NCF(y) is the largest weight of K4 mixtures
  through y; NCF(y) = 1 means y can carry the 1 of B₄.

**Structure of the H33 realization (exact, `ks/h33model.py`).**
* **Class subspaces.** ℂ²⁰ = V₁ ⊕ V₂ ⊕ V₃ ⊕ V₄, one 5-dimensional subspace per parallel class.
  * V_i is spanned by the class sum u_i and the two Fourier blocks of that class.
  * The lines of class i lie in V_i. Points and ψ spread over all four V_i with weight 1/4 each.
* **Natural completions.** In each Fourier block the six relevant directions are 30° apart. The block-pure completions
  B′ = {β_k} and C′ = {γ_k} make B₂ and B₃ **block-diagonal**: 5 rays in each V_i.
* **Ten equiangular lines.** Inside V_i, the components of ψ, of the points and of everything below run over 10
  equiangular lines with overlap 1/3. All their triple products are real (±1/27), so this is the maximal system of
  10 equiangular lines in ℝ⁵.
* **Duality W (exact).** There is a unitary involution W with W ψ = the point (0,0). W maps B₂ and B₃ onto themselves,
  fixing the 4 lines through (0,0) and exchanging the other 8 with the β's (resp. γ's). Consequently:
  * Each point a is itself a GHZ state for (B₂, B₃, B₄) whenever B₄ contains ψ and 8 specific vectors (its "dual
    points").
  * The neighbourhood of a is then an exact copy of H33. Its lines are the 4 lines through a plus the 8 β's (resp. γ's).
* **Orbit.** Γ = ⟨AGL(2,3), W⟩ has an orbit of **10000 = 10⁴ rays** through ψ (`ks/h33orbit.py`).
  * The orbit is a product of the four 10-line systems. All overlaps are multiples of 1/6.
  * Every orbit ray carries its own copy of H33: 9 orbit rays, 12 rays of B₂ and 12 of B₃. It is a GHZ state as soon
    as B₁ contains its 9 "points".
  * No orbit ray is orthogonal to a full 9-context. The 72 dual points of the 9 points are distinct, and no two of the
    dual 8-sets are mutually orthogonal.

**Attempts in d = 20. All of them fail.**

| design | GHZ rays | K4s | internal bases (≥ 3 pieces) | outcome |
|---|---|---|---|---|
| random completions, random B₄ | ψ | 80408 | – | every other ray of B₄ has NCF = 1; the uniform state on ψ⊥ is a K4 mixture |
| symmetric completion, random B₄ | ψ | 7448 | – | uniform state on ψ⊥ is a K4 mixture |
| block-aligned (B₄ = {e, β} in every block) | ψ | 512 | 111 (65) | all 512 K4s meet every internal basis once |
| W-symmetric, B₄ = W(B₁) | ψ, (0,0) | 1184 | 14 (0) | all K4s are colourings; uniform state on {ψ, a₀}⊥ is a K4 mixture |
| second copy of H33 with an orthogonal state in the same three bases (`ks/h33copies.py`) | – | – | – | not found (best penalty 2.7, 30 starts per placement) |

* Contextual states orthogonal to ψ do exist: NCF ≈ 0.80 for the symmetric completion and ≈ 0.96 for random ones
  (`ks/h33ctx.py`, smoothed search, checked on ψ itself). So the per-ray test alone does not rule out d = 20.
* The failure is structural. With block-diagonal B₂, B₃, a K4 must take its B₂- and B₃-rays from the same V_i. The
  only internal bases inside B₂ ∪ B₃ are the 14 "block switches", and those meet only two covering bases.
* What a construction in d = 20 still needs:
  * completions that couple the four class subspaces;
  * exact coincidences that create internal bases meeting 3–4 covering bases;
  * failure of the state-averaging test (Lemma 20) for every state.
* This is not an impossibility proof. Completions in special position are not exhausted.
* Next steps:
  * a universal state ρ for Lemma 20, which would close d = 20;
  * larger D, reusing the GHZ web of the orbit;
  * for Route A, only graphs whose cross overlaps are at most 1/3 (H36, H39, Perkel) remain usable.

### 3.12 Families of K4-free four-colourable graphs (overnight search; new)

Proposal: construct families of K4-free four-colourable graphs H, 0-1 colour the complement, and look for KS sets.

**The colouring step is automatic.** Take the colour classes to be the four bases. A 0-1 colouring of the
complement (the orthogonality graph) needs one ray per basis, pairwise non-orthogonal, which is a K4 in H. So
every K4-free H gives a non-colourable complement, and the SAT check below confirms this on every candidate.
* Complete classes lose nothing (Corollary 9).
* The whole question is therefore **realizability** (Route A). It needs d ≥ 18, and the neighbourhood of every
  ray must be a three-context GHZ graph with ϑ = 3.

**Families.** Four-partite Cayley graphs: Γ a group of order d, (k, g) ~ (l, h) iff g⁻¹h ∈ S_kl.
* These are exactly the non-orthogonality graphs of four bases on which Γ acts regularly.
* By Corollary 12a any four-context KS set with such a symmetry is K4-free. So "UNSAT for Γ" means that **no
  four-context KS set of either route has a symmetry group acting regularly on all four bases with that Γ**.
* Two families are searched:
  * **multiplier family:** Γ = Z_p and every S_kl a union of orbits of the e-th powers in Z_p^*;
  * **group family:** all Cayley graphs over each group of order 18–40 in `ks/groups.py` (every abelian
    group, dihedral, dicyclic and metacyclic groups, a few products; complete for orders 18–23).

**Method** (`ks/cayley4.py`): SAT over the connection sets with every necessary condition we have:
* (K4) no transversal K4;
* (DEG) ≥ 9 neighbours per basis (Lemma 14), ≥ 30 in total (Corollary 6), and the Theorem 5 refinement;
* (C2) and (C4), including (C4) for pairs inside one basis;
* (TRI) every three bases contain a triangle;
* normalisation (group family only, see the correction below), and lex-leader symmetry breaking.

New sound pruning — **quotient constraints**. A homomorphism Γ → Z_m maps every neighbourhood homomorphically into
the neighbourhood of a residue pattern. If that pattern neighbourhood is bipartite, ϑ = 2 < 3, so such patterns are
forbidden for m = 2, 3. These were the dominant failure: without them, every model over Z₁₈ and Z₃×Z₆ had all four
neighbourhoods bipartite, because each S_kl avoided a whole coset mod 3.

CEGAR on the local GHZ condition: an SDP computes ϑ(N̄(v)) for one vertex per basis. A failure is blocked by "N(v)
gains a vertex or an edge"; a bipartite N(v) is blocked by "an edge inside one side". Candidates passing all four
checks get:
* a brute-force recheck of every condition;
* the 0-1 colouring of the complement (SAT);
* a Douglas–Rachford realizability test on the 4d × 4d Gram matrix, which must be 4 × a rank-d projector with the
  prescribed zeros.

**Validation.**
* The SAT encoding agrees with an independent brute-force checker:
  * every model passes the checker;
  * on random assignments over D₁₈ and Z₃⋊Z₆ the solver and checker give the same verdict (40/40);
  * the orbit version agrees on 30/30 random invariant assignments over Z₃₁;
  * it reproduces an exhaustive enumeration of the cyclotomic family (`ks/cyclo4.py`, without the normalisation: no
    survivor for e = 2 and any prime 19 ≤ p ≤ 151, `ks/cayley4_mult/cyclo4_e2.log`; an earlier check with the
    normalisation covered e = 3, p ≤ 73).
* Realizability controls:
  * four bases of ℂ²⁰ with a planted zero pattern (31% zeros) are recovered to residual 2·10⁻¹²;
  * an impossible pattern stays at residual 1.18.

**Correction (29 Oct).** The pilot and the overnight run imposed the normalisation e ∈ S_3k on the multiplier
family too. There it is not without loss of generality: M fixes 0 in every basis, and a translation of one basis does
not commute with x ↦ mx. Their multiplier-family verdicts therefore cover only the subfamily with (3, 0) ~ (k, 0) for
k = 0, 1, 2. The encoder (and `ks/cyclo4.py`) now drop the normalisation for M ≠ 1 and use all four bases in the
lex-leader clauses instead. The 77 UNSAT instances were re-run: 56 stay UNSAT (below). The group family is unaffected:
there the normalisation is a right translation, which commutes with Γ.

**Daytime pilot (2 min per instance).**
* Multiplier family (normalised subfamily only, see the correction): UNSAT for p = 19 (every e), p = 23 (e = 2),
  p = 29 (e ≤ 7), p = 31 (e ≤ 6), p = 37 (e ≤ 6) and p = 67, 71 (e = 2, 3, 6 and e = 2).
  * Several of these became UNSAT only after the CEGAR had rejected 2–36 graphs whose neighbourhoods have
    ϑ ≈ 2.2–2.6 < 3.
  * Small orbits (size 2–3, close to the full group family) need longer.
  * No graph passed the GHZ check.
* Group family, order 18: models exist, but every one found had a bipartite or ϑ < 3 neighbourhood. Not yet decided.

**Overnight run (28–29 Oct).** `ks/cayley4_run.py` on 4 cores, resumable (state in `ks/cayley4_runs/`, log
`ks/cayley4.log`, table `ks/cayley4_results.jsonl`).
* Job order:
  * multiplier instances with orbit size h = (p − 1)/e ≥ 3, primes 19 ≤ p ≤ 151, 900 s each;
  * the group family for orders 18–24, 1800 s per group;
  * the multiplier instances with h = 2.
* It ran 20:15–01:25 UTC (19.8 CPU-hours). Two container restarts cut it short. The first cost one chunk's in-flight
  jobs. The second stopped the run at 01:25, the 01:30 safety check did not relaunch it, and the last 4.5 hours
  were lost.
* 151 of the 204 jobs ran: 3217 models examined, 11 802 ϑ-blocks and 939 bipartite blocks learned.
* **No graph passed the local GHZ test** (ϑ = 3 at all four bases), so the colouring and realizability stages were
  never reached.
* Multiplier family, normalised subfamily: 77 of the 145 instances with h ≥ 3 were UNSAT and 67 were undecided at
  900 s. The instance (151, 10) and the 29 instances with h = 2 were not reached.
* Group family: nothing decided.
  * Z₁₈ and Z₃×Z₆: no model at all in 1800 s each. With the quotient constraints the solver no longer finds the
    bipartite blow-ups mod 3, and finds nothing else.
  * D₁₈ (32 models, 1806 s), Dih(Z₃×Z₃) (149 models, 1636 s), Z₃⋊Z₆ (23 models, 395 s): in every model all four
    neighbourhoods were bipartite (ϑ = 2), although none is forced to be by a quotient pattern (those are excluded).
  * Z₁₉ (375 s) and Z₂₀ (334 s): no model yet. The other 23 groups of orders 20–24 were not reached.

**Corrected multiplier run (29 Oct).** The 77 instances above without the normalisation, 600 s each (4.1 CPU-hours;
`ks/cayley4_mult/`, per-instance table in `summary.txt`).
* **UNSAT: 56 instances.**
  * e = 2 and e = 3: every prime 19 ≤ p ≤ 151 (29 + 14 instances);
  * e = 4: p = 29, 37, 41, 101, 109, 149;
  * e = 5: p = 31, 41; e = 6: p = 19, 67, 79, 97; e = 7: p = 29.
* Undecided at 600 s: 21 instances.
  * e = 4: p = 53, 61, 73, 89, 97, 113, 137; e = 5: p = 61, 131, 151;
  * e = 6: p = 31, 37, 43, 73, 109, 139, 151; e = 7: p = 43, 71; e = 8: p = 97, 137.
* 15 of the 56 needed the ϑ-CEGAR (up to 107 rejected graphs, for p = 41, e = 4).
* The largest blocked ϑ was 2.966 (p = 89, e = 4), and 2.922 within the UNSAT instances. The 10⁻⁶ threshold is far
  away, so the blocks are not a solver-precision artefact.
* Meaning of UNSAT: for these (p, e) no four-context KS set has a non-orthogonality graph on which the affine group
  {x ↦ mx + t : m ∈ H_e, t ∈ Z_p} acts on every basis (H_e = the e-th powers in Z_p^*).
* For e = 2 the obstruction is pure combinatorics of quadratic residues. In the brute-force check every member fails
  (DEG) or contains a transversal K4. For 23 ≤ p ≤ 151 the counts do not depend on p: of 59 189 classes, 35 829 fail
  (DEG) and 23 360 contain a K4.
  * This should hold for every p. It would follow from the 3-e.c. property of Paley graphs (p ≡ 1 mod 4) and of
    quadratic-residue tournaments (p ≡ 3 mod 4), plus a finite check. Candidate theorem, not yet written out.

**A candidate: p = 137, e = 4.** The undecided instance (137, 4) produced the first graph, in either family, that
passes every test we have.
* Connection sets: S₀₁ = S₂₃ = H₄ (the fourth powers) and S₀₂ = S₀₃ = S₁₂ = S₁₃ = 3H₄. Each has 34 elements and is
  symmetric (−1 ∈ H₄ since p ≡ 1 mod 8).
* Checks (`ks/cayley4_candidate_p137.json`):
  * (K4), (DEG), (C2), (C4) and (TRI) hold (brute-force recheck);
  * every neighbourhood is a 102-vertex, 12-regular, triangle-free, non-bipartite three-partite graph with
    ϑ = 3.000000000000 (Clarabel at tolerance 10⁻¹²). By the criterion in the notation section these are
    three-context GHZ graphs (new ones, tri-Cayley over H₄ ≅ Z₃₄);
  * the complement is not 0-1 colourable, as it must be;
  * global invariants (§3.6, §3.10): ϑ(G) = 4 and ϑ(H) = 137 = d; the graph is vertex-transitive (D₄ permutes the
    bases);
  * Proposition 18: H is 102-regular with λ_min(H) = −34 = −D/3 (equality, multiplicity 3) and
    λ₂(H) = 12.76 ≥ 2.27. This is the spectral signature of a *flat* configuration with m = 34, i.e. transition
    matrices that are complex weighing matrices W(137, 34);
  * every non-orthogonal pair has exactly 6 + 6 common non-orthogonal rays, just above the flat bound √34 = 5.83 per
    basis.
* Realizability in ℂ¹³⁷ (548 rays) is open:
  * Douglas–Rachford on the 548 × 548 Gram matrix (`cayley4.realize` projections) does not converge. From four
    random starts (two run to 3000 iterations, two to 900) it plateaus at residual ≈ 15.7. With the overlaps forced
    flat (modulus 1/√34 on the pattern), two runs of 2000 iterations plateau at ≈ 27.3. A d = 137 control with a
    known solution and this sparsity is not available, so this is evidence, not proof.
  * A realization on which Z₁₃₇ acts by unitaries would need a function supported on a coset of H₄ with a flat
    Fourier transform. Alternating projections from 200 starts stay at rms deviation 0.11 from flat, so numerically
    there is no such symmetric realization.
  * Parameter count: with basis 3 fixed, three unitaries (3p² = 56 307 real parameters) must satisfy
    12p(p − 34) = 169 332 real zero conditions. Only special algebraic structure could beat that.
  * A rigid feature to exploit: (0, x) and (1, x) are orthogonal but have the same 34 + 34 neighbours in bases 2
    and 3. So for every x the spans of {(2, x + 3H₄)} and {(3, x + 3H₄)}, both 34-dimensional, meet in at least 2
    dimensions; likewise for bases 0, 1 with (2, x), (3, x).
* So the candidate is most likely not realizable, but this is not proved.

**The candidate is one member of a family.** I scanned every pattern in which each S_kl is a single class of index 4
(4⁶ patterns, up to multiplying by a generator) for the primes p ≡ 1 mod 4 with 41 ≤ p ≤ 229, through (DEG), (K4),
(C2), (C4), (TRI) and ϑ = 3 at all four bases (`ks/cayley4_mult/siblings_e4.py`, log `siblings_e4.log`). Survivors,
up to basis permutations and multipliers:

| p | rays | pattern (classes of S₀₁, S₀₂, S₀₃, S₁₂, S₁₃, S₂₃) | neighbourhood N(v) | common nbrs of a non-orth. pair | √h | Douglas–Rachford |
|---|---|---|---|---|---|---|
| 97 | 388 | 0, 1, 2, 2, 3, 0 | 72 vertices, 10-regular | 5 + 5 | 4.90 | no convergence (≈ 13.3; flat ≈ 22.9) |
| 137 | 548 | 0, 0, 3, 3, 0, 0 (a perfect matching gets one class) | 102, 12-regular | 6 + 6 | 5.83 | no convergence (≈ 15.7; flat ≈ 27.3) |
| 157 | 628 | 0, 1, 2, 0, 3, 1 | 117, 14-regular | 7 + 7 | 6.24 | not run |
| 173 | 692 | 0, 0, 0, 1, 3, 1 | 129, degrees 15–18 | 9 + 9 or 9 + 6 | 6.56 | not run; **no flat realization** (6 < √43) |
| 193 | 772 | 0, 0, 0, 2, 3, 3 | 144, degrees 16–19 | 11 + 8 or 8 + 8 | 6.93 | not run |
| 193 | 772 | 0, 0, 3, 3, 0, 0 (same shape as 137) | 144, 16-regular | 8 + 8 | 6.93 | not run |

* Class i is g^i H₄ for the smallest primitive root g (`cyclo4.classes`); h = (p − 1)/4.
* There are no survivors at p = 41, 53, 61, 73, 89, 101, 109, 113, 149, 181, 197, 229: most patterns contain a K4,
  and the rest fail ϑ = 3. p = 233 and 241 were not scanned.
* Every member was re-verified independently (`ks/cayley4_candidate_p*.json`):
  * the conditions hold by brute force, and the complement is not 0-1 colourable;
  * all four neighbourhoods have ϑ = 3 to 10⁻¹¹, so they are three-context GHZ graphs;
  * H is D-regular with D = 3h and λ_min(H) = −h = −D/3, the flat signature of Proposition 18, and λ₂(H) is well
    above its bound;
  * so a(I + A_H/h) is a feasible point with value 4, and ϑ(G) = 4 for every member.
* The last two columns test flat realizations, where every nonzero overlap has modulus 1/√h. A non-orthogonal pair
  with t common neighbours in some basis needs t ≥ √h. At p = 97, 137 and 157 this holds only just: t unit phases
  must sum to modulus √h, just below t.

Taken together, these graphs pass every necessary condition in this file, local and global, so those conditions do
not characterise four-context KS sets. Either some of these graphs are realizable (they would be the first
four-context KS sets, the smallest in d = 97), or there is an obstruction we have not yet found. The near-tightness
of the flat bound suggests looking first at flat realizations: complex weighing matrices W(p, h) with cyclotomic
support.

*Update, 30 Oct (§3.13):* the scan is now complete for every p. The all-fourth-powers pattern adds four members at
p = 281, 337, 353 and 457. The Douglas–Rachford column says nothing: DR also stalls on a planted realizable control.
Calibrated searches, and a proof that p = 97 has no flat realization, are in §3.13.

**What the runs show.**
* Every other graph the solver produced, in both families, failed at a single ray: its neighbourhood was bipartite or
  had ϑ < 3. The local GHZ condition is what kills these symmetric families.
* The decided instances are the ones with few orbits. For small orbits and for the group family the SAT search itself
  is the bottleneck (no model in 30 min for Z₁₈ or Z₃×Z₆).
* All of this is computer-assisted without certificates: no DRAT proofs yet, and the ϑ-blocks rest on floating-point
  SDP values (with margin ≥ 0.03).

**Next steps.** Items 1–2 of the earlier list are done or superseded in §3.13. Still open:
3. Resume the undecided instances, the h = 2 instances and the group family. The driver resumes from
   `ks/cayley4_runs/` or `ks/cayley4_mult/`, and `--only` selects jobs.
4. Certificates before anything goes into a paper: DRAT/LRAT proofs for the SAT part, rational dual certificates
   (ϑ ≤ 3 − δ) for each ϑ-block.

### 3.13 Multiplier families: no quadratic or cubic examples, a finite quartic list (new, 30 Oct)

F(p, e) is the multiplier family of §3.12:
* four-partite Cayley skeletons over Z_p, (k, x) ~ (l, y) iff y − x ∈ S_kl, with S_lk = −S_kl;
* every S_kl is a union of classes C_j = g^j H_e (H_e the e-th powers in Z_p^*, j ∈ Z_e), possibly with 0.

Equivalently, the affine group {x ↦ mx + t : m ∈ H_e} acts on all four bases. A four-context KS set with such a
graph is K4-free (Corollary 12a), so Route A applies: d = p ≥ 18 and |S_kl| ≥ 9 (Lemma 14). In particular every S_kl
contains a whole class once (p − 1)/e ≥ 9, and always for e ≤ 2.

**Lemma 23 (cyclotomic extension property).** Let a_1, …, a_m ∈ F_p be distinct, m ≤ 3, and j_1, …, j_m ∈ Z_e. The
number N of x ∉ {a_i} with x − a_i ∈ C_{j_i} for all i satisfies

    N ≥ (p − 3e² − (e − 1)²(2e + 1)√p) / e³        (m = 3),

and N > 0 already for smaller p when m ≤ 2. For e = 2 the sharper bound N ≥ (p − 15 − 2√p)/8 holds, positive for
p > 25.

*Proof.* Write 1_{C_j}(y) = (1/e) Σ_r ψ^r(y) ψ^{−r}(g^j), with ψ a character of order e, and expand the product.
The trivial term gives p − m. A term with s ≥ 1 nontrivial factors is a character sum of a polynomial with s distinct
roots, each of multiplicity < e:
* it vanishes for s = 1;
* it has modulus ≤ √p for s = 2;
* it has modulus ≤ 2√p for s = 3 (Weil);
* the excluded points x = a_i change each term by at most 1 per point.

For e = 2, Σ_x χ((x − a)(x − b)) = −1 exactly, which gives the sharper bound. ∎

**Theorem 24 (no quadratic-residue symmetry).** For no prime p does a four-context KS set have its
non-orthogonality graph in F(p, 2).

*Proof.* d ≥ 18 excludes p ≤ 17, so let p ≥ 19. Each S_kl contains a class, QR or QNR. A transversal K4 through
(3, 0) is a triple
* a ∈ S_30, b ∈ S_31, c ∈ S_32, with
* b − a ∈ S_01, c − a ∈ S_02 and c − b ∈ S_12.

Choose a in its class. Choose b with b ∈ S_31 and b − a ∈ S_01, which is Lemma 23 with m = 2. Choose c, which is
m = 3 and needs p > 25. So for p ≥ 29 every member has a transversal K4 and is not K4-free. For p = 19 and 23 the
exhaustive scan (`ks/multiplier_classify.py 2`) finds a K4 in every single-class pattern, hence in every member. ∎

**Proposition 25 (multiplier families are finite).** For every e, a member of F(p, e) satisfying (DEG) can be K4-free
only if p − (e − 1)²(2e + 1)√p − 3e² ≤ 0, that is, p < p₀(e):

| e | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| p₀(e) | 46 | 838 | 6657 | 31126 | 105841 |

Same proof as Theorem 24. Union patterns contain single-class patterns, so a single-class scan up to p₀(e) is
complete.

**Complete scans** (`ks/multiplier_classify.py`). These list the single-class patterns without a transversal K4,
up to multiplying by a generator:
* **e = 2:** only p ≤ 17 (excluded by d ≥ 18). This gives Theorem 24 a second time, by exhaustion.
* **e = 3:** p = 19, 31, 37, 43, and at p = 61, 67, 79, 127 only the pattern S_kl = H₃ for all pairs. There are none
  for 128 ≤ p ≤ 838.
  * Every instance (p, 3) with 19 ≤ p ≤ 151 is UNSAT in the corrected multiplier run (§3.12).
  * **Corollary 26: no four-context KS set has its graph in F(p, 3), for any p.** This is computer-assisted: it
    relies on the uncertified SAT verdicts for p ≤ 151.
* **e = 4:** many patterns for p ≤ 193; at p = 233, 281, 313, 337, 353, 457 exactly one, S_kl = H₄ for all pairs;
  none for 458 ≤ p ≤ 6657.
  * Of these, the patterns that also pass (C2), (C4), (TRI), (DEG) and ϑ = 3 at all four bases are the six of §3.12
    and the all-H₄ pattern at p = 281, 337, 353 and 457.
  * The all-H₄ pattern fails at p = 233 (ϑ = 2.653) and at p = 313 (ϑ = 2.918).
  * ϑ comes from the h 3×3 blocks of `search/ftheta.py`: a neighbourhood is a tri-Cayley graph over H₄ ≅ Z_h.
* **e = 4, unions.** A union pattern is K4-free only if every single-class choice inside it is, so its class sets
  form a box inside the K4-free set. Enumerating all boxes, with or without 0 in each S_kl
  (`multiplier_classify.py unions`):
  * at p = 29, 137 and 233–457, nothing passes beyond the single-class patterns above;
  * at the other primes up to 173, the boxes leave 1275 multi-class patterns (some with 0) that pass the exact
    checks: 36 at p = 37, 399 at 41, 372 at 53, 186 at 61, 189 at 73, 45 at 89, 24 at 97, 18 at 113 and 6 at 173.
  * **Every one of them fails ϑ = 3** at some basis. The test is the orbit-reduced `theta_orbits`: free H₄-orbits
    are circulant blocks, 0 is a fixed point, and it was checked against the dense SDP. At p = 37 and 41 this agrees
    with the SAT verdicts.

**Theorem 28 (computer-assisted). F(p, 4) has exactly ten skeletons, up to isomorphism, that pass every necessary
condition we have:** p = 97, 137, 157, 173, 193 (two), 281, 337, 353 and 457. The exact conditions are decided
exactly; ϑ = 3 is a numerical SDP value. The all-H₄ members are the most symmetric possible: all four bases are
equivalent under graph automorphisms.

**Beyond e = 4: the families grow** (night of 30 Oct; logs in `ks/multiplier_runs/`). Isomorphism is basis
permutations together with multiplying by a generator.
* **e = 5.** The single-class scan is complete up to p₀ = 31126; K4-free patterns end at p = 941.
  * Single-class survivors of every condition: **ten classes**, at p = 191, 211 (two), 241, 251 and 281, and the
    all-H₅ pattern at p = 631, 701, 761 and 941.
  * Union patterns add more. Each is a single-class pattern with one or two connection sets enlarged to two classes:
    p = 131 (two classes), 151, 191, 211 (two), 241, 251 and 271 (one each). At p = 131, 151 and 271 no single-class
    pattern survives.
  * No union pattern passes at p = 61, 71, 101 or 181. At p = 31 and 41 only unions can meet (DEG), and both are
    UNSAT in the corrected SAT run. **So the e = 5 list is complete: 10 single-class and 9 union-type skeletons,
    19 in all.**
* **e = 6.** The single-class scan is complete up to p₀ = 105841; K4-free patterns end at p = 3457.
  * Single-class survivors: **76 classes** over 42 primes from 223 to 3457. For p ≥ 769 there is one per prime.
  * Union patterns at p ≥ 283 add 35 classes over 13 primes from 307 to 661.
  * Union patterns for p < 283 were still being enumerated in the morning; the boxes at p = 277 alone run to
    hundreds of thousands. Where only unions can meet (DEG): p = 19 is UNSAT in the SAT run; p = 31, 37 and 43
    are undecided.
* So the multiplier families are empty for e ≤ 3, give ten skeletons for e = 4, and grow with e. All of these pass
  every necessary condition of §2–3, which therefore are far from sufficient.

**Realizability: calibration first** (`ks/realize_search.py`).
* The control is four bases of (ℂ²)⁶, basis k a Hadamard on a set Q_k of qubits with |Q_k Δ Q_l| = 4. It has
  d = 64 and density 1/4, the same density as the candidates, and is flat with h = 16. Results on it:
  * Douglas–Rachford and RRR (β = 0.5) stall at residual 2.35 and 2.85. **So the DR failures in §3.12 carry no
    information.**
  * A flat-phase search solves it in one start. That search minimises ‖S² − 2√h S − 3h I‖² over the edge phases of
    S with L-BFGS; G = I + S/√h is the Gram matrix of four bases iff S² = 2√h S + 3h I.
  * A general search solves it in 17 of 40 starts. Here bases 0–2 are free matrices, basis 3 is the standard basis,
    and L-BFGS minimises the orthonormality error plus the forbidden overlaps. The failed starts end at objective
    1–3.4.
* On the candidates, every start of every search stalls far from zero. The spread across starts is tiny, which
  points to a genuine obstruction rather than bad luck. Logs are in `ks/realize_runs/`.

| candidate | flat search: starts, ‖G² − 4G‖ best–worst | general search: starts, objective best (median) |
|---|---|---|
| e = 4, p = 97 | 30, 53.39–53.57 | 735, 58.98 (60.37) |
| e = 4, p = 137 | 30, 63.95–64.26 | 150, 81.87 (84.15) |
| e = 4, p = 157 | 10, 69.09–69.17 | 40, 96.92 (97.81) |
| e = 4, p = 193 (matching shape) | 8, 77.06–77.20 | 30, 119.0 (120.1) |
| e = 4, p = 193 (other) | 8, 76.48–76.64 | — |
| e = 5, p = 191 | (no flat realization) | 9, 117.9 (119.7) |
| e = 5, p = 211 | (no flat realization) | 10, 136.1 (137.5) |
| control (d = 64) | 1, 0 (solved) | 40, 5.6·10⁻¹⁶ (17 of 40 solved; failures 1–3.4) |

**Proposition 27 (no flat realization at p = 97 or 173).** Neither member has four bases in which every nonzero
overlap between different bases has modulus 1/√h.

*Proof.* For p = 173, some non-orthogonal pairs have only t = 6 common neighbours in a basis. Flatness needs
t/h ≥ 1/√h, i.e. t ≥ √43 ≈ 6.56. Contradiction.

For p = 97 (h = 24, and t = 5 common neighbours per basis for every non-orthogonal pair), fix a ray u and rephase
every ray in N(u) so that ⟨u|z⟩ = 1/√h > 0.
1. **Edges inside N(u) are almost real.** Take adjacent y, c ∈ N(u), and let B be the basis of c, which is neither
   u's nor y's. Completeness of B gives ⟨u|y⟩ = Σ_{c′ ∈ B} ⟨u|c′⟩⟨c′|y⟩. The five nonzero terms come from the
   common neighbours c′ of u and y in B, one of which is c. They make √h⟨c′|y⟩ five unit complex numbers that sum to
   √24 ≈ 4.899, so each has real part ≥ 0.899, i.e. phase within ±26° of 0. In particular ⟨c|y⟩ does.
2. **The obstruction.** Now take u = (0, 0), x = (3, 2) and x′ = (3, 95). Both are in N(u) and they are orthogonal
   (same basis). In basis 2 their only common neighbours are (2, 45) and (2, 52), both in N(u). Completeness of
   basis 2 gives
   0 = ⟨x|x′⟩ = Σ_{y ∈ {(2,45), (2,52)}} ⟨x|y⟩⟨y|x′⟩.
   By step 1 each term has phase within ±52° and modulus 1/h, so its real part is ≥ cos 52°/h > 0. Two such terms
   cannot sum to zero. ∎

**Flat realizations across the families.** The bound t ≥ √h alone (t = common neighbours of a non-orthogonal
pair in one basis, minimised over all pairs) rules out flat realizations for:
* p = 173 at e = 4;
* 4 of the 10 single-class classes at e = 5 (p = 191, 211 twice, 241);
* 45 of the 76 at e = 6.

Flat realizations remain possible only for:
* e = 4: p = 137, 157, 193 (both), 281, 337, 353, 457;
* e = 5: p = 251, 281, 631, 701, 761, 941;
* e = 6: 31 classes, at p = 277, 397, 421, 457, 487, 601, 661, 673 and the single-class primes ≥ 769.

(`ks/multiplier_runs/tmin_e.py`.)

`ks/flat_obstruction.py` checks this count for every u, every orthogonal pair x, x′ ∈ N(u) and every third basis; at
u = (0, 0) alone there are 12 contradictions. The same count gives no contradiction at p = 137, where the phase bound is ±34° and cos 68° is too small; for p ≥ 157 the bound
gives no information. Flatness was the natural guess from the spectrum (Proposition 18). **A realization of the
p = 97 member would therefore have to be non-flat, and the general search finds none.**

**The rank-constrained hierarchy** (references/, `ks/rankhier.py`). Level 2 with PPT (Eq. 25 of the paper) is
implemented in the complex case, as index permutations of one variable solved with SCS.
* Validation on C₅ (ϑ = 2.236, true dimension 3): k = 2 and 2.2 are infeasible; k = 2.5, 3 and 4 are feasible.
* On the Shrikhande pattern of §3.10 (16 rays, which are not four bases of ℂ⁴ because their local graphs have
  ϑ = 2), level 2 is **feasible at k = 4**. It stays feasible when the degree-2 identities of complete bases
  (Σ_{m∈B} Γ_im Γ_mj = Γ_ij) are added. As a sanity check, k = 3 and 3.9 (below ϑ = 4) are being tested for infeasibility.
* So level 2 is weaker than the local ϑ test, which our candidates pass. Level 2 sees Γ ⊗ Γ, while the local test
  uses degree-3 products ⟨v|b⟩⟨b|b′⟩⟨b′|v⟩.
* Level 3 lives on (ℂⁿ)^{⊗3} with n = 4p ≥ 388, out of reach even with the Z_p ⋊ H₄ symmetry. The hierarchy
  is the right framework for certificates, but not a practical one for these candidates.

**Next steps.**
1. Decide the p = 97 member in general, or find a non-flat realization. The flat case is closed.
2. Extend Proposition 27 to p = 137, perhaps by combining several rays u, and to the all-H₄ members.
3. e = 5, 6: run the scan up to p₀(e) = 31126 and 105841 (vectorised K4 test, feasible on a cluster).

## Files

* `search/verify_all.py`: exact verification of the 39- and 37-vertex examples (runs in under a second).
* `search/h33.py`, `search/h33.log`, `search/h33.g6`: H33 (AG(2,3)) with exact checks: rank-20 PSD Gram matrix, state,
  3-colourings, and the uniqueness step of Theorem 5.
* `search/vertex_identity.py`, `search/vertex_identity.log`: Lemma 3's identity and inequality on H33, H36, H37, H39.
* `search/bounds30.py`, `search/bounds30.log`: the exact spectral table of Corollary 6 (largest n for each d).
* `search/affine_ghz.py`, `search/affine_ghz.log`: §1.8–1.9 (the H(q, r) family, PSD criterion, dimensions;
  qudit phase-space dictionary for q = 3, 5, 7).
* `search/rigid_events.py`, `search/rigid_events.log`: Theorem 21 (exact LP, MILP, sharpness in H33, the identity
  of the remark on H36 and H37).
* `search/minimal37.py`: exact dual certificates showing H37 is vertex-critical.
* `search/ftheta.py`, `search/fsearch.py`, `search/tcount.py`: the Z_q tri-Cayley search.
* `search/groupsearch.py`: tri-Cayley search over non-cyclic groups.
* `search/gcert.py`, `search/dualcert.py`: rational primal and dual certificate generators.
* `search/shrink*.py`, `search/anneal.py`: vertex deletion and local search.
* `search/lowrank.py`, `search/check24.py`: low-dimensional realizations (gradient search).
* `search/face2.py`: optimal-face dimension and extreme-point ranks.
* `search/rigid.c`, `search/check.py`: exhaustive small-n search (n ≤ 13).
* `search/k4opt.py`: four-basis K4-free search.
* `ks/ghzgeom.py`, `ks/ghz_gadget.py`, `ks/ghz_gadget2.py`: GHZ-gadget completions and SAT check.
* `ks/covariant4.py`: covariant (phase-function) four-basis search.
* `ks/alphabet6.py`, `ks/ks4search.py`: d = 6 alphabet search with exact internal-basis SAT check.
* `ks/ks21.py`, `ks/ks21_cover5.py`: the 21-ray d = 6 KS set, its clique cover number 5, and the 5-context cover.
* `ks/fourbasis.py`: shared helpers (orthogonality graph, all bases, SAT colouring, K4s, the averaging LP of Lemma 12).
* `ks/reductions_check.py`: checks for Lemmas 7, 8 and 11 (deletion, monotonicity, sharpness) and criticality of the 21-ray set.
* `ks/pauli4.py`: Theorem 13 check: the Mermin star and random 4-tuples of three-qubit Lagrangian eigenbases.
* `ks/smalld_sat.py`, `ks/realize.py`: SAT relaxation of the four-context problem in small d, and numerical realizability of a pattern.
* `ks/cayley4.py`, `ks/cayley4_run.py`, `ks/groups.py`, `ks/cyclo4.py`: §3.12 (Cayley and multiplier families of
  K4-free four-colourable graphs, SAT + ϑ-CEGAR, complement colouring, realizability; resumable driver; group tables;
  brute-force cyclotomic cross-check). Overnight run: `ks/cayley4_runs/`, `ks/cayley4.log`, `ks/cayley4_results.jsonl`.
  Corrected multiplier run: `ks/cayley4_mult/` (states, log, `summary.txt`, `cyclo4_e2.log`). The p = 137 candidate
  and its checks: `ks/cayley4_candidate_p137.json`; the family members: `ks/cayley4_candidate_p{97,137,157,173,193a,193b}.json`, found by
  `ks/cayley4_mult/siblings_e4.py` (log `siblings_e4.log`).
* `ks/smalld_cegar.py`: second, independent encoding (selection variables, lazy completion/span clauses).
* `ks/smalld_soundness.py`: the geometric clauses hold for patterns of real configurations; `ks/smalld_d5.log`: run log.
* `ks/srg_check.py`, `ks/srg_check.log`: §3.10 checks (SRG parameter families, rook/Shrikhande, regular-graph spectral bounds).
* `ks/h33_attempt.py`, `ks/h33_attempt.log`: §3.11 driver (Theorem 19 check, structure, all attempts). Modules:
  `ks/h33model.py` (explicit coordinates), `ks/h33ctx.py` (NCF and contextual-state search), `ks/h33design.py`
  (block-aligned designs), `ks/h33dual.py` and `ks/h33w.py` (dual copies, the involution W), `ks/h33group.py`,
  `ks/h33orbit.py`, `ks/h33omega.py` (the group and its 10000-ray orbit), `ks/h33copies.py` (second-copy search).

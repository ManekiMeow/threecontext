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

## Summary

| question | result |
|---|---|
| (1) smaller graph than Perkel (57 vertices) | **36 vertices** (exact): an 8-regular, vertex-transitive tri-Cayley graph over Z2×Z6 with integer spectrum {8, 2^16, 0^9, −4^10}. It is certified by the **unweighted** Hoffman bound and has a real realization in dimension 26 in which every overlap is ±1/4 and every event has probability 1/12. Also found: a 37-vertex graph that has realizations of rank 32 (and approximately 24), and a 39-vertex graph over Z13. See Section 1.6 on dimensions. |
| (2) minimality | **Theorem: every three-context GHZ graph has at least 28 vertices**, and every realization needs dimension at least 17 (exact proofs, Section 2). So Perkel is not minimal, and n_min ∈ [28, 36]. The 36- and 37-vertex graphs are vertex-critical (exact). 36 is the minimum within every symmetric family searched: tri-Cayley graphs over all groups of order ≤ 12, all Cayley graphs of order 30 and 33. A full proof of minimality is still open. |
| (3) four-context KS set | Not found. Sufficient criterion: four orthonormal bases whose non-orthogonality graph is K4-free. Consequences in that setting: every ray is non-orthogonal to at least 28 rays in the other bases, and the dimension must be at least 17. Numerical searches found nothing. Suggested routes are listed in Section 3. |

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

### 1.2 The 36-vertex graph (new record)

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

So for Perkel the answer to the paper's open question ("we do not know if a lower-dimensional realization exists")
is no among real realizations. **37 is forced**, because the optimal Gram matrix is unique, as expected for an
edge-transitive graph where symmetrization loses nothing. The same holds for H36 with 26. Complex realizations were
not analysed: a complex Gram matrix could add antisymmetric imaginary parts on edges.

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

**What is still missing for a full minimality proof.** Perturbing Theorem 2 (n = 28, …) needs a quantitative
version of the integrality argument. The identity `n − 27 = 4 Σ_i ||E g_i||²` (with `E = F − 3/2 − 3/2·ψψ†`) gives
per-vertex bounds `||E g_i||² ≥ c_i (1/2 − 3/2·v)² / (v(1−v))`, where `v` is the weight of the neighbourhood of i
in a single other class. These vanish when every neighbourhood weighs exactly 1/3, so they do not by themselves
push the bound past 28. An exhaustive graph search is infeasible beyond about 14 vertices (`search/rigid.c`
handles n ≤ 13 using a rigidity condition derived below).

**Rigidity (used for pruning).** In a vertex-minimal example every vertex has positive weight. Hence in **every**
proper 3-colouring of H, every vertex has neighbours of both other colours. Otherwise moving it to another colour
class would force its weight to 0. No triangle-free graph on ≤ 11 vertices is rigid in this sense, and the rigid
ones on 12–13 vertices have theta ≤ 2.5.

---

## 3. Four-context Kochen–Specker sets

**Reduction.** Take four orthonormal bases `B_1, …, B_4` of `C^D` such that no 4-tuple (one vector from each) is pairwise
non-orthogonal. That is, their 4-partite non-orthogonality graph is K4-free. Then `B_1 ∪ … ∪ B_4` is a KS set.

*Proof:* a KS colouring must pick one vector per basis, and the picked vectors must be pairwise non-orthogonal.

The converse holds when the union contains no further bases. In general, extra bases inside the union could help,
and this route stays open (see below).

**Necessary local structure (K4-free case).** Fix `w ∈ B_4`. Its neighbours `N_k(w) ⊂ B_k` (k = 1, 2, 3) are
orthonormal, `w ∈ span N_k(w)` (since B_k is complete), and the non-orthogonality graph on `N(w)` is triangle-free.
So `(w, N(w))` is a **three-context GHZ paradox with state w**. By Section 2:

* every ray is non-orthogonal to **at least 28** rays of the other three bases;
* `D ≥ 17`.

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
   * The smallest KS sets have 18 rays in d = 4 (proven minimal) and at least 22–24 rays in d = 3 (computer-assisted lower bounds). Both exceed 4d, so **d ≥ 5**.
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

**Most promising routes now.**
1. Several GHZ states in the same three bases. One could look for a non-abelian symmetry (e.g. the multiplier
   group of Z13 acting on the H39 realization), mapping ψ to ψ' ⊥ ψ while permuting the completed contexts.
2. Cubic or higher-degree phase functions over larger or non-abelian groups (the covariant formulation above),
   searched with the K4 count as the objective.
3. (Route B) Larger alphabets ({0, ±1, ±i, ω^k}) in d = 6–8, using the exact SAT pipeline in `ks/ks4search.py`.

## Files

* `search/verify_all.py`: exact verification of the 39- and 37-vertex examples (runs in under a second).
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

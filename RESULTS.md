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
| (3) four-context KS set | Not found. Sufficient criterion: four orthonormal bases whose non-orthogonality graph is K4-free. Consequences in that setting: every ray is non-orthogonal to at least 28 rays in the other bases, and the dimension must be at least 17. Numerical searches found nothing. **Reductions (Section 3.7):** rays outside every basis are redundant, so the problem is exactly "four disjoint orthonormal bases with a KS union". Only internal bases meeting ≥ 3 of the four can block a K4. A symmetry group transitive on all four bases forces the K4-free case. Four Pauli eigenbases are never KS. In the K4-free case every cross overlap is at most 1/3. The path forward is in Section 3.8. |

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
  at least **9 rays of each other basis** (the per-basis version of the bound of 28).
* **Route B:** every cross pair with overlap above 1/3 lies in a K4.
* The same argument, with ψ in place of the fourth basis, applies to any three-context GHZ paradox: every event has
  probability ≤ 1/9, so every context has at least 9 events. The known examples comply: 1/12 (H36), 1/13 (H39), 1/19
  (Perkel), and a maximum of 0.094 for H37.

**Updated picture.**

| | Route A (K4-free) | Route B (K4s exist, all blocked) |
|---|---|---|
| internal bases | irrelevant | blockers must meet ≥ 3 covering bases (Lemma 11) |
| symmetry | allowed; covariant ⇒ KS iff one ray's neighbourhood is triangle-free | no group transitive on all four bases (Cor. 12a) |
| each ray | neighbourhood is a 3-context GHZ graph with ≥ 28 rays, ≥ 9 in each other basis; all cross overlaps ≤ 1/3 (Lemma 14) | lies in a K4 if d ≤ 16; p(w) ∉ STAB(G_w); pairs with overlap > 1/3 lie in a K4 |
| dimension | d ≥ 17 | d ≥ 5 |
| Pauli eigenbases (Lagrangians) | impossible (Thm 13) | impossible (Thm 13) |

Parity proofs are not excluded once internal bases are counted. With disjoint covering bases, each ray lies in
1 + m(r) bases, where m(r) counts internal bases. A parity proof therefore needs every m(r) odd and an odd number M of
internal bases. Double counting, Σ_r (1 + m(r)) = d(4 + M), then forces d even. Every ray lies in an internal basis, so
M ≥ 4, i.e. M ≥ 5.

### 3.8 Path forward

The two routes now live in different regimes. Route B is asymmetric and may sit in small dimension. Route A sits in
d ≥ 17 and allows symmetry, under which KS reduces to a single-ray test. Recommended order:

1. **Route B in d = 5 and d = 6: exhaustive computer-assisted classification.** After §3.7 the problem in fixed d is
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
   * **Pilot, d = 5.**
     * The instance has about 19k variables and 2.06M clauses.
     * Without the unique-completion and span clauses it is satisfiable in 0.4 s, by a model that completes the same
       4-clique twice, which is geometrically impossible. So the geometry carries the difficulty.
     * With all clauses, fix one K4 (WLOG the last ray of each basis) and an internal basis blocking it, made of the
       first rays of the bases it meets. Its piece sizes are 3+1+1, 2+2+1 or 2+1+1+1, the only ≥ 3-piece shapes in
       d = 5 (`--case 311 | 221 | 2111`, with lex-leader breaking of the remaining symmetry).
       * **Case 311 is UNSAT** (CaDiCaL, 493 s).
       * Cases 221 and 2111 were still running after 20 min (`ks/smalld_d5.log`).
       * If both are UNSAT, no four-context KS set exists in ℂ⁵, so d ≥ 6.
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
   * Outcome: either a four-context KS set with 20 or 24 rays, or a proof that every four-context KS set has d ≥ 7.
     Either is new; the current bound is d ≥ 5.
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
* `ks/fourbasis.py`: shared helpers (orthogonality graph, all bases, SAT colouring, K4s, the averaging LP of Lemma 12).
* `ks/reductions_check.py`: checks for Lemmas 7, 8 and 11 (deletion, monotonicity, sharpness) and criticality of the 21-ray set.
* `ks/pauli4.py`: Theorem 13 check: the Mermin star and random 4-tuples of three-qubit Lagrangian eigenbases.
* `ks/smalld_sat.py`, `ks/realize.py`: SAT relaxation of the four-context problem in small d, and numerical realizability of a pattern.
* `ks/smalld_cegar.py`: second, independent encoding (selection variables, lazy completion/span clauses).
* `ks/smalld_soundness.py`: the geometric clauses hold for patterns of real configurations; `ks/smalld_d5.log`: run log.

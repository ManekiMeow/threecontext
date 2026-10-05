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
| (1) smaller graph than Perkel (57 vertices) | **36 vertices** (exact): an 8-regular, vertex-transitive tri-Cayley graph over Z2×Z6 with integer spectrum {8, 2^16, 0^9, −4^10}. It is certified by the **unweighted** Hoffman bound and has a real realization in dimension 26 in which every overlap is ±1/4 and every event has probability 1/12. Also found: a 37-vertex graph that fits in dimension **24** (numerical), and a 39-vertex graph over Z13. |
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
* **Quantum realization.** Take the Gram matrix `I + A/4` (PSD of rank 36 − 10 = **26**). The state ψ is the normalized
  sum of the vectors of any one context. Then every event has probability exactly 1/12, every non-orthogonal pair has
  overlap 1/4, and each context sums to 1. Numerically, no real realization exists in dimension 24 or lower
  (`search/lowrank.py`).
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
* Numerically it realizes with real vectors in **dimension 24** (residual about 1e-9). Dimensions 22 and 23 fail.
  This is the lowest dimension we found, compared with 37 for Perkel.

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

**Most promising routes.**

1. **Clique-cover of known KS sets.** A KS set stays KS when vectors are added. So a KS set qualifies iff its
   orthogonality graph can be covered by 4 cliques, which requires |V| ≤ 4D. Among the known minimal sets only the
   21-ray, 7-basis set in D = 6 (Lisoněk–Badziąg–Portillo–Cabello 2014) passes the counting test (21 ≤ 24). If its orthogonality
   graph were exactly the line graph of K7 (each ray in exactly two bases and no other orthogonalities), the best
   4-clique cover would reach only 18 rays. So the question is whether its explicit vectors have **extra
   orthogonalities**. I could not download the paper from this environment (arXiv is blocked here).
2. **Symmetric "lifting" of the new tri-Cayley paradoxes.** Find a fourth basis of `C^d` all of whose vectors play
   the role of ψ for completions of `I_1, I_2, I_3`. This needs a group of unitaries preserving the three
   (completed) contexts that moves ψ to an orthonormal basis. In the Z13 examples the translations fix ψ, so a
   larger group (e.g. extending by the multiplier group of Z13) is needed.

## Files

* `search/verify_all.py`: exact verification of the 39- and 37-vertex examples (runs in under a second).
* `search/minimal37.py`: exact dual certificates showing H37 is vertex-critical.
* `search/ftheta.py`, `search/fsearch.py`, `search/tcount.py`: the Z_q tri-Cayley search.
* `search/groupsearch.py`: tri-Cayley search over non-cyclic groups.
* `search/gcert.py`, `search/dualcert.py`: rational primal and dual certificate generators.
* `search/shrink*.py`, `search/anneal.py`: vertex deletion and local search.
* `search/lowrank.py`: low-dimensional realizations.
* `search/rigid.c`, `search/check.py`: exhaustive small-n search (n ≤ 13).
* `search/k4opt.py`: four-basis K4-free search.

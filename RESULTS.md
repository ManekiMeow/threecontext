# Smaller three-context GHZ graphs, lower bounds, and the four-context Kochen–Specker problem

Follow-up to *Exploring the boundary of quantum correlations with a time-domain optical processor*
(arXiv:2208.07794; Sci. Adv. 2025). Everything numerical is reproducible from `search/`; every
claim marked **(exact)** is checked in rational arithmetic by `search/verify_all.py` or
`search/minimal37.py`.

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
| (1) smaller graph than Perkel (57 vertices) | **37 vertices** (exact). Also a highly symmetric 39-vertex graph with an exact certificate that uses only the weights 1/4 and 1/3. The 37-vertex graph has a real realization in dimension **25** (numerical), compared with 37 for Perkel. |
| (2) minimality | **Theorem: every three-context GHZ graph has at least 28 vertices**, and every realization needs Hilbert-space dimension at least 17 (both exact proofs below). Perkel is not minimal. The 37-vertex graph is **vertex-critical** (exact): deleting any vertex drops theta below 2.986. It is also the smallest example in all the symmetric families searched. The gap 28 ≤ n_min ≤ 37 is still open. |
| (3) four-context KS set | Not found. Sufficient criterion: four orthonormal bases whose non-orthogonality graph is K4-free. Consequences in that setting: every ray is non-orthogonal to at least 28 rays in the other bases, and the dimension must be at least 17. Numerical searches (real, D ≤ 24) found nothing. Suggested routes are listed below. |

---

## 1. Smaller examples

### 1.1 Tri-Cayley graphs

Write `T(q; S0, S1, S2)` for the graph with vertex set `Z3 x Zq`, where `(a, x) ~ (a+1, x+s)` for every `s` in `S_a`.

* It is automatically 3-partite, and the classes `a = 0, 1, 2` are the three contexts.
* It is triangle-free iff `0 ∉ S0 + S1 + S2 (mod q)`.
* The Perkel graph is `T(19; S, 7S, 49S)` with `|S| = 3`.

The translations `x -> x+1` act as automorphisms. So the Lovász SDP can be symmetrized and splits into
q Hermitian 3x3 blocks (`search/ftheta.py`), which makes the whole family cheap to scan.

**Uniform weights (|S_a| = 3).** If each `|S_a| = 3`, the optimal weights must be uniform. Then `theta(Hbar) = 3`
iff every Fourier block `M_k` (`k != 0`) has smallest eigenvalue `>= -1`. A full scan
(`search/tricayley333.py`) shows that **q = 19 is the smallest q that works**, so the Perkel graph is the
smallest *6-regular uniform* member of this family.

**General weights.** With non-uniform weights, smaller q works.

**The 39-vertex graph** `T(13; {0,1,3,9}, {0,1,10}, {1,6,8})`

* 130 edges, degrees 6 and 7, triangle-free, girth 4.
* Certificate **(exact)**: let `W = 1/4 A_{01} + 1/3 A_{12} + 1/3 A_{20}`, where `A_{a,a+1}` is the bipartite adjacency
  matrix between consecutive classes. Then:
  * every vertex has total W-weight exactly 1 into each other class, so `W 1 = 2·1`;
  * `W + I` is PSD (spectrum of W: `2, 0.7138^(12), 0.2438^(12), -0.9576^(12), -1^(2)`).

  By Lovász's theorem, `theta(Hbar) >= 1 - lambda_max/lambda_min = 1 + 2/1 = 3`.
* Quantum realization in dimension 37 (rank of `W+I`). Each of the 39 events has probability exactly 1/13,
  so each context sums to 1 and noncontextual models predict 0 in at least one context.

**Smaller tri-Cayley graphs.** Within the Z_q family:

* q = 13 also gives a second, non-isomorphic graph, `T(13; {0,1,3,9}, {0,1,4}, {1,2,5,7})` (143 edges).
* q = 14 gives several more (42 vertices).
* q ≤ 12 gives nothing for any connection-set sizes (see `search/fsall.log` once the sweep finishes).
* For the non-cyclic groups of order 9–12 (Z3², D5, Z2×Z6, D6, A4, Dic3), see `search/groups.log`.

### 1.2 The 37-vertex graph

Let `H37 = T(13; {0,1,3,9}, {0,1,4}, {1,2,5,7})` with the two vertices `(0,5)` and `(2,11)` removed.

* 37 vertices with context sizes 12, 13, 12. 129 edges, degrees 6/7/8, triangle-free.
* **(exact)** `theta(Hbar37) = 3`. Rational certificate: `search/shr_q13_1_0139_014_1257_cert.json`, with weights
  `c_i > 0` and `Y_ij` on edges. Here `B = Y + diag(c)` is PSD, `B_ij = 0` on non-edges, and `sum(B)/tr(B) = 3`.
  The event probabilities are `c_i`.
* **(exact)** It is vertex-critical: for every vertex v, a rational dual certificate shows
  `theta(complement(H37 - v)) <= 2.9857 < 3` (`search/minimal37.py`).
* For the 39-vertex graph, the edges are needed too: deleting any one of its 10 edge orbits drops theta to between 2.77 and 2.91.
* Numerically, H37 is realizable with real vectors in **dimension 25** (`search/lowrank.py`; residual about 1e-9).
  Dimension 22 fails, and 23/24 are under test. This is a large saving over the 37 dimensions used in the
  experiment.

Randomized greedy vertex deletion from the 39-vertex graph and from the q = 14 graphs never went below 37.
Annealing over edge-maximal tripartite triangle-free graphs on 36 vertices has not reached theta = 3 so far
(best 2.989).

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

* `search/k4opt.py` minimizes `Σ_{4-tuples} Π_{pairs} |<v_a|v_b>|²` over four real orthonormal bases.
  The minimum stays at about 40–50% of the Haar-random value for D = 8, 12, which is consistent with D ≥ 17.
  Larger D is still running.
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

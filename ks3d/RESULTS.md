# The smallest Kochen–Specker set in dimension 3

Question: what is the minimum number of rays in a Kochen–Specker (KS) set in R^3?

Definition used throughout (the standard one in the SAT literature): a finite set S of rays is a
KS set iff it has no 010-colouring, i.e. no map c: S -> {0,1} with no two orthogonal rays both 1
and every orthogonal basis (orthogonal triple) inside S containing a ray coloured 1. Size = |S|.

## 1. State of the art (checked October 2026)

| bound | value | source |
|---|---|---|
| smallest known KS set | **31 rays** | Conway–Kochen (unpublished, reported in Peres, *Quantum Theory: Concepts and Methods*, 1993) |
| proven lower bound | **24 rays** | Li, Bright, Ganesh, *A SAT Solver + Computer Algebra Attack on the Minimum Kochen–Specker Problem*, IJCAI-24 (arXiv:2306.13319), and independently Kirchweger, Peitl, Szeider, *Co-Certificate Learning with SAT Modulo Symmetries*, IJCAI-23 |
| earlier lower bounds | 22 (Uijlen–Westerbaan 2016), 18 (Arends–Ouaknine–Wampler 2011) | |

So the screenshot's "22" is out of date: the bound has been **24** since 2023, with a
40 TiB proof certificate for order 23 (Li–Bright–Ganesh). Li–Bright–Ganesh also show 24 holds for
complex rays. Szeider has estimated that pushing the bound to 25 with SAT modulo symmetries needs
roughly 100–200 CPU-years. I found no claim of 25 or more as of this date (searches of arXiv
listings via search engines; arxiv.org itself is not reachable from this environment). The
true minimum is open, in **[24, 31]**.

### 1.1 Li, Bright, Trandafir, Cabello, Ganesh, arXiv:2604.19947 (April 2026)

*SAT + nauty: orderly generation of small KS sets containing the smallest SI-C set.*
Does not move the bounds (still [24, 31]). What it adds:

* Exhaustive result: up to 33 rays, Schütte's 33-ray set is the **only** KS set containing the
  complete 25-ray SI-C set (the cross-product closure of the 13-ray Yu–Oh set); 1,641 CPU-hours,
  13 TiB DRAT certificate. Corollary: **no KS set with ≤ 30 rays contains that 25-ray core**.
* Recursive Canonical Labelling (RCL) on top of nauty: hereditary canonical form for orderly
  generation at ~27 ms per check independent of order, versus minutes for the lexicographic
  checks used by SAT+CAS and SMS at orders 23–26. This removes the isomorphism bottleneck of the
  graph-enumeration route.
* A geometric propagator: coordinates of new rays are derived by cross products during the search,
  and partial graphs that cannot be embedded are cut immediately with small clauses.
* Code: github.com/BrianLi009/SAT-nauty (reachable from this environment).

## 2. What is realistic here

This container has 4 cores. The exhaustive graph-enumeration route to 25 costs ~100 CPU-years, so
it is out of reach. Instead this work attacks the problem from the **upper-bound side** on
well-defined finite sub-problems: fix a natural universe V of rays (all rays with small
coordinates over Z or Z[√d]) and determine **exactly** the smallest KS set contained in V. A
result "min over V < 31" would be a new record; "min over V = 31" rules out a natural place where a
smaller set could hide.

## 3. Method (`ks3d/`)

* `rays.py`, `qrays.py`: exact ray universes (integers / Q(√d) with Fractions), orthogonal pairs,
  bases.
* Implicit hitting sets / CEGAR (`minks.py` with CP-SAT master, `cegar.py` with incremental
  CaDiCaL master). Key fact: for **any** 0/1 assignment c of the whole universe, every KS set
  S ⊆ V must contain a witness against c (an orthogonal pair coloured 1,1 or a basis coloured
  0,0,0). Each such clause is a valid cut. The master picks S with |S| ≤ T hitting all cuts; a SAT
  oracle either colours S (new cuts, from that colouring greedily extended to V) or proves S is a
  KS set. If the master becomes infeasible, V contains no KS set of size ≤ T (**exact**, modulo
  solver correctness).
* Extra master constraints, valid for every inclusion-minimal KS set (so no solution is lost):
  every ray of S lies in a basis inside S, and has ≥ 3 orthogonal rays inside S (a ray of degree
  ≤ 2 can always be coloured last).
* Symmetry: all 48 signed coordinate permutations act on V; every cut is added with all its images
  and the master carries lex-leader constraints. Sound because the master's feasible set is then
  invariant under the group.
* `verify.py`: independent re-check of found sets (exact orthogonality, Glucose instead of CaDiCaL,
  vertex-criticality).

Cross-checks: the CP-SAT and incremental-SAT masters agree on heights 2 and 3; height 2 reproduces
a 31-ray set and height 1 / Q(√2) reproduces Peres's 33-ray set.

## 4. Results

| universe V (rays spanned by vectors with coordinates ...) | rays in a basis | bases | result | time |
|---|---|---|---|---|
| integers in [-1, 1] | 9 | 4 | no KS set | <1 s |
| integers in [-2, 2] | 49 | 26 | **min = 31** (exact) | <1 s |
| integers in [-3, 3] | 97 | 50 | **min = 31** (exact) | 2 s |
| integers in [-4, 4] | 205 | 110 | **no KS set with ≤ 30 rays** (exact) | 9 min |
| integers in [-5, 5] | 421 | 266 | ≤ 30: **undecided**, master did not finish 20 iterations in ~5 h | |
| a + b√2, a, b ∈ {-1,0,1} | 201 | 88 | **min = 33** (exact; Peres's set is optimal here) | 15 min |
| a + b√2, a ∈ [-2,2], b ∈ {-1,0,1} | 481 | 258 | ≤ 30: **undecided**, same | |
| a + b√3, a, b ∈ {-1,0,1} | 177 | 84 | KS sets exist, **none with ≤ 33 rays** (exact) | 2 min |

The 31-ray sets found (`found_sets.json`) are uncolourable and vertex-critical by `verify.py`, with
71 orthogonal pairs and 17 bases. They are presumably the Conway–Kochen set (not yet checked for
isomorphism).

**Conclusion so far (exact, computer-assisted):** no KS set with fewer than 31 rays can be built
from rays with integer coordinates of absolute value ≤ 4, from {0,±1,±√2}-type coordinates
a+b√2 with |a|,|b| ≤ 1, or from a+b√3 with |a|,|b| ≤ 1.

**Conjecture (not proved):** 31 is the true minimum. The evidence above is consistent with it but
does not touch the general problem: a smaller KS set could use rays with larger or other algebraic
coordinates.

## 5. Next steps

* Height 5 and Q(√2) with |a| ≤ 2 are beyond the current master: the hitting-set SAT problem
  itself becomes hard. Needs case splitting on a fixed basis orbit plus parallel runs, or a
  stronger master (e.g. cuts from many diverse near-colourings up front).
* Run each universe at T = 23..30 split by a fixed basis orbit to parallelise.
* A different route that does touch the general problem: replicate the Li–Bright–Ganesh / SMS
  orderly generation to 24 on this machine as a calibration, and estimate the cost of 25 with the
  cube-and-conquer splitting; that number decides whether a cloud run is worth proposing.

## Reproduce

```
pip install python-sat ortools
cd ks3d
python3 cegar.py 4 30            # integer height 4, target <= 30
python3 fieldks.py 2 1 1 32      # Q(sqrt2), |a|,|b| <= 1, target <= 32
python3 verify.py found_sets.json
```

## 6. SAT + nauty timing test (October 9, 2026)

The paper's solver (github.com/BrianLi009/SAT-nauty, commit a9577f2) was built here (4 cores).
Each run is one exact order, single-threaded. Runner for arbitrary seeds: `ks3d/satnauty_seed.py`
(the seed must be in its own RCL canonical order, or the solver blocks it immediately and
reports a vacuous UNSAT. The runner relabels it with `verifiers/rcl_canon`.)

| seed | order | result | canonical partial graphs | wall time |
|---|---|---|---|---|
| 25-ray SI-C core | 26–29 | no KS graph | 0 | < 1 s |
| 25-ray SI-C core | 30 | no KS graph | 4 | 1 s |
| 25-ray SI-C core | 31 | no KS graph | 72 | 8 s |
| Yu–Oh 13 rays | 20 | no KS graph (3,993 colourable candidates rejected) | 31,842 | 56 s |
| Yu–Oh 13 rays | 21, 22, 23 | **unfinished after 2 h 20 min** | | > 8,000 s |

These reproduce the paper for orders 26–31. With the Yu–Oh seed the cost grows by more than
100× from order 20 to 21, so orders 24–30 with this seed are far beyond this machine. Even a cluster
would need a better split, e.g. Trandafir–Cabello's suggestion to seed Yu–Oh and forbid the full
25-ray completion.

Lovász theta screen (`theta.py`). ϑ(G) < ν(G) (ν = max number of disjoint bases) is a valid
sufficient test for KS. But on the 31-ray set α = 11, ϑ = 11.71, ν = 7, and the colourable
30-ray subsets have almost the same values, so the screen cannot separate KS from non-KS graphs.
Cost: ϑ SDP (cvxopt) 79 ms vs SAT colouring 0.3 ms on 31 vertices.
In fact the screen can never fire on a graph that has rays: ϑ(G) ≥ n/3 ≥ ν(G) (Theorem 7.1 (b)
with w = 1). Section 7 shows that no spectral test of this kind can work.

### 6.1 Cube-split timing for the Yu–Oh seed (`estimate.py`)

Cubes fix the seed neighbourhoods of the first two new vertices and the edge between them
(141 feasible neighbourhoods, 39,762 cubes). Times are CPU-seconds on this machine; each cube
costs ~0.015 s solver start-up, i.e. ~600 CPU-s of overhead per order.

| order | unsplit run | all cubes (CPU-s) | hardest cubes | result |
|---|---|---|---|---|
| 19 | 0.9 s | 590 | < 0.1 s | no KS set |
| 20 | 56 s | 684 | 3.3 s | no KS set |
| 21 | > 8,700 s (killed) | **2,168** | 140 s, 225 s, 267 s | no KS set (exact, all cubes UNSAT) |
| 22 | | **≥ 39,700** (4,563 for 39,744 finished cubes; 18 cubes unfinished: 15 at > 900 s, 3 at > 7,200 s) | > 2 h | 1 candidate graph, **not embeddable**; 18 cubes undecided (stopped) |

The hard cubes are the ones where both new vertices have **no** Yu–Oh neighbour: the cost lives
in structures far from the seed, so seeding helps less and less as the order grows. For the same
reason, forbidding the full 25-ray completion (`yuoh13x`) cannot change orders < 25 at all (25
vertices are needed to contain it). From 25 on it removes only the branch the paper showed is cheap
(25-core extensions to order 31 take 8 s).

Order-22 candidate (cube 286): 22 vertices, 49 edges, 12 bases, min degree 3, not 010-colourable.
Forcing coordinates by cross products from the fixed Yu–Oh vectors determines all 22 rays, and 7
required orthogonalities then fail (e.g. vertex 13 = (1,-3,-2) is adjacent to vertex 15 = (1,-2,-3)
but not orthogonal to it). So the graph has no orthogonal representation extending Yu–Oh. The
solver runs with `--unembeddable-check 0` and reports such graphs for a post-check, as in the paper
(43 of its 44 order-33 candidates failed the same check). Checker: `check_solution.py`.

### 6.2 Cost projection for "no KS set with ≤ 30 rays contains Yu–Oh"

Measured growth per added ray: 20 → 21 about 30× (56 s → 2,168 s), 21 → 22 at least 18×.
Orders ≤ 23 are already implied by the global bound 24, so the new content is orders 24–30, the
most expensive ones. Projecting order 30 alone from order 21 (2.2·10³ s) with a constant
per-order factor r:

| r | order 30 (CPU-hours) |
|---|---|
| 18 (current) | ~10^11 |
| 10 | ~6·10^8 |
| 5 | ~10^6 |

Even the optimistic r = 5 means about a million CPU-hours for order 30 alone, plus orders 24–29. A
~20,000-core cluster would need weeks for that, and r = 10 would need decades. **Conclusion: a
brute-force Yu–Oh-seeded proof up to 30 rays is not a realistic cluster project with this
solver.** The cost sits in cubes whose new rays have no Yu–Oh neighbour, i.e. the effectively
unseeded part of the problem. Forbidding the 25-ray completion does not touch that part.
Caveat: r was measured at only two steps and may fall at higher orders. Finishing orders 22–24
on a cluster (hours to days) would measure that before committing a large allocation.

## 7. Graph spectra and the minimum (October 27, 2026)

Question: does the spectrum of the orthogonality (exclusivity) graph G, or of its complement Ḡ
(the non-orthogonality graph), say anything about the smallest KS set?

**Answer: not about the KS property itself, provably; a little about embeddability.** Every
spectral bound, and every degree-2 semidefinite bound, is satisfied by the rays themselves, for
every set of rays in ℂ³, KS or not (Theorem 7.1, Proposition 7.2). The KS contradiction first becomes visible to algebraic certificates of
degree 4 (Theorem 7.3), and spectral bounds have degree 2. Spectral tests can only help on the
other side of the problem, ruling out graphs that have no rays (§7.4). Even there they are weak:
* unweighted spectra catch none of the 17 minimal obstructions;
* a degree-2 SDP with propagation catches 4;
* a degree-3 SDP catches 8 of the 13 that are obstructions over ℂ.

Notation: R is a set of pairwise non-collinear rays in ℂ³ with orthogonality graph G, T the set
of its orthogonal bases (= triangles of G), m = |T|, Π_v the projector onto ray v, and
w_T(v) = number of bases through v. R is KS iff α(G, w_T) < m. (An independent set meets each
basis at most once, and w_T(S) counts the bases it meets, so a 010-colouring is exactly an
independent set of w_T-weight m.)

### 7.1 The barrier

**Theorem 7.1.** Let m ≥ 1.

(a) *Exclusivity graph.* ϑ(Ḡ) = ξ_f(G) = 3, where ξ_f is the projective rank.
* Consequently, for every Hermitian matrix B supported on the edges of G and every p ≥ 0,
  `1 + max(E⁺_p(B)/E⁻_p(B), E⁻_p(B)/E⁺_p(B))^{1/|p−1|} ≤ 3`, where E^±_p are the positive and
  negative p-energies.
* This family covers p = 0 (the inertia bound), p = 2 (Ando–Lin) and p → ∞ (Hoffman). The
  inequality is Elphick–Tang–Zhang (arXiv:2504.01295, Theorem 1.8 and Remark 1.9).
* By Lovász, the best Hoffman bound over all edge weightings B of G is exactly ϑ(Ḡ) = 3.
* A KS set, on the other hand, has χ(G) ≥ 4 and χ_f(G) ≥ 3m/(m−1) > 3.

(b) *Complement / independence side.* For every weight w ≥ 0, Schrijver's ϑ′(G, w) ≥ w(V)/3. For
the basis weights, ϑ(G, w_T) = ϑ′(G, w_T) = m, whereas α(G, w_T) = m if R is 010-colourable and
α(G, w_T) ≤ m − 1 if R is KS.

(c) Hence no bound of the following kinds can certify that R is KS, or tell a KS set from a
non-KS set:
* Hoffman, inertia or p-energy lower bounds on χ(G), with any weights.
* Upper bounds on α(G, w) that dominate ϑ′ (Lovász, Schrijver, Hoffman's ratio bound, Delsarte's LP
  bound).
* Cvetković inertia upper bounds on α(G, w), with any weights.

  Here w is any basis-generated weight w = Σ_T λ_T 1_T, λ ≥ 0. Every 010-colouring S has
  w(S) = w(V)/3 exactly, and all these bounds are ≥ w(V)/3.

  On the complement this says, by Lovász's formula
  `ϑ(G, w) = min { λ_max(A) : A_ii = w_i, A_ij = √(w_i w_j) on the edges of Ḡ }`, that every
  matrix carrying the weights √w_T √w_Tᵀ on Ḡ and the diagonal has λ_max ≥ m. The classical value
  of a KS set is at most m − 1.

*Proof.*
* (a) The projectors give a 3/1-representation, so ξ_f ≤ 3. The vectors vec(Π_v) with handle
  vec(I)/√3 form an orthonormal representation of Ḡ with value 3, so ϑ(Ḡ) ≤ 3. Both quantities are
  at least ω(G) = 3.
* The χ_f bound: χ_f(G) ≥ w_T(V)/α(G, w_T) ≥ 3m/(m−1).
* (b) Let Y_ij = |⟨v_i|v_j⟩|² and D = diag(√w). The matrix X = DYD / w(V) is psd (Schur product of
  psd Gram matrices) and entrywise ≥ 0, has trace 1, and vanishes on edges. Its value is
  ‖Σ w_i Π_i‖_F² / w(V) ≥ (tr Σ w_i Π_i)² / (3 w(V)) = w(V)/3.
* For w_T, Σ w_T(v) Π_v = Σ_T Σ_{t∈T} Π_t = m I gives equality. The m bases form a fractional
  clique cover of weight m, so ϑ(G, w_T) ≤ m.
* (c) Ratio-type bounds dominate ϑ′. For the inertia bound, apply Cvetković to the blow-up G_w,
  where each v becomes w_v non-adjacent copies; this is again a multiset of rays, so (a) applies.
  The p = 0 case of (a) gives n⁺ ≤ 2n⁻ and n⁻ ≤ 2n⁺, hence α-bound = |G_w| − max(n⁺, n⁻) ≥ |G_w|/3.
  ∎

**Proposition 7.2 (the unifying reason: degree 2).** For any state ψ, define the functional L by
L(1) = 1, L(x_i) = ⟨ψ|Π_i|ψ⟩ and L(x_i x_j) = Re⟨ψ|Π_iΠ_j|ψ⟩. Then L satisfies every constraint of
degree ≤ 2 of the 010-colouring system x_i² = x_i, x_i x_j = 0 (ij ∈ E), Σ_{t∈T} x_t = 1, because
Σ_{t∈T} Π_t = I. Its moment matrix [[1, p], [p, X]] is psd, since X − ppᵀ is the real part of the
Gram matrix of the vectors (Π_i − p_i)ψ. So no degree-2 Nullstellensatz or sum-of-squares
certificate of non-colourability exists. Ratio-type eigenvalue bounds are certificates of exactly
this kind; the inertia bounds are not, and are handled by (a).

The numbers bear this out (`spectral.py`, `spectral.log`):

| set | KS | χ_f | Hoffman | inertia | p = 2 | ϑ(Ḡ) | ϑ(G, w_T) | m | α(G, w_T) |
|---|---|---|---|---|---|---|---|---|---|
| 31 rays (§4, height 2) | yes | 16/5 | 2.869 | 2.067 | 2.242 | 3 | 17 | 17 | **16** |
| Peres 33 | yes | 13/4 | 2.801 | 2.200 | 2.185 | 3 | 16 | 16 | **15** |
| Yu–Oh 13 | no | 35/11 | 2.702 | 2.167 | 2.101 | 3 | 4 | 4 | 4 |
| Yu–Oh closure, 25 rays | no | 35/11 | **2.879** | 2.200 | 2.412 | 3 | 16 | 16 | 16 |
| 31 rays minus ray 0 (3-colourable) | no | 3 | 2.746 | 2.143 | 2.240 | 3 | 14 | 14 | 14 |
| 31 rays minus ray 10 (not 3-colourable) | no | 16/5 | 2.842 | 2.143 | 2.242 | 3 | 16 | 16 | 16 |

* Every spectral bound is ≤ 3.
* The colourable 25-ray set has a larger Hoffman bound than the 31-ray KS set.
* ϑ(G, w_T) = m in every row, to within 3e-7. The KS property shows only in the exact integer
  α(G, w_T).
* Random Hermitian edge weightings (2000 per set) never exceed 2.66.
* The complement's unweighted spectrum is no better: its Hoffman bound for the clique cover number
  is 7.99 (31 rays) and 7.70 (Peres), against true values 13 and 12.

### 7.2 Where the KS contradiction becomes visible: degree 4

Work in R = ℚ[x]/(x_i² − x_i, x_i x_j : ij ∈ E). Its monomials are x_S for independent sets S. A
Nullstellensatz refutation of degree d is an identity 1 = Σ c_{S,T} x_S g_T with |S| ≤ d − 1 and
g_T = Σ_{t∈T} x_t − 1.

**Theorem 7.3.** For every finite set of rays in ℂ³, the 010-colouring system has no
Nullstellensatz refutation of degree ≤ 3 over any field of characteristic 0.

*Proof.* Take the maximally mixed state: L(x_S) = Re tr(Π_{s₁}⋯Π_{s_k})/3 for |S| = k ≤ 3.
* For k ≤ 3 the real part of tr(Π_aΠ_bΠ_c) = ⟨a|b⟩⟨b|c⟩⟨c|a⟩ is symmetric in a, b, c. Cyclic
  shifts leave it unchanged, and a transposition conjugates it. So L is well defined on monomials.
* L vanishes on dependent sets, because the cyclic product contains a zero inner product.
* For |S| ≤ 2 and a basis T disjoint from S, Σ_{t∈T} L(x_{S∪t}) = Re tr(Π_S Σ_t Π_t)/3 = L(x_S),
  and the dropped terms (t orthogonal to S) are 0 on both sides.
* Rows with S ∩ T ≠ ∅ are 0 in R.
* So L kills every row of degree ≤ 3, while L(1) = 1. ∎

For four factors the trace is no longer symmetric, and the bound is attained (`nsdegree.py`,
`nsdegree.log`):

| set | over ℚ (exact rational elimination) | GF(2) | GF(3) | GF(p), two 31-bit primes |
|---|---|---|---|---|
| 31 rays | **degree 4** (none of degree 3) | 4 | 3 | 4 |
| Peres 33 | **degree 4** (none of degree 3) | 3 | 4 | 4 |
| Yu–Oh closure 25 (colourable) | none (any degree) | none | none | none |

Consequences:
* A degree-4 Nullstellensatz refutation implies a degree-4 sum-of-squares refutation (Lasserre
  level 2), and level 1 never refutes.
* The KS contradiction in dimension 3 is therefore exactly a degree-4 phenomenon for both known
  small sets, and at least degree 4 for every KS set, including any unknown one with 24–30 rays.
* The characteristic-0 assumption matters. The proof uses 1/3 and real parts, and over GF(2) or
  GF(3) degree 3 sometimes suffices.

### 7.3 What the complement does see: the dimension

The Gram matrix of R is I + W, where W is a weighting of the edges of Ḡ (non-zero exactly there)
with eigenvalue −1 of multiplicity n − 3. In other words, R exists iff Ḡ has minimum semidefinite
rank ≤ 3, realised by pairwise non-collinear vectors. This is the only place where "3" enters, and
it is a rank condition. Convex (spectral) relaxations see it only through cliques: ω(G) ≤ 3 and
ϑ(Ḡ) ≤ 3. A 4-cycle, which no set of rays contains (two rays have exactly one common orthogonal
ray), passes every spectral test.

### 7.4 Embeddability tests on the 17 minimal unembeddable graphs

Li–Bright–Ganesh block every KS candidate that contains one of 17 minimal unembeddable graphs,
with 2, 5 and 10 graphs on 10, 11 and 12 vertices. `mugs.json` is the published list from
github.com/BrianLi009/SAT-nauty, `log_archive/min_nonembed_graph_10-12.txt`.
`unembeddable.py` (nauty geng, numerical embedding, then Z3) independently reproduces the two
order-10 graphs and finds none below order 10. Z3 refutes `I?qa``eeO`; for `ICOedPKL?` it times
out and the hand proof below applies.

The degree-2 "Born" test (`born.py`) applies to rays in ℂ³ as well as ℝ³. It uses
X = (3Y − J)/2 with Y_ij = |⟨v_i|v_j⟩|². X is the Gram matrix of the normalised traceless parts of
the projectors, so it satisfies:
* X ⪰ 0 and X_ii = 1;
* X_ij = −1/2 on edges and X_ij ≥ −1/2 everywhere;
* X_ia + X_ib ≤ 1/2 for every edge ab, because two orthonormal rays carry total probability ≤ 1.

If every feasible X has X_ij = −1/2 on a non-edge, the rays i, j are orthogonal in every
embedding, so the edge is added. If a 4-cycle appears, the graph is refuted.

| test | refutes |
|---|---|
| Hoffman, inertia, p = 2 energy (unweighted) | 0 of 17 (values 2.00–2.63) |
| ϑ(Ḡ) ≤ 3 | 0 of 17 |
| Born SDP, distinct rays only | 0 of 17 |
| Born SDP + forced-orthogonality propagation | **4 of 17**: `I?qa``eeO`, `` JsOH?`JH_i_ ``, `` Kko_GC`C?b`q ``, `` K`?KA@se_YCX `` |
| degree-3 Born SDP + propagation | **8 of 17**: the four above and `` J@`ICC[Gv_? ``, `JRPCSGoAgJ_`, `KqCaC?[I_J?Z`, `` K`?LAQOP@KbK `` (three of the eight by forced equal rays) |
| anything valid over ℂ | at most 13 of 17: the other 4 (`` J?GX?dH`f_? ``, `K_GTAQOP@KbK`, `` K_GP`GWBACoK ``, `Kt?G?DKOoqCo`) are embeddable in ℂ³ (numerical complex embeddings) |

The degree-3 test (`born.py --degree3`, `born3.log`) adds Z_ijk = Re⟨i|j⟩⟨j|k⟩⟨k|i⟩, fully symmetric.
* For each ray i, the Gram matrix of |i⟩ and the vectors Π_j|i⟩ is
  [[1, Y_i·], [Y_i·, W_i]] with Re W_i[j, k] = Z_ijk and W_i[j, j] = Y_ij. It must be psd, vanish on
  edges jk, and satisfy Σ_{t∈T} Z_ijt = Y_ij.
* This is the degree-2 moment matrix of the state v_i. The symmetry of Z ties the states together;
  it is the same cyclic-trace symmetry as in Theorem 7.3.

Full tables: `born.log`, `born3.log`. Every refuted graph is among the 13 listed in the solver
source (`cadical-rcl/src/unembeddable_graphs.h`, marked * in the logs). The other 4 are embeddable
in ℂ³, so no test valid over ℂ can refute them. Of the 13, degree 2 refutes 4 and degree 3
refutes 8. Neither refutes the cubic graph below.

Example: the cubic order-10 graph `ICOedPKL?` passes the degree-2 and degree-3 tests. Its
obstruction is a holonomy (proof over ℝ).
* Vertex 7 = v has neighbours a = 0, b = 4, c = 2 in the plane v^⊥, one from each of the bases
  {0, 3, 6}, {4, 1, 8} and {2, 9, 5}.
* Put v = e₃ and a′ = e₃ × a. Then rays 3 and 6 are cos θ₁ e₃ + sin θ₁ a′ and
  −sin θ₁ e₃ + cos θ₁ a′; likewise θ₂ for rays 1, 8 (with b) and θ₃ for rays 5, 9 (with c).
* Let c_ab = ⟨a, b⟩ = ⟨a′, b′⟩. The cross edges 6–1, 8–5, 3–9 give tan θ₁ = c_ab tan θ₂,
  tan θ₂ = c_bc tan θ₃ and tan θ₃ = c_ac tan θ₁.
* If some tan θ is 0 or ∞, a ray coincides with v. Otherwise multiplying gives
  c_ab c_bc c_ac = 1, so |c| = 1 and a, b, c are the same ray.
* The other order-10 graph `I?qa``eeO` is the same mechanism with two planes 8^⊥ and 9^⊥ and a
  cycle of four cross edges between them (c⁴ = 1, c the cosine between the planes). Degree 2
  happens to see that case.

### 7.5 Why the GHZ3 method does not transfer

In the three-context GHZ problem (main `RESULTS.md` §2) spectral arguments were decisive:
* The non-orthogonality graph H is triangle-free, so tr(Π₁Π₂Π₃) = 0.
* That gives the cubic identity Σ f(f−1)(f−2) = 0 for F = Π₁ + Π₂ + Π₃.
* The state ψ is a common eigenvector (Fψ = 3ψ).

There the classical contradiction is free (triangle-freeness) and the hard part is quantum
realizability, a Gram-matrix question, which is exactly what spectra measure. In 3D KS the roles
are swapped:
* Realizability is free (any rays).
* The hard part is the classical side, which is invisible below degree 4.
* Every context is a complete basis, so the analogue of F is Σ_T Π_T = m I, which carries no
  information.
* No state is distinguished.
* Ḡ is dense and full of triangles, so there is no trace identity.

### 7.6 Consequences for the bounds [24, 31]

* Spectral or SDP relaxations of degree ≤ 3 can never replace the colourability check of the SAT
  searches.
* At degree 2 and 3 they can only pre-filter embeddability, and only partly: 4 and 8 of the 13
  complex minimal obstructions. Holonomy-type obstructions such as `ICOedPKL?` pass both.
* A non-computational lower bound would need certificates of degree ≥ 4 (for example Lasserre
  level 2 on candidate graphs, or a structural degree-4 argument). Alternatively it needs
  algebraic or geometric embeddability arguments of the holonomy type above.

Files: `spectral.py` / `spectral.log` (Theorem 7.1 checks, tables), `nsdegree.py` /
`nsdegree.log` (Theorem 7.3, degree computations), `born.py` / `born.log` / `born3.log` /
`mugs.json` (§7.4),
`unembeddable.py` / `unembeddable.log` / `mugs_recomputed.json` (order ≤ 10 recomputation, 533
graphs at order 10; needs nauty's `geng` via `GENG=...`).

## Related: the three-context GHZ problem

The minimum here, between 24 and 31 rays, resembles the three-context GHZ minimum of between 30 and 33 events in the
main `RESULTS.md`. Section 1.9 there compares the two problems.
* There is no reduction either way: GHZ paradoxes need dimension ≥ 18.
* Peres's 33 rays split 9 + 12 + 12 like H33, but with a different orthogonality pattern.
* The real link of H33 is to *complex* dimension 3. H33 is the operator lift of the qutrit Hesse configuration, the
  same 21-ray SIC-plus-MUB structure that gives state-independent contextuality in Bengtsson–Blanchfield–Cabello
  (arXiv:1109.6514).

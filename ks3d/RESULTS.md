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

### 6.1 Cube-split timing for the Yu–Oh seed (`estimate.py`)

Cubes fix the seed neighbourhoods of the first two new vertices and the edge between them
(141 feasible neighbourhoods, 39,762 cubes). Times are CPU-seconds on this machine; each cube
costs ~0.015 s solver start-up, i.e. ~600 CPU-s of overhead per order.

| order | unsplit run | all cubes (CPU-s) | hardest cubes | result |
|---|---|---|---|---|
| 19 | 0.9 s | 590 | < 0.1 s | no KS set |
| 20 | 56 s | 684 | 3.3 s | no KS set |
| 21 | > 8,700 s (killed) | **2,168** | 140 s, 225 s, 267 s | no KS set (exact, all cubes UNSAT) |
| 22 | | running | | |

The hard cubes are the ones where both new vertices have **no** Yu–Oh neighbour: the cost lives
in structures far from the seed, so seeding helps less and less as the order grows. For the same
reason, forbidding the full 25-ray completion (`yuoh13x`) cannot change orders < 25 at all (25
vertices are needed to contain it). From 25 on it removes only the branch the paper showed is cheap
(25-core extensions to order 31 take 8 s).

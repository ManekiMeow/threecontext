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
| integers in [-5, 5] | 421 | 266 | running | |
| a + b√2, a, b ∈ {-1,0,1} | 201 | 88 | **min = 33** (exact; Peres's set is optimal here) | 15 min |
| a + b√2, a ∈ [-2,2], b ∈ {-1,0,1} | 481 | 258 | ≤ 30: running | |
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

* Finish height 5 and Q(√2) with |a| ≤ 2; then mixed fields Q(√2, √3) and Q(√5).
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

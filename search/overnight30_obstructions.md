# Exact screening of 30-event candidates

`overnight30_screen.py` stores self-contained certificates and verifies them
using standard Python. Discovery may use SAT or SDP; numerical failure is
always inconclusive. These conditions apply to a triangle-free support graph
with a prescribed proper partition into three contexts of ten vertices.

## Dimension-independent exclusions

* **Disconnected support.** A theta-three Gram matrix on a disjoint union
  restricts to theta three on at least one component: the matrix is block
  diagonal, so its objective is a trace-weighted average of component
  objectives. Every component has fewer than 30 vertices, contradicting the
  proven event bound.
* **False twins in one context.** Merge their orthogonal Gram vectors by
  addition. Their squared norms add, the total vector sum is unchanged, and
  every required zero remains a zero because their neighborhoods coincide.
  The resulting graph is triangle-free and three-colorable on 29 vertices.
* **A same-context pair without a length-three path.** This strictly
  strengthens the false-twin test. Identifying two vertices in one context
  can create a triangle only if their neighborhoods contain adjacent
  vertices, equivalently if the original vertices are joined by a
  length-three path through the other two contexts. If no such path exists,
  the same Gram-vector merge produces a triangle-free 29-event graph.
* **Any unbalanced proper three-coloring.** Equality in the three-color theta
  bound makes each color class a complete context for the GHZ state. No event
  can have zero probability, since deleting it would contradict the 30-event
  lower bound. Each class therefore has at least nine vertices. A class of
  exactly nine forces at least twelve vertices in each other class by
  Theorem 5, exceeding 30. Consequently every proper three-coloring must
  have class sizes exactly `(10,10,10)`.
* **Unequal sides in a pair-support component.** Swapping the two colors on
  that component supplies an explicit unbalanced proper coloring. This is
  the cheap special case of the preceding SAT-discovered obstruction.
* **Too few third-context bridges.** A balanced pair component with `r`
  vertices on each side needs at least `r` third-context vertices adjacent
  to both sides. A complete bipartite component has no such bridges by
  triangle-freeness and is therefore excluded, including a single isolated
  edge in a pair support.
* **An event bridging no pair component.** Every event in the third context
  must touch both sides of at least one component of the other two contexts.
  Otherwise its state probability is zero, contradicting event minimality.
* **Insufficient probability mass from cross-degrees.** The exact strict
  probability bounds `U(q,r)` described below must sum to more than one in
  each context. A sum of one or less is an exact global exclusion, with no
  symmetry or uniform-probability assumption.
* **Rational theta dual below three.** A saved rational symmetric matrix
  `D` has diagonal entries one and entries one on support edges. The verifier
  checks `t I - D` positive semidefinite over the rationals with `t < 3`.
  Thus `theta(complement(H)) <= t < 3` independently of all discovery solvers.

These obstructions remain valid when an allowed edge becomes zero in the
realization. The degree-mass argument uses the orthogonality of the two
entire allowed neighborhoods, without requiring every allowed overlap to
be nonzero.

## Strict probability bounds from the two cross-degrees

For cross-degrees `q,r`, let `a,b` be `(q-1)/9,(r-1)/9`, ordered so that
`a >= b`. Positive probabilities, context sums one and `c_i < 1/9` imply
that the two neighborhood masses strictly exceed their respective lower
bounds. The relevant degree pairs have `1 <= q,r <= 9`; a missing neighbor
context contradicts completeness, and a complete neighbor context leaves
no positive third sector.

Lemma 14 gives `c_i <= 1/(1/alpha + 1/beta + 1/rho)`, with
`rho = 1-alpha-beta > 0`. Minimizing the strictly convex reciprocal sum
over the relaxed closed bounds yields the exact upper bound

* `U = 0` if `a+b >= 1`, indicating an impossible event;
* `U = 1/9` if `a <= 1/3`;
* `U = a(1-a)/(1+3a)` if `a > 1/3` and `b <= (1-a)/2`;
* `U = 1/(1/a + 1/b + 1/(1-a-b))` otherwise.

In the last two cases the unique relaxed optimizer has the larger
neighborhood mass exactly `a`, while the actual mass is strictly larger.
Thus `c_i < U`; in the `1/9` case Corollary 32 supplies the same strictness.
Consequently each context necessarily has `sum_i U(q_i,r_i) > 1`.

The `degree_mass` certificate stores the context index and exact rational
sum. Its verifier recomputes every degree and bound with the pure-standard-
Python `event30_orbit_degree.relaxed_event_upper` helper. It rejects the
previously screened asymmetric seed because one context has upper sum
`313/315 < 1`; the C5-double leader has upper sums `(10/9,10/9,10/9)` and
passes. The built-in controls also check exclusion at exact equality:
six bounds `5/54` plus four bounds `1/9` sum to one, which is insufficient.

## Pair-component bridge inequality

For a component with sides `X,Y`, let `u = sum_X f_i` and `v = sum_Y f_j`,
where the orthogonal weighted Gram vectors in each context sum to the unit
state `psi` and have squared norms `c_i`. There are no cross-support edges
leaving this component. Hence `<u,v> = sum_X c_i = sum_Y c_j`, and
`||u-v||^2 = 0`. Write their common vector as `v` and its squared norm as
`m > 0`. In particular `<psi,v> = m`.

Let `W` be the third-context vertices adjacent to both `X` and `Y`. A vertex
outside `W` is orthogonal to at least one of the two equal expansions of
`v`. Thus the projection `Pi_C v` is supported on the orthonormal rays in
`W`. Since `psi` belongs to the third context's span, Cauchy-Schwarz gives

`m^2 = |<psi,v>|^2 <= (sum_W c_i) ||Pi_C v||^2 <= (sum_W c_i) m`.

Therefore `m <= sum_W c_i`. This argument does **not** assume that `v`
belongs to the third context's span.

The earlier Kempe condition makes `|X| = |Y| = r`. With ten events in each
context and all `0 < c_i < 1/9`, one has `m >= (r-1)/9`; the inequality is
strict when `r < 10`, while `m=1` when `r=10`. If `0 < |W| < r`, then
`sum_W c_i < |W|/9 <= (r-1)/9 <= m`, a contradiction. The case `W` empty
contradicts `m>0` directly. Consequently `|W| >= r`.

Finally the common vectors of all pair components sum to `psi`. A third-
context vertex outside every bridge set is orthogonal to every common
vector, hence to `psi`, which would give zero probability. Thus the bridge
sets must cover the entire third context.

The exact certificate methods are `bridge_shortage`, storing the pair,
component and recomputed bridge set, and `unbridged_vertex`, storing the
pair and uncovered vertex. Their verifier uses standard Python only. The
`pair_bridge_data` helper computes components and bridge sets without
applying the thirty-event cardinality restriction. The isolated B-C edge
`[10,21]` in the C2 prototype supplies a `bridge_shortage` control with empty
bridge set; H33 supplies positive bridge-coverage controls.

## Additional conditions only in dimension 19

Theorem 35 excludes the exceptional shape in Lemma 33. Every remaining pair
of reduced context spaces has nine principal angles strictly between zero
and pi/2, and the full context spaces share exactly the state direction.

* **Connected pair support.** Every bipartite support component supplies a
  nonzero common vector obtained by completing the state inside that
  component. Distinct components supply orthogonal common vectors, so there
  can be only one component.
* **Perfect matching in each pair support.** The full 10-by-10 cross-context
  Gram matrix has singular values one and the nine positive principal-angle
  cosines. It is invertible. Its determinant expansion therefore contains
  a nonzero term, which is a perfect matching in its support. A Hall-deficient
  vertex subset is a compact exact obstruction.

The program stores these separately under `dimension19_certificate` and
`dimension19_eligible`. They do **not** change the global status to `excluded`.
Passing them does not prove any realization or rank bound.

## Contraction lemma and a single-transposition action

Let `f_i` be weighted Gram vectors for a theta-three realization, normalized
so that the sum in each context is the same unit state `psi`. Within each
context the vectors are orthogonal. Partition each context into arbitrary
nonempty groups, and replace a group `O` by `F_O = sum_{i in O} f_i`.
Orthogonality gives `||F_O||^2 = sum_{i in O} ||f_i||^2`, so total trace is
unchanged. The three context sums and their total squared norm are also
unchanged. The new Gram matrix is positive semidefinite and supported on
the quotient graph, where two groups are adjacent when some original pair
between them is adjacent. If this quotient is triangle-free, it is another
three-context GHZ graph. Thus a quotient with fewer than thirty vertices is
impossible. No averaging or projection onto a symmetry eigenspace is needed.

In particular, suppose an involution swaps exactly one pair in each context
and fixes the other twenty-four vertices (`c2_onepair`). Contract its three
moving pairs. The quotient has twenty-seven vertices, so it must contain a
triangle. A quotient triangle involving any fixed vertex lifts to a triangle
in the original graph: adjacency to a moving pair is all-or-none, and any
edge joining the other two quotient vertices then completes the lift.
Consequently all three moving-pair quotient vertices must be mutually
adjacent.

Each invariant bipartite block between two moving pairs is empty, one of
the two perfect matchings, or the full `K_{2,2}`. None is empty by the
preceding paragraph. A full block, together with the other two nonempty
invariant blocks, gives an original triangle: choose either vertex in the
third moving pair and one of its neighbors in each endpoint pair. Therefore
all three blocks are perfect matchings. Their union is either two triangles
or one six-cycle; triangle-freeness forces the six-cycle.

Hence **every possible `c2_onepair` realization must induce exactly `C6` on
its six moved vertices**. This is a necessary support condition, not an
exclusion of the whole action. The three moving pairs are the opposite
pairs of that cycle. The broader contraction test applies without symmetry
and requires a length-three path for every same-context vertex pair.

An immediate simpler corollary is that **each pair-context adjacency block
has distinct rows and distinct columns**. If two vertices `u,v` in one
context have identical neighborhoods in another context, any proposed
length-three path between them uses one vertex of that other context.
Its adjacency to both `u` and `v` closes an original triangle, contradicting
triangle-freeness. Thus the required path cannot exist.

Consequently the color-preserving automorphism group acts faithfully on
each context. An automorphism fixing one context pointwise leaves every
vertex's neighborhood in it unchanged. Distinct rows force every vertex
in both remaining contexts to be fixed as well. In particular, the order
of any color-preserving automorphism equals the order of its restriction
to each of the three contexts.

The cheap contraction test excludes seven of the eight initial weak-action
family bests; only the `c5_double` best passes it. Those candidates already
have independent rational theta exclusions, but the new witnesses avoid SDP.
Six of those seven failures already have identical neighborhoods in just
one other context. The `c3_fixed` best fails the full length-three-path test
without any repeated row or column in any pair block, so the full test is
strictly stronger on these samples.

## Initial checks and examples

The built-in controls check false twins, a SAT recoloring witness, a
dimension-19 Hall witness, compact graph6 decoding, a rational dual witness,
and the exact H33 boundary primal. H33 is deliberately outside the 30-event
structural screen; its exact theta-three certificate verifies normally.

A smoke scan of 33 previous saved candidates produced two Kempe-recoloring
exclusions and 23 rational theta exclusions. Eight candidates were left
inconclusive by the selected numerical-score threshold; this was not an
exhaustive run. The asymmetric seed711 candidate with numerical theta
2.9122784989 has a proper coloring of sizes `(9,10,11)`, so its exclusion does
not depend on a numerical SDP. Its B-C support also has four components.

The existing rational-primal discovery helper requires a strictly positive
margin and a positive probability floor. It can miss an exact feasible
boundary face, as the H33 control illustrates. Its failure is never treated
as nonexistence.

## API and command line

From another script in `search/`, import `fast_screen` for inexpensive
structural checks or `screen_graph` for the optional SAT and rational
certificate discovery stages. Both accept `edges` and `parts`, or compact
`graph6` / `g6` records. Compact 30-vertex records default to contiguous
contexts `0..9`, `10..19`, and `20..29`.

```
python search/overnight30_screen.py self-test
python search/overnight30_screen.py scan POOL --output-dir CERTS --seconds 60
python search/overnight30_screen.py scan TOP --output-dir CERTS --exact-dual --min-theta 2.85 --recheck
python search/overnight30_screen.py verify CERT.json
```

The scan resumes by labeled graph hash. `--recheck` revisits prior records,
for example to upgrade a fast scan with exact dual certificates. Only saved
decisive witnesses support exclusion; all other statuses are inconclusive.
Archived `.jsonl.gz` pools are discovered recursively and streamed directly;
compressed JSON and graph6 files are supported as well. Incomplete trailing
JSON records or archives being replaced by the supervisor are skipped.

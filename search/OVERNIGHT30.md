# Overnight construction search: 4–5 November 2026

User-authorized objective: construct as many promising three-context GHZ support
graphs with context sizes (10,10,10) and weak or absent symmetry as possible,
and turn recurring failures into exact obstructions. Continue until **5 November
2026, 06:00 Europe/Copenhagen = 05:00 UTC**, then stop new searches and report.
The user has explicitly requested that **every commit be pushed automatically**.

## Repositories and persistence

- Code/results branch: `research/ghz-gaps-four-context` in
  `https://github.com/ManekiMeow/threecontext.git`.
- Dedicated experiment branch: `experiments/ghz30-20261104` in the same repository.
  This is a separate worktree so the background publisher cannot accidentally
  commit an agent's in-progress source edits.
- Original code workspace: `/workspace/threecontext`.
- Experiment worktree: `/workspace/ghz30-overnight`.
- Run directory within the experiment worktree:
  `search/overnight30_runs/2026-11-04/`.
- Python environment: `/workspace/threecontext-venv/bin/python`. If a scheduled
  run has a fresh environment, recreate a virtual environment and install
  `numpy scipy sympy networkx cvxpy python-sat`; use the repository's existing
  cloud network instructions and proxy configuration.
- Use `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, and `MKL_NUM_THREADS=1`.
  The current machine has a four-core CPU quota and 16 GiB RAM; the supervisor
  uses three worker processes, leaving capacity for verification and analysis.

The supervisor commits only experiment outputs in its dedicated worktree and
immediately pushes each commit. Never force-push. If a push fails, retain the
commit and resolve/retry it before making further automatic data commits.
At completion, merge the data branch into the code branch, resolve ordinary
documentation conflicts, verify and push. Do not discard unrelated user edits.

## Scientific baseline and priorities

The current global bounds remain `30 <= n_min <= 33`, `19 <= d_min <= 20`.
Keep Claude's theorem numbering in `RESULTS.md`; Theorem 38 was the last result
at campaign setup. The latest numbered result is now Theorem 57 and Corollaries 57a–b. Results
39–57 cover contraction, orbit counting, moving cycles, large-orbit exclusions,
pair-component bridges, degree probability budgets, regular pair supports and
the excluded six-orbit quotient of the current C5 numerical leader, followed
by affine-plane mixed-action exclusions, iterated probability budgets, and
the exclusion of a nine-vertex symmetry orbit in any single context, global
pair-support matchings, disjoint-neighborhood probability budgets, and complete
phase exclusions for nine specified C5 quotients, and the exclusion of even
one isoclinic reduced pair in a uniform-probability dimension-19 realization,
followed by strict weighted-degree bounds, one-sided regularity restrictions,
and zero rectangles forced by two uniform contexts in every dimension, followed
by a closed-small-neighborhood obstruction with arbitrary probabilities.
Read Sections 1.10 and 2.1–2.3 before adding claims.

Already excluded and not worth repeating as standalone searches:

- The rotated-tail dimension-19 case (Theorem 35).
- Three pairwise-isoclinic reduced context spaces (Theorem 36).
- Even one isoclinic reduced pair in dimension 19 when all thirty event
  probabilities are 1/10 (Theorem 55). The nonuniform case is not covered.
- Regular unweighted Hoffman constructions on thirty events (Proposition 37).
- Identical copies of any linear 10_3 configuration (Theorem 38, exhaustive).
- Regular tri-Cayley families over groups of order at most twelve (earlier work).
- `c5_mixed`, with orbit patterns `(5,5)`, `(5,1,1,1,1,1)`, `(5,1,1,1,1,1)`
  (Corollary 40). Five fixed vertices cannot have distinct neighborhoods that
  are unions of only two orbits. SAT independently found the new constraints
  unsatisfiable, but the exclusion rests on this counting proof.
- A common orbit pattern `(m,1,...,1)` in all three contexts for m=7,8,9,
  for any context-preserving symmetry group (Theorem 42). No unitary symmetry
  or averaging assumption is used. Order-7 and order-9 color-preserving
  automorphisms are excluded. Mixed order-eight actions involving a two-cycle
  are **not** excluded by that theorem.
- The six-orbit quotient in Proposition 46, including all cyclic shifts of
  the current 95-edge C5-double leader. The averaging argument there may
  increase rank; it is an event-number exclusion in every dimension.
- Two `(9,1)` context orbits with any fixed vertex in the third context
  (Theorem 47). This includes all C3×C3 RRT actions.
- C3×C3 RTT actions with three distinct kernel lines in both T contexts
  (Theorem 48). An exact structural reduction leaves three canonical graphs,
  all excluded by rational theta duals.
- Both two-kernel RTT cases (Theorem 50). Opposed multiplicities reduce to
  at most 27 canonical supports with exact duals. Aligned multiplicities have
  a complete finite proof: 299 structural certificates, 65 assignments covered
  by 14 rational duals, and an independently checked 222-step LRAT proof.
- Any context-preserving group with a nine-vertex orbit in even one context
  (Theorem 51). A transitive degree-nine group contains a regular C9 or C3²
  subgroup; every resulting action is now excluded. This does not imply that
  the full automorphism group has order indivisible by nine.
- The 150 C5 quotient classes selected in the 23:00 review, now all exactly
  excluded. This does not exhaust the C5-double family. Theorem 53 adds ten
  probability certificates; Theorem 54 covers every phase of the final nine
  quotients with 1,328 exact graph certificates.

Search directions, allowing nonuniform event probabilities and unequal overlap
weights unless a run explicitly records a narrower ansatz:

1. Weak C2, C3 and C5 actions with multiple orbits and fixed points, including
   different orbit decompositions in the three contexts. These are not the
   previously exhausted regular tri-Cayley actions. The initial eleven families
   include one transposition per context (only six moved vertices), two
   transpositions per context, and one three-cycle per context.
2. Asymmetric tripartite triangle-free graphs: population search, neighbourhood
   replacements, multi-edge moves, crossovers and independent restarts. Start
   from the strongest existing graphs, not only the incidence constructions.
3. Apply exact combinatorial screens and mine surviving graph patterns for a
   new obstruction. In particular, a thirty-event GHZ graph must be connected,
   have no nonadjacent false twins, and every proper three-colouring must have
   classes (10,10,10). Every same-context pair must also have a length-three
   path through the other contexts (Proposition 39). In the one-transposition
   per context family, the six moved vertices must induce C6. These are exact
   exclusions, not numerical score heuristics.
   Theorem 43 additionally requires enough third-context bridges for every
   pair component, and every event must bridge some component. Theorem 44
   requires the strict rational degree-probability budget to exceed one in
   every context. Proposition 49 propagates strict probability intervals for
   three rounds. Theorem 52 makes Hall-deficient pair supports globally
   impossible. Theorem 53 adds exact budgets for partitions into subsets whose
   cross-neighborhoods are pairwise disjoint. These are in the workers' exact
   screening/cut pipeline.
4. Promote promising supports to full Lovasz SDP checks. Preserve candidates
   near theta=3 and their numerical matrices. Seek exact rational primal
   certificates; if a valid example is found, investigate its realization rank
   and dimension 19 separately. A rational-certificate failure is not an
   exclusion, especially on a rank-deficient face.

Pair-support connectivity is necessary in dimension 19, but is NOT a universal
event-minimality filter for all dimensions. Uniform-probability degree restrictions
must likewise not be applied to unrestricted candidates. For thirty positive
events, the general bounds c_i<1/9 and neighbourhood orthogonality imply total
degree at most ten; weak-symmetry workers may use this necessary pruning rule.
Theorem 43 is a global restriction on disconnected pair supports, but does not
require full pair connectivity. Its common component vector need not belong
to the third context span; the proof uses its projection there.

## 22:00 continuation findings and next priorities

The immutable review is `search/overnight30_reviews/2026-11-04T2200/`: 69,112
scored records, 68,486 distinct labeled supports, and 18 representatives after
exact isomorphism checks within a selected set of 89. The C5-double leader
has numerical theta 2.9134341612 and exact dual upper bound 1821/625. Proposition
46 now also excludes its whole quotient class by forced averaged probabilities
and inconsistent saturated-edge overlaps. No global bound has changed.

The workers load `overnight30_shared_seeds.json` on every asymmetric chunk;
injection history prevents repeats. Conflict-clearing switches replace the
previously stalled two-edge switches. All seeds remain construction guides,
including those with known individual exclusions. The earlier asymmetric
2.8876788 seed is now excluded by the exact degree budget 313/315, and the C2
2.8797224 seed by empty component bridges; their metadata records these facts.

C5-double scoring uses the exact Fourier reduction in `overnight30_c5.py`.
The worker learns at most 100 new rational-dual subgraph cuts per chunk, with
gzip certificates verified on discovery and resume by `overnight30_fast_dual.py`.
A 20-second smoke produced 397 scores and 90 certificates. Numerical failure
never creates a cut. The full SDP remains a fallback.

## 23:00 continuation findings and next priorities

The immutable review `search/overnight30_reviews/2026-11-04T2300/` adds 107,777
labeled supports, including 93,632 newly scored labels, against the preceding
completed pool. Selected 169 labels give 125 representatives relative to this
selection and previously selected classes. Thirty-three received exact
exclusions. The strongest new selected support is asymmetric, theta about
2.9130510892, with exact upper bound 29141/10000. The overall historical leader
is unchanged. All 36 certificates, 44 isomorphism mappings and input hashes
verify with the review's `verify.py --inputs`. Two substantive numerical
failures were independently excluded; deadline endings remain ordinary bounded
search outcomes. The global bounds are unchanged.

Asymmetric chunks now use bounded, checkpointed hot excursions with clustered
partial rewires; the resume control accepted four downhill transitions and
preserved best/pool. The shared inbox has eight diverse mutation seeds.

The exact probability propagation helper is `event30_probability_propagation.py`.
Three rounds add 16 exclusions among the selected 125 representatives, including
theta 2.90776942 via budget 993319/1000000. Proofs and timing controls are in
`overnight30_probability_reviews/2026-11-04T2300/`. Rounding is upward and all
certificate replay uses rational arithmetic; no numerical bound is trusted.

The C5 quotient review in `overnight30_quotient_reviews/2026-11-04T2300/`
examines 150 of 5,273 published quotient classes and excludes 126 with exact
orbit-probability Farkas certificates. `event30_probability_lp.py` separates
numerical discovery from exact verification. `overnight30_quotient.py` verifies
this catalog and Proposition 46, and also checks regular pair supports from
Theorem 45. A separate saturation catalog in
`overnight30_saturation_reviews/2026-11-04T2300/` supplies five exact equality
transcripts, one reproducing Proposition 46 and four excluding additional
classes. The combined 131/150 selected classes expand to 5,756 forbidden
context/orbit signatures. `event30_quotient_saturation.py` uses exact affine
elimination and saturated-edge consistency; nineteen selected classes remain
unresolved. The asymmetric worker also applies this shared verified catalog.
The 60 new base auxiliaries preserve replay by rebasing legacy path cuts. The former
C5 leader is retired as a live best (historical evidence remains), and a
2.899889 guide outside the catalog is loaded. Numerical score reduction after
retiring an excluded best is expected, not lost progress.

The affine-plane investigation now completely excludes every C3² action with
a regular nine-orbit in any context. Theorem 51 extends this to any transitive
nine-orbit of any context-preserving group. The aligned proof bundle is
`event30_c3square_aligned_certificates/`; run
`python search/overnight30_c3square_search.py --verify-coverage search/event30_c3square_aligned_certificates`
for standard-library replay, adding `--rebuild-base` in the research environment
to regenerate the reviewed CNF with python-sat 1.9.dev20. Opposed and three-kernel
proofs have their own helpers and exact duals. The global event and dimension
bounds remain unchanged.

The broader auxiliary two-kernel experiment (16,973 assignment records before
the refined finite proof) is archived as a 1.7 MB tar.gz under host chunk
`00053_asymmetric/auxiliary_c3square_review/`. This immutable archive was handed
to the existing supervisor publisher without changing its Git index or control
state and was published on the experiment branch in commit `92f8064`.
It is auxiliary research data; do not count its records as new distinct
labels without deduplication against the main pools. The host chunk's normal
search remains separate and unchanged.

Next useful research steps, preserving the running supervisor:

- The nineteen residual selected quotients are now excluded by Theorems 53–54.
  Explore further published quotients or unclassified asymmetric supports.
  Preserve the actual group-action requirement and verify every new class
  certificate. Existing derivation formats and immutable proof catalogs must
  retain reproducible verification.
- If pursuing C3² further, use only three-orbit/fixed-point actions:
  `T=(3,3,3,1)` and `D=(3,3,1,1,1,1)` with faithful kernel-line sets. All actions
  involving `R=(9,1)` are excluded. General asymmetric and weak C2/C3/C5
  constructions remain the main search priority.
- If changing runner configuration, skip `c5_mixed` without reassigning family
  seeds, and consider favoring C5-double after reviewing progress. Preserve
  diversity and the absolute deadline. Source commits capture dependency
  versions; expanding per-job hashes to the Fourier/screen/seed dependencies
  would improve provenance at the next runner restart.

## 00:00 continuation findings (5 November)

The immutable review `overnight30_reviews/2026-11-05T0000/`, pinned to data
commit `5484ce3`, adds 311,477 labeled supports and 60,905 newly scored labels.
Among 130 selected leaders, 92 are new representatives relative to this
selection and earlier reviewed classes; 38 explicit mappings identify the
others. Twenty-two selected representatives have exact exclusions. The new
asymmetric leader has theta about 2.9162207091 and exact upper bound 2.9173.
The four strongest selected supports have trivial context-preserving automorphism
groups. The review verifies 24 local certificates, 14 referenced proofs,
isomorphism mappings and immutable input hashes with `verify.py --inputs`.
One numerical timeout is exactly excluded; eight deadline endings are ordinary
bounded outcomes. Global event and dimension bounds remain unchanged.

Theorem 52 upgrades the pair Hall condition to every dimension. Its certificate
method is `pair_hall`; historical `pair_hall_d19` certificates still replay.
Pair disconnection remains dimension-specific. `event30_hall.py` checks the
new filter, a genuine 8-to-7 witness, boundary cases and dominance-cut replay.

Theorem 53 uses a 3-by-3 and 4-by-4 moment relaxation to bound total probability
on subsets with disjoint neighborhoods in another context. The certificates
are rational and their trivial dual matrices are positive definite, so the
bounds are strict. In particular, five events have total mass below 15/32.
`event30_disjoint_partition.py` searches context partitions by an integer
subset DP; equality budgets are accepted only with independently verified
strictness. The ten additional quotient Farkas certificates and independent
principal-minor audit are in `overnight30_disjoint_reviews/2026-11-05T0000/`.

Theorem 54 completely covers the other nine selected quotient classes.
Singleton-tree and general proper-subset-tree translation gauges reduce their
phase constructions to 5,300 normalized assignments and 1,328 representatives
under common units. Every representative has an exact theta upper bound below
three. Proof bundles and independent gauge/enumeration checks are in
`overnight30_phase_reviews/2026-11-05T0000/` and
`overnight30_general_phase_reviews/2026-11-05T0000/`. All 150 selected classes
are now excluded; this is not an exhaustive classification of C5-double graphs.

The asymmetric worker preserves full best records across chunk resumes and
rechecks active records against new filters, retiring newly excluded records
while retaining historical graphs and scores. A targeted finite recombination
of four asymmetric leaders scored twenty previously unscored labels; the best
was about 2.91470491, below the current leader. The auxiliary archive is in
host chunk `00066_asymmetric/auxiliary_asymmetric_review/`, published in data
commit `992cbf1`. It must not be counted as twenty globally new
isomorphism classes. Numerical optimizer ranks are not realization dimensions.

The remaining C3² actions using only T=(3,3,3,1) and D=(3,3,1,1,1,1) reduce to
seven action types, checked by `event30_td_actions.py`. The saved 99-edge DDD
prism passes the earlier local screens but has exact theta bound 6383/2500.
It is an excluded individual control, not a family-wide DDD obstruction.

Theorem 55 strengthens the dimension-19 exclusion under uniform probabilities:
even one reduced isoclinic pair is impossible. Its overlap matrix would have
only three or four entries per row/column. Counting single intersections and
using exact norm identities forces four-regular support of one of two types.
One type has Petersen intersection parameters and no required 3-by-3 zero
rectangle; the other forces a cubic graph on ten vertices to decompose into
K4 components. Both are impossible. The complete analytic proof and exact
real/complex two-context positive controls are in `event30_isoclinic_uniform.py`.
Independent proof review passed. This does not assume an actual graph symmetry,
and it does not exclude nonuniform realizations with an isoclinic pair.

Theorem 56 proves strict weighted degree between two and five in every
oriented context pair. A side with constant degree must therefore have degree
three or four; the opposite side need not be regular. Every pair support has
at least twenty edges. `event30_one_sided_regular.py` verifies the Jensen,
projection and common-neighbor arguments. The global certificate method
`one_sided_regular` preserves the older `regular_pair` certificate format.

The separate bounded expansion `overnight30_quotient_expansion/2026-11-05T0000/`
uses four published C5 pools containing 5,503 quotient classes. In sixty seconds
it reviews 3,560 further classes and exactly excludes 3,378 (502 sector LP,
1,299 saturation, 1,577 disjoint-moment LP). Combined with the earlier 150,
3,528 published classes are excluded; 182 reviewed classes remain inconclusive
and 1,793 were not reached. Proofs and representative metadata are compressed
JSONL, and all exact certificates and actual C5 actions replay in standard
Python. This remains a published-pool review, not a full family classification.

The expanded catalog has 163,348 permuted signatures and 72 exact-count
auxiliaries. Startup benchmarks measured about four seconds for the new proof
catalog and under one second for CNF generation/solver bootstrap. Replaying
47,691 historical C5 models and 11,070 cuts took about 69 seconds and exhausted
a short smoke budget before new models; this is setup cost, not a solver
failure. A fresh twenty-second control generated 336 models and 81 scores.
Its only error record was a deadline-before-solve skip, with no substantive
solver failure. The raw fresh/resume controls and logs are archived under
chunk `00076_asymmetric/auxiliary_quotient_expansion_benchmarks/`, published
in experiment commit `a40db76`; the source stores only proof and benchmark
summaries.
Legacy cut rebasing from the earlier 60 auxiliaries was verified. The refreshed
2.8894779990 C5 guide has an unexcluded quotient, but its individual graph has
an exact bound below three and is explicitly marked as a mutation guide.

The 01:00 continuation below completes this published quotient census and
supersedes these residual counts. The nonuniform dimension-19 case and the
uniform case with all three pairs non-isoclinic remain open. Do not repeat the
completed 150-class phase review or weaken the actual-group-action requirement.

## 01:00 continuation findings (5 November)

The immutable review `overnight30_reviews/2026-11-05T0100/` is pinned to data
commit `5a2b454` and adds 421,977 labeled supports and 46,053 newly scored
labels. Among 128 selected leaders, 111 are new representatives relative to
that selection and previously reviewed classes, with seventeen explicit
isomorphism mappings. Twenty selected representatives have exact exclusions.
The verifier checks twenty-three certificates including supplemental checks,
five referenced proofs, input hashes and the persisted historical graph.

The historical numerical leader improved to 2.9166384145, but is now directly
excluded by a disjoint-neighborhood partition budget 15/16. Its full
automorphism group has order two; only the context-preserving group is trivial.
The active guide at this cutoff is fully asymmetric, theta about 2.9073614218,
and has an individual exact theta upper bound 2.9084. The worker correctly
retired the historical graph without losing its graph, score or certificate.
One substantive numerical timeout is exactly excluded; seven deadline skips
are ordinary budget endings. No global event or dimension bound changed.

Theorem 57 proves that two uniform ten-event contexts must have a 3-by-3 zero
rectangle in their allowed pair support, in every dimension. Otherwise the
third context becomes uniform with saturated degree types (2,4)/(4,2).
Opposite types have impossible integer intersection equations; one type alone
would require too many disjoint two-element neighborhoods. This proof is
analytic and needs no moment certificate. Related uniform restrictions allow
only one saturated type per context, at most one (4,4) event, and at most three
degree-two events (the latter uses Theorem 53's strict four-event bound).

Corollary 57a excludes an isoclinic pair with allowed closed-Petersen support
in every dimension, initially allowing nonuniform probabilities. The pair's
Gram identities force both contexts uniform, contradicting the zero-rectangle
theorem. `event30_uniform_rectangle.py` contains the reviewed proof and exact
controls. These statements do not exclude arbitrary nonuniform graphs or
arbitrary isoclinic pair supports, and no unrestricted worker filter was added.

Corollary 57b does give a global filter without uniformity. A subset S containing
an entire context is impossible if every vertex of S has a whole allowed
neighborhood of size at most two, in another context, contained in S. The maximum
probability in S obeys M≤h(2M), hence M≤1/10; the contained context becomes uniform
and weighted saturation compatibility forces impossible disjoint pairs. A
deterministic peeling procedure finds the maximal such subset. All degrees and
neighborhoods refer to the ORIGINAL support, never the graph induced after
peeling. In particular, every realization has an event with at least three
neighbors in both other contexts, giving a 3-by-3 zero rectangle in at least one
pair. This does not extend Theorem 57 to every specified nonuniform pair.
The helper `event30_degree_two_cover.py` supplies exact `degree_two_closure`
certificates in the existing global format. Both workers use the filter;
symmetry runs also save and replay compact dominance cuts, with a separate
rejection counter compatible with old checkpoints. Controls cover a full
30-vertex closure, a proper 20-vertex closure, an altered graph with empty
closure, false witnesses and exact H33/H36 positive primals. No new numerical
search was needed to establish this obstruction.

The new quotient expansion reviews all 1,879 previously unexamined classes in
the pinned published pool: 34,721 scored records give 5,589 quotient classes.
It certifies 1,843 additional exclusions (534 sector LP, 502 saturation,
807 disjoint-moment LP), leaving 218 combined residual classes. The five
highest-scoring residual classes selected from the midnight expansion are
then excluded by complete phase enumeration: 9,062 normalized assignments,
2,269 unit representatives, largest exact upper bound 2.967026795602. The
combined catalogs exclude 5,376 of 5,589 classes, leaving 213 unresolved.
Every published class was reviewed; the full C5-double family was not enumerated.
The new immutable catalogs are `overnight30_quotient_expansion/2026-11-05T0100/`
and `overnight30_phase_reviews/2026-11-05T0100/`. Their verifiers check input
bindings, exact certificates and independent fixed-order phase coverage.

The worker catalog now has 248,946 permuted signatures, 251,250 SAT clauses
and 72 count auxiliaries. Cold verification measured 17.7 seconds and CNF
construction 3.3 seconds. The old 2.889477999 C5 guide's entire quotient is
now excluded. The replacement guide scores 2.8521658155 and comes from pinned
chunk 55, record 45891; it passes current screens but is numerical only.
The separate `event30_strict_probability_lp.py` checks exact strict Farkas
identities, including zero-margin contradictions. A bounded trial on all 218
pre-phase residuals produced no new certificates. Keep this as research
infrastructure; no feasibility conclusion follows from its failure to certify.

A witness-guided asymmetric repair batch adds shared neighbors to break the
retired leader's disjoint-partition certificate, minimally clearing each new
triangle. It processes 223 previously unscored labels: 118 exact pre-score
exclusions and 105 scores, with 55 isomorphism classes within this batch and
50 explicit mappings. The best graph is fully asymmetric, theta 2.9092320745,
with exact individual upper bound 29103/10000. It is injected once as an
excluded mutation guide. The worker attempts a retained partition-witness
repair every 128 ordinary proposals, without changing other proposal weights,
screening or acceptance. Exact rejection replay and 32 deterministic repair
controls passed; natural chunk 88 execution confirmed valid proposals and
accepted moves. The source summary is `overnight30_asymmetric_reviews/2026-11-05T0100/`.
Raw repair and phase logs are immutable auxiliary archives in chunk 88,
`auxiliary_asymmetric_repair_0100/`, `auxiliary_asymmetric_repair_singleton_0100/`
and `auxiliary_phase_review_0100/`, published in experiment commit `434f719`.
The singleton archive adds 32 controls for a real 5+4+1 partition certificate.
Do not count auxiliary records as main-worker coverage.

Next priorities are the 213 remaining published quotient classes, new weakly
symmetric/asymmetric supports, and structural restrictions on nonuniform
dimension-19 realizations. Inspect the repair move's natural performance
before changing its cadence. Do not rerun the completed published-pool census.

## Operation and scheduled continuation

`overnight30_supervisor.py` runs bounded worker chunks, publishes checkpoints,
and enforces the absolute deadline. It has an exclusive lock and a status
heartbeat in the run directory. The workers are:

- `overnight30_symmetry.py`: weak-action SAT generation and full theta scoring.
- `overnight30_asymmetric.py`: adaptive asymmetric construction search.
- `overnight30_screen.py`: exact combinatorial and rational certificate checks.

The initial runner uses three slots, a 60-second first batch, then 600-second
chunks, with a 256 MiB archive ceiling. Completed JSONL pools are compressed.
The lifecycle, orphan-process cleanup and publication-order regression checks
are reproducible with `python search/test_overnight30_supervisor.py`.
For recovery in the original environment (after checking that no runner is
active), launch the following with detached process/session handling and
network access for the authorized Git pushes:

```sh
/workspace/threecontext-venv/bin/python search/overnight30_supervisor.py \
  --until 2026-11-05T05:00:00Z --data-repo /workspace/ghz30-overnight \
  --code-root /workspace/threecontext --python /workspace/threecontext-venv/bin/python \
  --workers 3 --first-chunk-seconds 60 --chunk-seconds 600 --max-archive-mb 256
```

The launch log is `/tmp/ghz30-supervisor.log`; the durable status and job history
are `status.json` and `manifest.json` inside the run directory. Each job records
its code commit and worker source hashes. A first-batch restart preserves
completed pools and advances chunk IDs rather than overwriting them.

Hourly scheduled continuations are research work, not merely status reports:

1. Fetch the two authorized branches and inspect the latest status, candidate
   pools, solver failures, coverage and best scores.
2. If the current supervisor is alive, let it continue. Do not launch a duplicate
   supervisor. In a different environment, a recent remote experiment checkpoint
   indicates an active publisher; inspect and analyze it instead of duplicating
   its work. If the heartbeat/checkpoints are stale and the original runner is
   unavailable, resume from the published experiment branch with new seeds and
   the same absolute deadline.
3. Work on the best new candidates and possible structural obstructions. Review
   and, where justified, improve search moves, filters or family allocation.
   Validate new filters against known constructions and distinguish a heuristic
   choice from a necessary condition. Push source/proof changes immediately.
   Useful next experiments include symmetry-preserving extensions and mutations
   of the best supports, and a separate dimension-19 mode learning pair-component
   cuts. Hall-deficiency cuts are global by Theorem 52; pair-disconnection cuts
   must not restrict the general event-number search. Initial best supports already had 95–122 edges;
   simply raising a minimum edge count did not reliably improve theta.
   The intersection helper has now generated a C2 seed with theta about 2.87972;
   the SAT worker loads it and mixes current-best neighborhood exploration with
   unrestricted models. On the next supervisor configuration change, skip the
   now-excluded `c5_mixed` family in selection while preserving the family-to-seed
   numbering and the other families' resume chains.
4. Keep original examples and exact witnesses. Store pools compactly, archive
   finished chunks when needed, and distinguish labelled graphs from proven
   nonisomorphic classes. A graph6 string is not a canonical isomorphism label.
5. Before the deadline, continue useful research without sending routine hourly
   notifications. Notify only on a credible near-GHZ candidate, an exact new
   result, a meaningful failure or required user action.
6. At or after 05:00 UTC, launch no new search. Check that workers have stopped,
   finish checkpoint publication, verify the main findings, merge experiment
   outputs into the research branch and push. Give the user a concise final
   report of coverage, best candidates, exact results and what remains open.

Do not turn numerical failures, interrupted runs, or an inability to find a
construction into nonexistence claims. Archive no chat unless the user requests
it. Do not send messages through email or other external services.

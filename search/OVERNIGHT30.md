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
at campaign setup. The latest numbered result is now Theorem 48. Results
39–48 cover contraction, orbit counting, moving cycles, large-orbit exclusions,
pair-component bridges, degree probability budgets, regular pair supports and
the excluded six-orbit quotient of the current C5 numerical leader, followed
by two-nine-orbit and affine-plane mixed-action exclusions.
Read Sections 1.10 and 2.1–2.3 before adding claims.

Already excluded and not worth repeating as standalone searches:

- The rotated-tail dimension-19 case (Theorem 35).
- Three pairwise-isoclinic reduced context spaces (Theorem 36).
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
  all excluded by rational theta duals. Two-kernel RTT remains open.

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
   every context. Both are now in the workers' exact screening/cut pipeline.
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

Next useful research steps, preserving the running supervisor:

- Add a verified quotient-class SAT cut for Proposition 46 (including context
  and orbit relabelings), stronger than excluding one support and its subgraphs.
- Consider C3×C3 mixed actions with orbit types `R=(9,1)` and
  `T=(3,3,3,1)`, using affine-plane translation actions on points and three
  parallel classes of lines. RRR and RRT are excluded, as is RTT with three
  distinct kernel lines. RTT with two shared kernel lines and multiplicities
  `(2,1)` in each T-context remains a useful direction. Necessary reductions
  and complete proofs are in `event30_c3square.py`.
- If changing runner configuration, skip `c5_mixed` without reassigning family
  seeds, and consider favoring C5-double after reviewing progress. Preserve
  diversity and the absolute deadline. Source commits capture dependency
  versions; expanding per-job hashes to the Fourier/screen/seed dependencies
  would improve provenance at the next runner restart.

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
   and Hall-deficiency cuts. Such dimension-specific cuts must not restrict the
   general event-number search. Initial best supports already had 95–122 edges;
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

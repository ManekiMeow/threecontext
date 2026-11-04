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
Keep Claude's theorem numbering in `RESULTS.md`; the last numbered result at
campaign setup is Theorem 38. Read Sections 1.10 and 2.1–2.3 before adding claims.

Already excluded and not worth repeating as standalone searches:

- The rotated-tail dimension-19 case (Theorem 35).
- Three pairwise-isoclinic reduced context spaces (Theorem 36).
- Regular unweighted Hoffman constructions on thirty events (Proposition 37).
- Identical copies of any linear 10_3 configuration (Theorem 38, exhaustive).
- Regular tri-Cayley families over groups of order at most twelve (earlier work).

Search directions, allowing nonuniform event probabilities and unequal overlap
weights unless a run explicitly records a narrower ansatz:

1. Weak C2, C3 and C5 actions with multiple orbits and fixed points, including
   different orbit decompositions in the three contexts. These are not the
   previously exhausted regular tri-Cayley actions.
2. Asymmetric tripartite triangle-free graphs: population search, neighbourhood
   replacements, multi-edge moves, crossovers and independent restarts. Start
   from the strongest existing graphs, not only the incidence constructions.
3. Apply exact combinatorial screens and mine surviving graph patterns for a
   new obstruction. In particular, a thirty-event GHZ graph must be connected,
   have no nonadjacent false twins, and every proper three-colouring must have
   classes (10,10,10). An unbalanced colouring is an exact exclusion witness.
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

## Operation and scheduled continuation

`overnight30_supervisor.py` runs bounded worker chunks, publishes checkpoints,
and enforces the absolute deadline. It has an exclusive lock and a status
heartbeat in the run directory. The workers are:

- `overnight30_symmetry.py`: weak-action SAT generation and full theta scoring.
- `overnight30_asymmetric.py`: adaptive asymmetric construction search.
- `overnight30_screen.py`: exact combinatorial and rational certificate checks.

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

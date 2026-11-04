# Realizability searches, night of 30–31 Oct 2026 (RESULTS.md §3.13)

Logs of `ks/realize_search.py` (flat-phase and general L-BFGS searches; DR/RRR on the control).
`SUMMARY.txt` lists starts, best, median and worst objective per run, produced by `search_summary.py`, which reads
the session's scratch directory. Two kinds of log:
* `control_m6`: the planted realizable control. It was solved by the flat search, and by the general search in
  17 of 40 starts.
* `cayley4_candidate_*`: the candidate skeletons. No run came near zero.

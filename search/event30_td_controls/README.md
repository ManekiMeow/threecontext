This directory preserves a DDD support graph that passes the structural and
local probability checks but is excluded by a rational theta dual.

Run `python search/event30_td_actions.py` from the repository root to regenerate
the graph in memory, compare it with `ddd_prism.json`, check the seven T/D action
edge-orbit counts, and replay `ddd_prism_certificate.json` exactly. The replay
uses the standard library and does not run an SDP solver.

The graph has 99 edges, three connected pair supports with perfect matchings,
distinct neighborhoods into each other context, and a length-three path for
each of the 135 same-context event pairs. The probabilities `21/200` on moving
vertices and `37/400` on fixed vertices satisfy all local inequalities strictly.
The exact dual bound is `6383/2500 < 3`, excluding this graph globally.

The seven types classify group actions, not support graphs. Neither this control
nor its certificate excludes the full DDD family. Proof notes and the explicit
construction are in `../event30_td_actions.py`.

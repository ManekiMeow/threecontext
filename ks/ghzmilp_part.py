# Same exact test as ghzmilp_run.py, for a range of normalisation rays r* (prints each status).
import sys, numpy as np
from ghzmilp_run import graph
from ghzstates import realisation, completions, symmetry_unitary, triangle_tensors, ghz_check
from ghzmilp import solve
name, mode, r0, r1, tl = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
P, cls, shift = graph(name)
X, psi0 = realisation(P, cls)
B = completions(X, cls, mode, np.random.default_rng(0), symmetry_unitary(X, shift))
E = triangle_tensors(B)
names = {0: "feasible", 1: "time limit", 2: "INFEASIBLE"}
for r in range(r0, r1):
    st, psi = solve(B, [psi0], r, tl)
    ok = ghz_check(psi, B, *E)[0] if psi is not None else None
    print(f"{name} {mode} r*={r}: {names.get(st, st)} {'(GHZ verified)' if ok else ''}", flush=True)

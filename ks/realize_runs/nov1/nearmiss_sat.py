"""Near misses of the CEGAR search: K4-free Cayley skeletons that satisfy every combinatorial condition of the SAT
encoding (DEG, C2, C4, TRI, quotient constraints) but have a neighbourhood with theta < 3, hence no realization.
For each SAT model, theta is computed at all four bases (the search stops at the first failure). Keeps the model
whose smallest local theta is largest.
usage: nearmiss.py P E SECONDS OUT.json      (multiplier family F(p, e))
"""
import os
import json
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))  # ks/
import warnings
warnings.filterwarnings("ignore")
import cayley4 as c4
import groups
from pysat.solvers import Solver

p, e, limit, out = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
G, orbits = groups.cyclic(p), c4.multiplier_orbits(p, e)
enc = c4.Cayley4(G, orbits=orbits)
enc.symmetry_breaking()
enc.quotient_constraints()
t0 = time.time()
best, models = None, 0
with Solver(name="cadical195", bootstrap_with=enc.cls) as solver:
    while time.time() - t0 < limit:
        solver.conf_budget(20000)
        res = solver.solve_limited()
        if res is None:
            continue
        if res is False:
            print("UNSAT (all models seen)", flush=True)
            break
        models += 1
        S = enc.sets(solver.get_model())
        ths, blocks = [], []
        for k in range(4):
            V, A = enc.neighbourhood(S, k)
            side = c4.two_colouring(A)
            th = 2.0 if side is not None else c4.theta_complement(V, A)
            ths.append(th)
            if th < 3 - 1e-6:
                blocks.append(enc.block_clause(S, k, V, A, side))
        key = (sum(t > 3 - 1e-6 for t in ths), min(ths))
        print(f"model {models}: sizes {[len(S[(a, b)]) for a, b in c4.PAIRS]} theta {np.round(ths, 4).tolist()}",
              flush=True)
        if best is None or key > best[0]:
            best = (key, S, ths)
            json.dump(dict(p=p, e=e, thetas=ths, sets={f"{a}{b}": S[(a, b)] for a, b in c4.PAIRS}),
                      open(out, "w"))
        if not blocks:
            print("*** passes all four theta tests ***", flush=True)
            blocks = [sorted({-enc.mem(a, b, s) if s in S[(a, b)] else enc.mem(a, b, s)
                              for a, b in c4.PAIRS for s in range(G.n)})]
        for cl in blocks:
            solver.add_clause(cl)
print(f"done: {models} models in {time.time() - t0:.0f} s; best {best[0] if best else None}", flush=True)

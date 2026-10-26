# Independent checker: rebuilds the graph from the model using only canonical (k < l) keys, no swapping code.
import sys, itertools, io, contextlib, importlib.util
d, NMIN, PMIN = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
sys.argv = ["skeleton_sat.py", str(d), str(NMIN), str(PMIN)]
spec = importlib.util.spec_from_file_location("sk", "skeleton_sat.py"); sk = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()): spec.loader.exec_module(sk)
if not sk.res: print("UNSAT"); raise SystemExit
model = set(v for v in sk.s.get_model() if v > 0); R = range(d)
adj = {(k, i): set() for k in range(4) for i in R}
for k, l in itertools.combinations(range(4), 2):
    for i, j in itertools.product(R, R):
        if sk.pool.obj2id[("e", k, i, l, j)] in model: adj[(k, i)].add((l, j)); adj[(l, j)].add((k, i))
V = list(adj)
k4 = sum(1 for i, j, m, n in itertools.product(R, repeat=4) if (1, j) in adj[(0, i)] and (2, m) in adj[(0, i)] and (3, n) in adj[(0, i)]
         and (2, m) in adj[(1, j)] and (3, n) in adj[(1, j)] and (3, n) in adj[(2, m)])
c2 = c4 = 0
for a, b in itertools.combinations(V, 2):
    for m in range(4):
        if m in (a[0], b[0]): continue
        common = sum(1 for c in R if (m, c) in adj[a] and (m, c) in adj[b])
        if b in adj[a] and common == 0: c2 += 1
        if b not in adj[a] and common == 1: c4 += 1
pb = min(sum(1 for x in adj[v] if x[0] == l) for v in V for l in range(4) if l != v[0])
dens = {(k, l): round(sum((l, j) in adj[(k, i)] for i in R for j in R) / d**2, 3) for k, l in itertools.combinations(range(4), 2)}
print(f"d={d}: K4 {k4}, C2 violations {c2}, C4 violations {c4}, min per-basis nbs {pb}, total degree {min(len(adj[v]) for v in V)}-{max(len(adj[v]) for v in V)}")
print("densities", dens)

"""Estimate the cost of a seeded SAT+nauty search by splitting it into cubes and timing them.

Cubes fix the seed-neighbourhoods of the first two new vertices (p+1, p+2) and the edge between
them. A new vertex adjacent to >= 2 seed rays is their cross product, so its seed-neighbourhood
is a subset of {seed rays orthogonal to c} for a ray c; every other neighbourhood is geometrically
impossible (the solver's orthogonality propagator refutes it at once). The feasible cubes are
therefore all products of: empty, a single seed ray, or a >= 2-subset of such a set.

usage: python3 estimate.py SATNAUTY_DIR SEED ORDER {all|SAMPLE_SIZE|idx:i,j,..} TIMEOUT_S [WORKERS]
"""
import itertools
import os
import random
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

from satnauty_seed import build, edge_var


def neighbourhoods(vecs):
    p = len(vecs)
    dot = lambda u, v: sum(a * b for a, b in zip(u, v))
    cross = lambda u, v: (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
    out = {frozenset()} | {frozenset([i]) for i in range(p)}
    for a, b in itertools.combinations(range(p), 2):
        c = cross(vecs[a], vecs[b])
        if not any(c):
            continue
        full = [i for i in range(p) if dot(vecs[i], c) == 0]
        for r in range(2, len(full) + 1):
            out |= {frozenset(s) for s in itertools.combinations(full, r)}
    return sorted(out, key=lambda s: (len(s), sorted(s)))


def main():
    root, seed, order, mode, tlim = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], float(sys.argv[5])
    workers = int(sys.argv[6]) if len(sys.argv) > 6 else 4
    wd, cnf, vecs, cmd = build(root, seed, order)
    p = len(vecs)
    N = neighbourhoods(vecs)
    cubes = []
    for n1, n2, e in itertools.product(N, N, (0, 1)):
        lits = [(edge_var(i + 1, p + 1) if i in n1 else -edge_var(i + 1, p + 1)) for i in range(p)]
        lits += [(edge_var(i + 1, p + 2) if i in n2 else -edge_var(i + 1, p + 2)) for i in range(p)]
        lits.append(edge_var(p + 1, p + 2) if e else -edge_var(p + 1, p + 2))
        cubes.append(lits)
    print(f"{seed} order {order}: {len(N)} neighbourhoods, {len(cubes)} cubes", flush=True)
    random.seed(1)
    if mode == "all":
        idx = list(range(len(cubes)))
    elif mode.startswith("idx:"):
        idx = [int(x) for x in mode[4:].split(",")]
    else:
        idx = random.sample(range(len(cubes)), int(mode))
    base = open(cnf).read().splitlines()
    nv, nc = map(int, base[0].split()[2:])
    cdir = os.path.join(wd, "cubes")
    os.makedirs(cdir, exist_ok=True)

    def run(k):
        f = os.path.join(cdir, f"c{k}.cnf")
        with open(f, "w") as fh:
            fh.write(f"p cnf {nv} {nc + len(cubes[k])}\n")
            fh.write("\n".join(base[1:]) + "\n")
            fh.write("".join(f"{l} 0\n" for l in cubes[k]))
        t0 = time.time()
        try:
            r = subprocess.run([cmd[0], f] + cmd[2:], capture_output=True, text=True, timeout=tlim, cwd=wd)
            dt, done = time.time() - t0, r.returncode == 20
            sols = [l for l in r.stdout.splitlines() if l.startswith("Number of solutions")]
            if not sols or not sols[0].rstrip().endswith(": 0"):
                open(os.path.join(cdir, f"c{k}.log"), "w").write(r.stdout)
        except subprocess.TimeoutExpired:
            dt, done, sols = tlim, False, []
        os.remove(f)
        return k, dt, done, sols

    t0 = time.time()
    res = []
    with ThreadPoolExecutor(workers) as ex:
        for k, dt, done, sols in ex.map(run, idx):
            res.append((dt, done))
            print(f"cube {k}: {'done' if done else 'TIMEOUT'} {dt:.2f}s {sols[0] if sols else ''}", flush=True)
    tot = sum(d for d, _ in res)
    nto = sum(1 for _, d in res if not d)
    scale = len(cubes) / len(idx)
    print(f"sampled {len(idx)}/{len(cubes)} cubes, {nto} timeouts, sum {tot:.1f}s, wall {time.time()-t0:.1f}s; "
          f"estimated total {'>= ' if nto else ''}{tot*scale:.0f} CPU-s", flush=True)


if __name__ == "__main__":
    main()

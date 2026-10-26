"""Run the SAT+nauty solver (github.com/BrianLi009/SAT-nauty) on a KS search seeded by an
arbitrary fixed ray set, e.g. the 13-ray Yu-Oh set.

usage: python3 satnauty_seed.py SATNAUTY_DIR SEED ORDER
SEED in {yuoh13, yuoh13x, sic25}. yuoh13x = Yu-Oh seed with the full 25-ray completion forbidden:
at least one of the 12 rays u x w (u, w non-orthogonal Yu-Oh rays) must be absent, encoded as
"no new vertex is adjacent to both u and w" behind a selector per completion ray. The clauses are
invariant under relabelling the non-seed vertices, so orderly generation stays sound. Writes seed files and the CNF under SATNAUTY_DIR/seed_<SEED>_<ORDER>/.
"""
import os
import subprocess
import sys
import time

SEEDS = {
    # Yu-Oh 2012: 3 axes, 6 face diagonals, 4 body diagonals
    "yuoh13": [(1, 0, 0), (0, 1, 0), (0, 0, 1),
               (0, 1, 1), (0, 1, -1), (1, 0, 1), (1, 0, -1), (1, 1, 0), (1, -1, 0),
               (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)],
}


def canon(v):
    from math import gcd
    g = 0
    for x in v:
        g = gcd(g, abs(x))
    v = tuple(x // g for x in v)
    return v if next(x for x in v if x) > 0 else tuple(-x for x in v)


def edge_var(i, j):
    """1-based vertices i < j; edge variables are numbered column by column."""
    return (j - 1) * (j - 2) // 2 + i


def forbid_completion(cnf, vecs, order):
    p = len(vecs)
    cross = lambda u, v: (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
    seedrays = {canon(v) for v in vecs}
    gen = {}
    for a in range(p):
        for b in range(a + 1, p):
            c = cross(vecs[a], vecs[b])
            if any(c) and canon(c) not in seedrays:
                gen.setdefault(canon(c), (a + 1, b + 1))
    assert len(gen) == 12, len(gen)
    lines = open(cnf).read().splitlines()
    _, _, nv, nc = lines[0].split()
    nv, nc = int(nv), int(nc)
    new = []
    sel = []
    for k, (c, (u, w)) in enumerate(sorted(gen.items())):
        a = nv + 1 + k
        sel.append(a)
        for v in range(p + 1, order + 1):
            new.append(f"{-a} {-edge_var(u, v)} {-edge_var(w, v)} 0")
    new.append(" ".join(map(str, sel)) + " 0")
    lines[0] = f"p cnf {nv + len(sel)} {nc + len(new)}"
    open(cnf, "w").write("\n".join(lines + new) + "\n")


def build(root, seed, order):
    """Write seed files and the seeded CNF; return (workdir, cnf, seed vectors, solver command)."""
    root = os.path.abspath(root)
    wd = os.path.join(root, f"seed_{seed}_{order}")
    os.makedirs(wd, exist_ok=True)
    if seed == "sic25":
        vecs = [tuple(map(int, l.split())) for l in open(os.path.join(root, "sic-25-vectors.txt")) if l.strip()]
    else:
        vecs = SEEDS[seed.rstrip("x")]
    p = len(vecs)
    dot = lambda u, v: sum(a * b for a, b in zip(u, v))
    # The solver checks every prefix, the seed itself included, for RCL canonicity, so the seed
    # must be listed in its own RCL canonical order; relabel with the verifier's rcl_canon.
    canon = os.path.join(root, "verifiers/rcl_canon")
    for _ in range(10):
        bits = "".join("1" if dot(vecs[i], vecs[j]) == 0 else "0" for j in range(p) for i in range(j))
        out = subprocess.run([canon], input=f"{p} {p} {bits}\n", capture_output=True, text=True).stdout.split()
        if out[1] == "1":
            break
        pos = list(map(int, out[2:]))
        new = [None] * p
        for i, q in enumerate(pos):
            new[q] = vecs[i]
        vecs = new
    else:
        raise RuntimeError("seed did not reach RCL canonical form")
    # edge variables are numbered column by column: (1,2), (1,3), (2,3), (1,4), ...
    lits, var = [], 0
    for j in range(p):
        for i in range(j):
            var += 1
            lits.append(var if dot(vecs[i], vecs[j]) == 0 else -var)
    open(os.path.join(wd, "seed.vars"), "w").write(" ".join(map(str, lits)) + "\n")
    open(os.path.join(wd, "seed-vectors.txt"), "w").write("\n".join(" ".join(map(str, v)) for v in vecs) + "\n")
    base = os.path.join(wd, f"constraints_{order}_0_no_lex")
    if not os.path.exists(base):
        link = os.path.join(wd, "gen_instance")
        if not os.path.exists(link):
            os.symlink(os.path.join(root, "gen_instance"), link)
        subprocess.run(["python3", "gen_instance/generate.py", str(order), "0", "no-lex"],
                       cwd=wd, check=True, capture_output=True)
    cnf = os.path.join(wd, "seeded.cnf")
    subprocess.run(["python3", os.path.join(root, "utils/add_vars_to_cnf.py"), base,
                    os.path.join(wd, "seed.vars"), cnf], check=True)
    if seed.endswith("x"):
        forbid_completion(cnf, vecs, order)
    cmd = [os.path.join(root, "cadical-rcl/build/cadical"), cnf, "--order", str(order),
           "--partition", str(p), "--vectors-file", os.path.join(wd, "seed-vectors.txt"),
           "--binary=false", "--unembeddable-check", "0", "--ortho"]
    return wd, cnf, vecs, cmd


def main():
    root, seed, order = sys.argv[1], sys.argv[2], int(sys.argv[3])
    wd, cnf, vecs, cmd = build(root, seed, order)
    t0 = time.time()
    with open(os.path.join(wd, "solver.log"), "w") as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, cwd=wd)
    dt = time.time() - t0
    out = open(os.path.join(wd, "solver.log")).read()
    keep = [l for l in out.splitlines() if l.startswith(("Number of solutions", "Colourable models",
                                                          "Canonical subgraphs", "Solution"))]
    print(f"seed={seed} order={order} exit={r.returncode} wall={dt:.1f}s")
    print("\n".join(keep[:20]), flush=True)


if __name__ == "__main__":
    main()

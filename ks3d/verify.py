"""Independent check of a 3D KS set given as integer vectors: recompute orthogonality exactly,
check uncolourability with a second SAT solver (Glucose), and check vertex-criticality."""
import json
import sys

from pysat.solvers import Glucose4


def colourable(vs):
    n = len(vs)
    d = lambda a, b: sum(x * y for x, y in zip(vs[a], vs[b]))
    E = [(i, j) for i in range(n) for j in range(i + 1, n) if d(i, j) == 0]
    T = [(i, j, k) for i, j in E for k in range(j + 1, n) if d(i, k) == 0 and d(j, k) == 0]
    cl = [[-(i + 1), -(j + 1)] for i, j in E] + [[i + 1, j + 1, k + 1] for i, j, k in T]
    with Glucose4(bootstrap_with=cl) as s:
        return s.solve(), len(E), len(T)


def check(vs):
    ok, ne, nt = colourable(vs)
    crit = all(colourable(vs[:i] + vs[i + 1:])[0] for i in range(len(vs)))
    return {"rays": len(vs), "orthogonal_pairs": ne, "bases": nt,
            "uncolourable": not ok, "vertex_critical": crit}


if __name__ == "__main__":
    for name, vs in json.load(open(sys.argv[1])).items():
        print(name, check([tuple(v) for v in vs]))

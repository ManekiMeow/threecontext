"""Heuristic construction search with independently labelled 10_3 line systems.

The point context has ten events. B and C are independently chosen linear 10_3
configurations on those points; edges between B,C join disjoint lines. Unlike
construction30_configurations.py, this search is not exhaustive. It optimizes
the relative point labelling of C and occasionally changes its incidence type.
The Lovasz SDP permits arbitrary event probabilities and edge overlaps.

  OPENBLAS_NUM_THREADS=1 python search/construction30_mixed.py search
  python search/construction30_mixed.py verify

Verification uses exact rational dual certificates saved for each run's best
graph, and does not need an SDP solver. It excludes those graphs only.
"""
import argparse
from fractions import Fraction as Q
from importlib.metadata import version
import itertools
import json
import math
from pathlib import Path
import random
import time

from construction30_configurations import positive_definite


def validate_lines(lines):
    lines = [frozenset(line) for line in lines]
    assert len(lines) == len(set(lines)) == 10
    assert all(len(line) == 3 and line <= set(range(10)) for line in lines)
    assert all(sum(p in line for line in lines) == 3 for p in range(10))
    assert all(len(a & b) <= 1 for a, b in itertools.combinations(lines, 2))
    return lines


def support(B, C):
    B, C = validate_lines(B), validate_lines(C)
    adj = [set() for _ in range(30)]

    def edge(a, b):
        adj[a].add(b)
        adj[b].add(a)

    for offset, lines in ((10, B), (20, C)):
        for j, line in enumerate(lines):
            for p in line:
                edge(p, offset+j)
    for j, k in itertools.product(range(10), repeat=2):
        if B[j].isdisjoint(C[k]):
            edge(10+j, 20+k)
    assert all(len(adj[p]) == 6 for p in range(10))
    assert all(1 <= len(adj[10+j] & set(range(20, 30))) <= 4 for j in range(10))
    assert all(v//10 != u//10 and not (adj[v] & adj[u]) for u in range(30) for v in adj[u])
    return adj


def as_graph(B, C):
    import networkx as nx
    adj = support(B, C)
    graph = nx.Graph()
    graph.add_nodes_from(range(30))
    graph.add_edges_from((u, v) for u in range(30) for v in adj[u] if u < v)
    return graph


def labelled(lines, permutation):
    return tuple(sorted(tuple(sorted(permutation[p] for p in line)) for line in lines))


def verify_run(run):
    adj = support(run["B_lines"], run["C_lines"])
    cert = run["dual_certificate"]
    bound, denominator = Q(cert["bound"]), cert["denominator"]
    assert bound < 3
    free = iter(cert["free_numerators"])
    dual = [[Q(1) for _ in range(30)] for _ in range(30)]
    for u in range(30):
        for v in range(u+1, 30):
            if v not in adj[u]:
                dual[u][v] = dual[v][u] = Q(next(free), denominator)
    sentinel = object()
    assert next(free, sentinel) is sentinel
    assert positive_definite([[(bound if u == v else 0)-dual[u][v] for v in range(30)]
                              for u in range(30)])


def certify_run(run):
    from dualcert import dual_upper
    graph = as_graph(run["B_lines"], run["C_lines"])
    _, _, matrix, numeric = dual_upper(graph)
    denominator = 10000
    # One hundredth slack safely exceeds entry rounding errors. Exact elimination
    # checks the result; no spectral-tolerance inference is used as a proof.
    bound = Q(math.ceil((numeric+0.01)*1000), 1000)
    if bound >= 3:
        raise RuntimeError("near-GHZ candidate needs a primal certificate, not exclusion")
    run["dual_certificate"] = {
        "bound": str(bound), "denominator": denominator,
        "free_numerators": [round(float(matrix[u, v])*denominator)
                            for u in range(30) for v in range(u+1, 30) if not graph.has_edge(u, v)]}
    verify_run(run)


def search(path, steps, seeds):
    from construction30_configurations import enumerate_configurations
    from theta import theta_comp
    _, configurations = enumerate_configurations()
    start = time.monotonic()
    results = {"scope": "heuristic mixed line-system search; not an exhaustive exclusion",
               "packages": {p: version(p) for p in ("numpy", "networkx", "cvxpy", "clarabel")},
               "steps_per_run": steps, "seeds_per_B_type": seeds, "runs": []}
    for btype, B in enumerate(configurations):
        cache = {}

        def score(C):
            if C not in cache:
                cache[C] = float(theta_comp(as_graph(B, C))[0])
            return cache[C]

        for seed in range(seeds):
            rng = random.Random(20261104 + 1000*btype + seed)
            ctype = btype if seed == 0 else rng.randrange(len(configurations))
            perm = list(range(10))
            if seed:
                rng.shuffle(perm)
            C = labelled(configurations[ctype], perm)
            current = score(C)
            best, best_C, best_type = current, C, ctype
            for step in range(steps):
                trial_type, trial_perm = ctype, perm[:]
                if rng.random() < 0.12:
                    trial_type = rng.randrange(len(configurations))
                else:
                    a, b = rng.sample(range(10), 2)
                    trial_perm[a], trial_perm[b] = trial_perm[b], trial_perm[a]
                trial_C = labelled(configurations[trial_type], trial_perm)
                value = score(trial_C)
                temp = 0.025*(1-step/max(steps, 1)) + 0.0005
                if value >= current or rng.random() < math.exp((value-current)/temp):
                    ctype, perm, C, current = trial_type, trial_perm, trial_C, value
                if value > best:
                    best, best_C, best_type = value, trial_C, trial_type
                if best > 3-1e-6:
                    path.write_text(json.dumps({"candidate": {"B_lines": B, "C_lines": best_C,
                                                              "theta_numeric": best}}, indent=2)+"\n")
                    raise RuntimeError("near-GHZ candidate found; inspect saved candidate immediately")
            run = {"B_type": btype, "C_type": best_type, "seed": seed,
                   "best_theta_numeric": best, "B_lines": B, "C_lines": best_C}
            certify_run(run)
            results["runs"].append(run)
            path.write_text(json.dumps(results, indent=2)+"\n")
            print("B type", btype, "seed", seed, "best", round(best, 9),
                  "exact upper", run["dual_certificate"]["bound"],
                  "elapsed", round(time.monotonic()-start, 1), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("search", "verify"))
    parser.add_argument("--steps", type=int, default=70)
    parser.add_argument("--seeds", type=int, default=2)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("construction30_mixed_results.json"))
    args = parser.parse_args()
    if args.mode == "search":
        if args.steps < 0 or args.seeds < 1:
            parser.error("steps must be nonnegative and seeds positive")
        search(args.output, args.steps, args.seeds)
    else:
        results = json.loads(args.output.read_text())
        for run in results["runs"]:
            verify_run(run)
        print(len(results["runs"]), "saved best graphs exactly excluded; mixed family remains open")

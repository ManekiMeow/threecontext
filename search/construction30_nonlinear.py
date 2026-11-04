"""Bounded construction search allowing repeated point pairs in the line systems.

Each B/C system has ten distinct three-point lines, with every one of the ten
points on three lines. Two lines may share TWO points, unlike a linear 10_3
configuration. Independent degree-preserving incidence switches mutate B and C;
their cross edges join disjoint lines, so triangle-freeness is automatic.
Duplicate lines within one context would be nonadjacent false twins in H.
Merging their primal Gram vectors preserves trace and objective, giving a
29-event theta=3 graph, which the existing lower bound excludes. Thus requiring
distinct lines loses no 30-event GHZ example in this construction family.

Start from the best saved mixed-linear candidate. This is a heuristic search,
not an exhaustive exclusion of nonlinear configurations. Exact duals exclude
only the saved best graph from each run.

  OPENBLAS_NUM_THREADS=1 python search/construction30_nonlinear.py search
  python search/construction30_nonlinear.py verify
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


def validate(lines):
    lines = [frozenset(line) for line in lines]
    assert len(lines) == len(set(lines)) == 10
    assert all(len(line) == 3 and line <= set(range(10)) for line in lines)
    assert all(sum(p in line for line in lines) == 3 for p in range(10))
    return lines


def canonical(lines):
    return tuple(sorted(tuple(sorted(line)) for line in lines))


def support(B, C):
    B, C = validate(B), validate(C)
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
    assert all(v//10 != u//10 and not (adj[v] & adj[u]) for u in range(30) for v in adj[u])
    return adj


def graph(B, C):
    import networkx as nx
    adj = support(B, C)
    H = nx.Graph()
    H.add_nodes_from(range(30))
    H.add_edges_from((u, v) for u in range(30) for v in adj[u] if u < v)
    return H


def mutate(lines, rng):
    """A bipartite 2-switch preserves all line sizes and all point degrees."""
    for _ in range(50):
        i, j = rng.sample(range(10), 2)
        a, b = set(lines[i]), set(lines[j])
        p, q = rng.choice(sorted(a-b)), rng.choice(sorted(b-a))
        a.remove(p)
        a.add(q)
        b.remove(q)
        b.add(p)
        trial = list(lines)
        trial[i], trial[j] = tuple(sorted(a)), tuple(sorted(b))
        trial = canonical(trial)
        if trial != lines and len(set(trial)) == 10:
            return trial
    return None


def repetitions(lines):
    counts = [sum(set(pair) <= set(line) for line in lines)
              for pair in itertools.combinations(range(10), 2)]
    return {"repeated_point_pairs": sum(c > 1 for c in counts), "maximum_pair_multiplicity": max(counts)}


def verify_run(run):
    adj = support(run["B_lines"], run["C_lines"])
    assert run["B_repetitions"] == repetitions(run["B_lines"])
    assert run["C_repetitions"] == repetitions(run["C_lines"])
    cert = run["dual_certificate"]
    bound, den = Q(cert["bound"]), cert["denominator"]
    assert bound < 3
    free = iter(cert["free_numerators"])
    dual = [[Q(1) for _ in range(30)] for _ in range(30)]
    for u in range(30):
        for v in range(u+1, 30):
            if v not in adj[u]:
                dual[u][v] = dual[v][u] = Q(next(free), den)
    marker = object()
    assert next(free, marker) is marker
    assert positive_definite([[(bound if u == v else 0)-dual[u][v] for v in range(30)]
                              for u in range(30)])


def certify_run(run):
    from dualcert import dual_upper
    H = graph(run["B_lines"], run["C_lines"])
    _, _, matrix, optimum = dual_upper(H)
    bound, den = Q(math.ceil((optimum+0.01)*1000), 1000), 10000
    if bound >= 3:
        raise RuntimeError("near-GHZ candidate requires a primal certificate")
    run["dual_certificate"] = {
        "bound": str(bound), "denominator": den,
        "free_numerators": [round(float(matrix[u, v])*den)
                            for u in range(30) for v in range(u+1, 30) if not H.has_edge(u, v)]}
    verify_run(run)


def search(output, start_path, steps, seeds):
    from theta import theta_comp
    source = max(json.loads(start_path.read_text())["runs"], key=lambda x: x["best_theta_numeric"])
    initial_B, initial_C = canonical(source["B_lines"]), canonical(source["C_lines"])
    results = {"scope": "heuristic nonlinear degree-three incidence search, not exhaustive",
               "source": start_path.name, "source_best_theta_numeric": source["best_theta_numeric"],
               "packages": {p: version(p) for p in ("numpy", "networkx", "cvxpy", "clarabel")},
               "steps_per_run": steps, "seeds": seeds, "runs": []}
    start = time.monotonic()
    cache = {}

    def score(B, C):
        key = tuple(sorted((B, C)))  # exchanging the two line contexts is an isomorphism
        if key not in cache:
            cache[key] = float(theta_comp(graph(B, C))[0])
        return cache[key]

    for seed in range(seeds):
        rng = random.Random(20261104+seed)
        B, C = initial_B, initial_C
        current = score(B, C)
        best, best_B, best_C = current, B, C
        proposed = accepted = failed = 0
        cache_start = len(cache)
        for step in range(steps):
            side = rng.randrange(2)
            trial = mutate(B if side == 0 else C, rng)
            if trial is None:
                failed += 1
                continue
            trial_B, trial_C = (trial, C) if side == 0 else (B, trial)
            value = score(trial_B, trial_C)
            proposed += 1
            temp = 0.02*(1-step/max(steps, 1)) + 0.0005
            if value >= current or rng.random() < math.exp((value-current)/temp):
                B, C, current = trial_B, trial_C, value
                accepted += 1
            if value > best:
                best, best_B, best_C = value, trial_B, trial_C
            if best > 3-1e-6:
                output.write_text(json.dumps({"candidate": {"B_lines": best_B, "C_lines": best_C,
                                                            "theta_numeric": best}}, indent=2)+"\n")
                raise RuntimeError("near-GHZ candidate found; inspect saved candidate immediately")
        run = {"seed": seed, "proposed_mutations": proposed, "accepted_mutations": accepted,
               "failed_mutations": failed, "new_SDP_solves": len(cache)-cache_start,
               "best_theta_numeric": best, "B_lines": best_B, "C_lines": best_C,
               "B_repetitions": repetitions(best_B), "C_repetitions": repetitions(best_C)}
        certify_run(run)
        results["runs"].append(run)
        output.write_text(json.dumps(results, indent=2)+"\n")
        print("seed", seed, "best", round(best, 9), "exact upper", run["dual_certificate"]["bound"],
              "proposals", proposed, "repeated pairs", run["B_repetitions"]["repeated_point_pairs"],
              run["C_repetitions"]["repeated_point_pairs"], "elapsed", round(time.monotonic()-start, 1), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("search", "verify"))
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--seeds", type=int, default=4)
    parser.add_argument("--start", type=Path,
                        default=Path(__file__).with_name("construction30_mixed_results.json"))
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("construction30_nonlinear_results.json"))
    args = parser.parse_args()
    if args.mode == "search":
        if args.steps < 0 or args.seeds < 1:
            parser.error("steps must be nonnegative and seeds positive")
        search(args.output, args.start, args.steps, args.seeds)
    else:
        results = json.loads(args.output.read_text())
        for run in results["runs"]:
            verify_run(run)
        print(len(results["runs"]), "saved best graphs exactly excluded; nonlinear family remains open")

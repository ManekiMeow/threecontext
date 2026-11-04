"""Bounded, reproducible search in the four remaining event-size partitions.

This is a heuristic, never an exhaustive nonexistence test. Start from H33,
move one line into the point context, delete surplus vertices, then anneal
over triangle-free graphs while keeping the new contexts fixed.
"""
import argparse
from importlib.metadata import version
import json
import math
from pathlib import Path
import platform
import random
import time

import networkx as nx
from anneal import maximize
from theta import theta_comp


def initial(shape, rng):
    source = nx.read_graph6(Path(__file__).with_name("h33.g6"))
    groups = [list(range(9)), list(range(9, 21)), list(range(21, 33))]
    moved = rng.choice(groups[1])
    groups[1].remove(moved)
    groups[0].append(moved)
    groups = [rng.sample(g, size) for g, size in zip(groups, shape)]
    col = {v: k for k, g in enumerate(groups) for v in g}
    H = source.subgraph(col).copy()
    H.remove_edges_from([(u, v) for u, v in H.edges() if col[u] == col[v]])
    return maximize(H, col, rng), col


def search(shape, seed, steps):
    rng = random.Random(seed)
    H, col = initial(shape, rng)
    value = theta_comp(H)[0]
    best, graph = value, H.copy()
    start = value
    for step in range(steps):
        trial = H.copy()
        if step % 7 == 0:
            trial.remove_edges_from(list(trial.edges(rng.choice(list(trial)))))
        else:
            edges = list(trial.edges())
            trial.remove_edges_from(rng.sample(edges, min(len(edges), 3)))
        maximize(trial, col, rng)
        candidate = theta_comp(trial)[0]
        temp = 0.012*(1-step/max(steps, 1)) + 0.0001
        if candidate >= value or rng.random() < math.exp((candidate-value)/temp):
            H, value = trial, candidate
        if candidate > best:
            best, graph = candidate, trial.copy()
    assert all(col[u] != col[v] for u, v in graph.edges())
    assert not any(nx.triangles(graph).values())
    nodes = sorted(graph)
    idx = {v: i for i, v in enumerate(nodes)}
    return {"shape": shape, "seed": seed, "steps": steps,
            "initial_theta_numeric": float(start), "best_theta_numeric": float(best),
            "parts": [[idx[v] for v in nodes if col[v] == k] for k in range(3)],
            "edges": sorted([sorted((idx[u], idx[v])) for u, v in graph.edges()])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steps", type=int, default=150)
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--verify", action="store_true",
                        help="generate and exactly check dual bounds for the best saved graph per shape")
    args = parser.parse_args()
    if args.verify:
        from dualcert import exact_upper
        runs = json.loads(args.output.read_text())["runs"]
        shapes = sorted({tuple(run["shape"]) for run in runs})
        for shape in shapes:
            run = max((r for r in runs if tuple(r["shape"]) == shape),
                      key=lambda r: r["best_theta_numeric"])
            parts = run["parts"]
            assert tuple(map(len, parts)) == shape
            col = {v: k for k, part in enumerate(parts) for v in part}
            assert sorted(col) == list(range(sum(shape)))
            H = nx.Graph()
            H.add_nodes_from(col)
            H.add_edges_from(run["edges"])
            assert set(H) == set(col)
            assert all(col[u] != col[v] for u, v in H.edges())
            assert not any(nx.triangles(H).values())
            upper, numeric = exact_upper(H)
            if upper is None or upper >= 3:
                raise RuntimeError(f"no exact exclusion certified for {shape}")
            print(shape, "best saved seed", run["seed"], "exact theta upper", upper,
                  "numeric dual", numeric, flush=True)
        return
    if args.steps < 0 or args.seeds < 1:
        parser.error("steps must be nonnegative and seeds must be positive")
    results = []
    environment = {"python": platform.python_version(),
                   "packages": {name: version(name) for name in ("numpy", "networkx", "cvxpy", "clarabel")}}
    start = time.monotonic()
    for shape in ((10, 10, 10), (10, 10, 11), (10, 10, 12), (10, 11, 11)):
        for seed in range(args.seeds):
            result = search(shape, seed, args.steps)
            results.append(result)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps({
                "scope": "heuristic; no exhaustive exclusion or realization certificate",
                "environment": environment,
                "runs": results}, indent=2) + "\n")
            print(shape, "seed", seed, "best theta", result["best_theta_numeric"],
                  "elapsed", round(time.monotonic()-start, 1), flush=True)


if __name__ == "__main__":
    main()

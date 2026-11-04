"""Seed weak C2 supports from high-theta asymmetric graphs.

Choose one transposed pair per context, intersect a source with its image, and
greedily add complete C2 edge orbits while preserving triangle-freeness and the
necessary degree bound ten. This is a construction heuristic, never an
exhaustive exclusion. All event weights remain free in the full theta SDP.

Example:
  python search/overnight30_symmetry_intersection.py --source best.json \
      --seconds 180 --seed 1 --output-dir intersection_run

Sources must have three contiguous contexts of ten vertices. Every compact
candidate has graph6 and provenance; best/near-3 files also contain edges,
parts, numerical diagnostics and the exact structural-screen result.
"""
import argparse
from collections import Counter
import heapq
import itertools
import json
from pathlib import Path
import random
import time

from overnight30_symmetry import (ThetaProcess, atomic_json, build, check_graph,
                                  generator, records)
from overnight30_screen import normalize_graph, screen_graph


def intersection(edges, pairs):
    """Return the invariant intersection, relabelled to standard C2_onepair."""
    permutation = list(range(30))
    relabel = {}
    for k, (u, v) in enumerate(pairs):
        permutation[u], permutation[v] = v, u
        order = [u, v]+[w for w in range(10*k, 10*k+10) if w not in (u, v)]
        relabel.update({w: 10*k+i for i, w in enumerate(order)})
    transformed = {tuple(sorted((permutation[u], permutation[v]))) for u, v in edges}
    return {tuple(sorted((relabel[u], relabel[v]))) for u, v in edges & transformed}


def ranked_pairs(edges, deadline, keep=256):
    """Rank all 45^3 transposition triples by retained edges, budget permitting."""
    pairs = [list(itertools.combinations(range(10*k, 10*k+10), 2)) for k in range(3)]
    heap, examined = [], 0
    for triple in itertools.product(*pairs):
        if examined % 256 == 0 and time.monotonic() >= deadline:
            break
        permutation = list(range(30))
        for u, v in triple:
            permutation[u], permutation[v] = v, u
        count = sum(tuple(sorted((permutation[u], permutation[v]))) in edges for u, v in edges)
        item = count, triple
        if len(heap) < keep:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)
        examined += 1
    return sorted(heap, reverse=True), examined


def maximal_extension(edges, orbits, rng, preferred):
    """One randomized pass is maximal: additional edges only obstruct later additions."""
    adj = [set() for _ in range(30)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    if max(map(len, adj)) > 10:
        return None
    order = list(orbits)
    rng.shuffle(order)
    # Alternate source-guided and completely random extension orders.
    if preferred is not None:
        order.sort(key=lambda orbit: -sum(edge in preferred for edge in orbit))
    added = set(edges)
    for orbit in order:
        if orbit[0] in added:
            continue
        inserted, valid = [], True
        for u, v in orbit:
            if len(adj[u]) >= 10 or len(adj[v]) >= 10 or adj[u] & adj[v]:
                valid = False
                break
            adj[u].add(v)
            adj[v].add(u)
            inserted.append((u, v))
        if valid:
            added.update(inserted)
        else:
            for u, v in inserted:
                adj[u].remove(v)
                adj[v].remove(u)
    return added


def run(args):
    import networkx as nx
    start = time.monotonic()
    deadline = start+args.seconds
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    state_path, pool_path = output/"state.json", output/"models.jsonl"
    old = json.loads(state_path.read_text()) if state_path.exists() else {}
    if old and old["seed"] != args.seed:
        raise ValueError("resume seed differs")
    state = dict(seed=args.seed, batches=old.get("batches", 0)+1,
                 family="c2_onepair", scope="heuristic source-intersection and maximal-orbit-extension construction",
                 sources=[str(p.resolve()) for p in args.source], context_sizes=[10, 10, 10],
                 pair_triples_ranked=0, generated=0, scored=0, screened_out=0,
                 duplicate_labelled=0, solver_errors=0)
    seen, statuses = set(), Counter()
    for record in records(pool_path):
        seen.add(record["graph6"])
        statuses[record["status"]] += 1
    state.update(generated=sum(statuses.values()), scored=statuses["scored_numerically"],
                 screened_out=statuses["excluded_structurally"], solver_errors=statuses["solver_error"])
    rng = random.Random(args.seed+1000003*state["batches"])
    solver, _, orbits, permutation = build("c2_onepair")
    solver.delete()
    sources = []
    for path in args.source:
        data = json.loads(path.read_text())
        normalized = normalize_graph(data)
        if normalized["n"] != 30 or normalized["parts"] != [list(range(10*k, 10*k+10)) for k in range(3)]:
            raise ValueError(f"source does not use three contiguous ten-event contexts: {path}")
        edges = {tuple(edge) for edge in normalized["edges"]}
        graph = nx.Graph()
        graph.add_nodes_from(range(30))
        graph.add_edges_from(edges)
        if any(nx.triangles(graph).values()) or any(u//10 == v//10 for u, v in edges):
            raise ValueError(f"source is not triangle-free and context-partite: {path}")
        ranked, count = ranked_pairs(edges, deadline, args.ranked_pairs)
        state["pair_triples_ranked"] += count
        if ranked:
            sources.append((str(path.resolve()), edges, ranked))
            print(path.name, "ranked", count, "best retained edges", ranked[0][0], flush=True)
    best_path = output/"best.json"
    best = json.loads(best_path.read_text()) if best_path.exists() else None
    near = output/"near3"
    near.mkdir(exist_ok=True)
    process = ThetaProcess()
    batches_generated = state["generated"]

    def checkpoint():
        state.update(elapsed_seconds=time.monotonic()-start, updated_unix=time.time(),
                     distinct_labelled_supports=len(seen), best_theta_numeric=best["theta_numeric"] if best else None,
                     batch_generated=state["generated"]-batches_generated)
        atomic_json(state_path, state)

    with pool_path.open("a", buffering=1) as stream:
        def evaluate(edges, provenance):
            nonlocal best
            graph = nx.Graph()
            graph.add_nodes_from(range(30))
            graph.add_edges_from(sorted(edges))
            assert not any(nx.triangles(graph).values())
            assert all(graph.has_edge(permutation[u], permutation[v]) for u, v in graph.edges())
            g6 = nx.to_graph6_bytes(graph, header=False).decode().strip()
            if g6 in seen:
                state["duplicate_labelled"] += 1
                return
            seen.add(g6)
            state["generated"] += 1
            record = dict(graph6=g6, family="c2_onepair", seed=args.seed, context_sizes=[10, 10, 10],
                          edges_count=len(edges), provenance=provenance)
            screen = screen_graph(record, sat_conflicts=2000)
            record["screen_status"] = screen["status"]
            record["dimension19_eligible"] = screen.get("dimension19_eligible", False)
            if screen["status"] == "excluded":
                state["screened_out"] += 1
                record.update(status="excluded_structurally", screen_reason=screen["reason"],
                              certificate_method=screen["certificate"]["method"],
                              certificate_witness=screen["certificate"]["witness"])
            else:
                remaining = deadline-time.monotonic()
                result = process.solve(30, sorted(edges), min(args.solve_timeout, remaining)) if remaining > .1 else {"error": "batch deadline"}
                if "error" in result:
                    state["solver_errors"] += 1
                    record.update(status="solver_error", error=result["error"])
                else:
                    state["scored"] += 1
                    record.update(status="scored_numerically", theta_numeric=result["theta_numeric"],
                                  solver=result["solver"], solver_status=result["solver_status"])
                    full = dict(record, edges=[list(e) for e in sorted(edges)],
                                parts=[list(range(10*k, 10*k+10)) for k in range(3)], screen=screen,
                                certificate_status="numerical only; exact theta certification required")
                    full.update(result)
                    if best is None or full["theta_numeric"] > best["theta_numeric"]:
                        best = full
                        atomic_json(best_path, best)
                        print("best theta", best["theta_numeric"], "edges", len(edges), "generated", state["generated"], flush=True)
                    if full["theta_numeric"] >= args.near_threshold:
                        atomic_json(near/f"candidate_{state['generated']}.json", full)
            stream.write(json.dumps(record, separators=(",", ":"))+"\n")
            checkpoint()

        try:
            iteration = 0
            while sources and time.monotonic() < deadline-.1:
                source, original, ranked = sources[iteration % len(sources)]
                ordinal = iteration//len(sources)
                if ordinal < min(32, len(ranked)):
                    pairs = ranked[ordinal][1]
                elif rng.random() < .8:
                    pairs = rng.choice(ranked)[1]
                else:
                    pairs = tuple(tuple(sorted(rng.sample(range(10*k, 10*k+10), 2))) for k in range(3))
                base = intersection(original, pairs)
                provenance = dict(source=source, pairs=pairs, iteration=iteration, retained_edges=len(base), method="intersection")
                evaluate(base, provenance)
                # Relabel the original source as well, to guide orbit additions.
                relabel = {}
                for k, pair in enumerate(pairs):
                    order = list(pair)+[v for v in range(10*k, 10*k+10) if v not in pair]
                    relabel.update({v: 10*k+i for i, v in enumerate(order)})
                preferred = {tuple(sorted((relabel[u], relabel[v]))) for u, v in original}
                for extension in range(args.extensions):
                    if time.monotonic() >= deadline-.1:
                        break
                    enlarged = maximal_extension(base, orbits, rng, preferred if extension % 2 == 0 else None)
                    if enlarged is not None:
                        evaluate(enlarged, dict(provenance, method="maximal_orbit_extension", extension=extension,
                                                order="source_guided" if extension % 2 == 0 else "random"))
                iteration += 1
            state["termination"] = "batch time budget" if sources else "no source ranked within budget"
        finally:
            process.stop()
            checkpoint()
    print(json.dumps(state), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seconds", type=float, default=180)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--ranked-pairs", type=int, default=256)
    parser.add_argument("--extensions", type=int, default=4)
    parser.add_argument("--solve-timeout", type=float, default=15)
    parser.add_argument("--near-threshold", type=float, default=2.999)
    args = parser.parse_args()
    if min(args.seconds, args.ranked_pairs, args.solve_timeout) <= 0 or args.extensions < 0:
        parser.error("positive budgets/ranked-pair count and nonnegative extensions required")
    run(args)


if __name__ == "__main__":
    main()

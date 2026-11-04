"""Resumable asymmetric (10,10,10) construction search, with arbitrary SDP weights.

This is a bounded heuristic, never an exhaustive nonexistence test. It searches
triangle-free tripartite supports using edge kicks, vertex rewiring, bipartite
2-switches, population crossover, and fresh starts. No symmetry, uniform Born
probability, regular incidence, or prescribed overlap alphabet is imposed.

  OPENBLAS_NUM_THREADS=1 python search/overnight30_asymmetric.py \
      --seconds 300 --seed 0 --output-dir search/overnight30/asymmetric

The same output directory and seed resume the RNG, current graph and population.
Use --resume-from CHECKPOINT to continue in another output directory; a different
--seed imports the previous population while starting new counters and RNG state.
The shared inbox overnight30_shared_seeds.json is checked on every invocation,
including resumed chunks; newly seen labelled graphs join the pending seed queue.
Every successfully scored labelled support is written compactly as graph6 in a
JSONL archive. Only new bests and numerical near-hits get full graph records.
The bounded elite population removes exact colour-preserving isomorphs within
WL hash buckets; the full archive still counts labelled supports, not classes.
Short checkpointed excursions from near-leading elites use partial clustered
rewires and higher temperature to explore temporary score losses. They preserve
historical bests and the elite pool and return to ordinary search automatically.
CLARABEL has an individual solve time limit; SCS is a bounded fallback. H33 and
H36 are scored as positive controls on each invocation. Numerical near-hits are
saved for subsequent exact certification and are NEVER called realizations here.

Actual positive-probability 30-event supports have at most ten neighbours per
vertex: in each size-10 context every probability is <1/9 and every k-subset has
weight >(k-1)/9; the two mutually orthogonal neighbourhoods have weight <1.
Thus their sizes sum to at most ten. Searching these actual supports loses no
30-event realization, even though supergraphs with extra allowed zeros may have
higher degree. The required one neighbour in each other context and connectedness
are likewise necessary for a minimal positive-support realization. Cross-context
bipartite connectedness is NOT imposed, to allow dimensions above 19 as well.
"""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import random
import time
import warnings

import cvxpy as cp
import networkx as nx
import numpy as np

from overnight30_screen import fast_screen, screen_graph


N = 30
PARTS = [list(range(10*k, 10*k+10)) for k in range(3)]
PART_MASKS = [sum(1 << v for v in part) for part in PARTS]
PAIRS = [(u, v) for u in range(N) for v in range(u+1, N) if u//10 != v//10]
ROOT = Path(__file__).resolve().parent


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def atomic_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    os.replace(temporary, path)


def tuple_tree(value):
    return tuple(tuple_tree(v) for v in value) if isinstance(value, list) else value


def from_edges(edges, n=N):
    adj = [0]*n
    for u, v in edges:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


def edges_of(adj):
    return [(u, v) for u in range(len(adj)) for v in range(u+1, len(adj)) if adj[u] >> v & 1]


def graph_of(adj):
    graph = nx.Graph()
    graph.add_nodes_from(range(len(adj)))
    graph.add_edges_from(edges_of(adj))
    return graph


def encode(adj):
    return nx.to_graph6_bytes(graph_of(adj), header=False).decode().strip()


def decode(code):
    graph = nx.from_graph6_bytes(code.encode())
    if len(graph) != N:
        raise ValueError("checkpoint support does not have thirty vertices")
    return from_edges(graph.edges())


def validate(adj, require_minimum=True):
    if len(adj) != N:
        return False
    for u, neighbours in enumerate(adj):
        if neighbours & PART_MASKS[u//10] or neighbours.bit_count() > 10:
            return False
        if neighbours >> N:
            return False
        if require_minimum and any(not neighbours & PART_MASKS[k] for k in range(3) if k != u//10):
            return False
        for v in range(u+1, N):
            if bool(neighbours >> v & 1) != bool(adj[v] >> u & 1):
                return False
            if neighbours >> v & 1 and neighbours & adj[v]:
                return False
    if require_minimum:
        reached, frontier = 1, 1
        while frontier:
            new = 0
            for u in range(N):
                if frontier >> u & 1:
                    new |= adj[u]
            frontier = new & ~reached
            reached |= frontier
        if reached.bit_count() != N:
            return False
    return True


def addable(adj, u, v):
    return (u//10 != v//10 and not (adj[u] >> v & 1)
            and adj[u].bit_count() < 10 and adj[v].bit_count() < 10
            and not (adj[u] & adj[v]))


def add(adj, u, v):
    adj[u] |= 1 << v
    adj[v] |= 1 << u


def remove(adj, u, v):
    adj[u] &= ~(1 << v)
    adj[v] &= ~(1 << u)


def fill(adj, rng, preferred=()):
    # Prefer covering an absent context-neighbourhood before adding density.
    preferred = list(preferred)
    rng.shuffle(preferred)
    missing = [(u, v) for u, v in PAIRS
               if not adj[u] & PART_MASKS[v//10] or not adj[v] & PART_MASKS[u//10]]
    rng.shuffle(missing)
    rest = PAIRS[:]
    rng.shuffle(rest)
    for u, v in missing + preferred + rest:
        if addable(adj, u, v):
            add(adj, u, v)
    return adj


def kick_edges(current, rng):
    trial = current[:]
    edges = edges_of(trial)
    for u, v in rng.sample(edges, min(len(edges), rng.choice((1, 2, 3, 4, 6, 9)))):
        remove(trial, u, v)
    return fill(trial, rng)


def cluster_rewire(current, rng):
    """Partially rewire two or three adjacent vertices, retaining coverage."""
    for _ in range(8):
        root = rng.randrange(N)
        neighbours = [v for v in range(N) if current[root] >> v & 1]
        cluster = {root, *rng.sample(neighbours, min(len(neighbours), rng.choice((1, 2))))}
        incident = [(u, v) for u, v in edges_of(current) if u in cluster or v in cluster]
        rng.shuffle(incident)
        trial = current[:]
        budget = min(14, max(4, round(len(incident)*rng.uniform(0.35, 0.65))))
        removed = 0
        for u, v in incident:
            if ((trial[u] & PART_MASKS[v//10]).bit_count() > 1
                    and (trial[v] & PART_MASKS[u//10]).bit_count() > 1):
                remove(trial, u, v)
                removed += 1
                if removed >= budget:
                    break
        trial = fill(trial, rng)
        if trial != current and validate(trial):
            return trial, "cluster_partial_rewire"
    return kick_edges(current, rng), "cluster_rewire_fallback_edge_kick"


def proposal(current, pool, rng, iteration, diagnostics=None):
    diagnostics = diagnostics if diagnostics is not None else {}

    def count(name):
        diagnostics[name] = diagnostics.get(name, 0)+1

    r = rng.random()
    if r < 0.42:
        return kick_edges(current, rng), "edge_kick"
    if r < 0.66:
        trial = current[:]
        vertices = rng.sample(range(N), 2 if rng.random() < 0.25 else 1)
        for u in vertices:
            for v in range(N):
                if trial[u] >> v & 1:
                    remove(trial, u, v)
        return fill(trial, rng), "vertex_rewire"
    if r < 0.81:
        trial = current[:]
        a, b = rng.sample(range(3), 2)
        block = [(u, v) for u in PARTS[a] for v in PARTS[b] if trial[u] >> v & 1]
        if len(block) >= 2:
            for _ in range(20):
                count("two_switch_attempts")
                (u, v), (x, y) = rng.sample(block, 2)
                if u == x or v == y or trial[u] >> y & 1 or trial[x] >> v & 1:
                    continue
                candidate = trial[:]
                remove(candidate, u, v)
                remove(candidate, x, y)
                cleared, failed = 0, False
                # An ordinary A--B switch leaves every third-context triangle
                # blocker unchanged. Clear those blockers explicitly first.
                for left, right in ((u, y), (x, v)):
                    common = candidate[left] & candidate[right]
                    for third in range(N):
                        if not common >> third & 1:
                            continue
                        choices = [(left, third), (right, third)]
                        rng.shuffle(choices)
                        choices = [(s, t) for s, t in choices
                                   if (candidate[s] & PART_MASKS[t//10]).bit_count() > 1
                                   and (candidate[t] & PART_MASKS[s//10]).bit_count() > 1]
                        if not choices:
                            failed = True
                            break
                        remove(candidate, *choices[0])
                        cleared += 1
                    if failed or not addable(candidate, left, right):
                        failed = True
                        break
                    add(candidate, left, right)
                if failed:
                    count("two_switch_conflict_failures")
                    continue
                candidate = fill(candidate, rng)
                if not validate(candidate):
                    count("two_switch_invalid_trials")
                    continue
                if candidate == current:
                    count("two_switch_unchanged_trials")
                    continue
                count("two_switch_conflict_edges_removed_" + str(cleared))
                return candidate, "conflict_clearing_two_switch" if cleared else "two_switch"
        count("two_switch_fallbacks")
        return kick_edges(current, rng), "two_switch_fallback_edge_kick"
    if r < 0.96 and pool:
        second = decode(rng.choice(pool)["g6"])
        trial = [a & b for a, b in zip(current, second)]
        union = [(u, v) for u, v in PAIRS if (current[u] | second[u]) >> v & 1]
        return fill(trial, rng, union), "population_crossover"
    # No empty-graph restarts were accepted after per-move counters were added;
    # spend that mutation share near the current geometry instead.
    return cluster_rewire(current, rng)


class ThetaScorer:
    def __init__(self, n):
        self.n = n
        self.B = cp.Variable((n, n), symmetric=True)
        self.mask = cp.Parameter((n, n))
        self.problem = cp.Problem(cp.Maximize(cp.sum(self.B)),
                                  [self.B >> 0, cp.trace(self.B) == 1,
                                   cp.multiply(self.mask, self.B) == 0])

    def solve(self, adj, limit=5.0):
        mask = np.ones((self.n, self.n))
        np.fill_diagonal(mask, 0)
        for u, v in edges_of(adj):
            mask[u, v] = mask[v, u] = 0
        self.mask.value = mask
        errors = []
        begin = time.monotonic()
        for solver in ("CLARABEL", "SCS"):
            remaining = limit-(time.monotonic()-begin)
            if remaining <= 0:
                break
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    if solver == "CLARABEL":
                        self.problem.solve(solver=solver, warm_start=True, max_iter=150,
                                           time_limit=max(0.05, remaining))
                    else:
                        self.problem.solve(solver=solver, warm_start=True, eps=2e-6,
                                           max_iters=12000, time_limit_secs=max(0.05, remaining))
                value = self.problem.value
                if value is None or not math.isfinite(value) or self.B.value is None:
                    errors.append(solver + ":" + str(self.problem.status))
                    continue
                B = np.asarray(self.B.value)
                diagnostics = {
                    "minimum_eigenvalue": float(np.linalg.eigvalsh(B)[0]),
                    "trace_error": float(abs(np.trace(B)-1)),
                    "zero_error": float(np.max(np.abs(mask*B))),
                }
                if (self.problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE)
                        or not np.isfinite(B).all()
                        or diagnostics["minimum_eigenvalue"] < -2e-5
                        or diagnostics["trace_error"] > 2e-5
                        or diagnostics["zero_error"] > 2e-5
                        or value > 3.0001):
                    errors.append(solver + ":unreliable_residuals_or_status")
                    continue
                return {"theta": float(value), "solver": solver,
                        "status": str(self.problem.status), "diagnostics": diagnostics,
                        "B": B, "seconds": time.monotonic()-begin}
            except (cp.error.SolverError, ValueError, ArithmeticError) as exc:
                errors.append(solver + ":" + type(exc).__name__)
        return {"error": errors, "seconds": time.monotonic()-begin}


def positive_controls():
    from verify_all import tricayley_z2z6
    H33 = nx.read_graph6(ROOT / "h33.g6")
    S = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)),
         ((0, 1), (1, 0), (1, 1), (1, 3)))
    vertices, edges = tricayley_z2z6(S)
    index = {v: i for i, v in enumerate(vertices)}
    configurations = {"H33": from_edges(H33.edges(), 33),
                      "H36": from_edges([(index[u], index[v]) for u, v in edges], 36)}
    values = {}
    for name, adj in configurations.items():
        result = ThetaScorer(len(adj)).solve(adj, limit=8)
        if "error" in result or abs(result["theta"]-3) > 5e-5:
            raise RuntimeError(f"positive control failed: {name}: {result}")
        values[name] = {"theta": result["theta"], "diagnostics": result["diagnostics"]}
    return values


def incidence_edges(B, C):
    B, C = list(map(set, B)), list(map(set, C))
    edges = [(p, 10+j) for j, line in enumerate(B) for p in line]
    edges += [(p, 20+j) for j, line in enumerate(C) for p in line]
    edges += [(10+j, 20+k) for j in range(10) for k in range(10) if not B[j] & C[k]]
    return edges


def seed_graphs(seed):
    rng = random.Random(seed ^ 0x10A10B10C)
    result = []
    source = ROOT / "gap_search_results.json"
    if source.exists():
        for row in json.loads(source.read_text())["runs"]:
            parts = row["parts"]
            for _ in range(1 if all(len(p) == 10 for p in parts) else 3):
                selected = [rng.sample(p, 10) for p in parts]
                mapping = {v: 10*k+i for k, part in enumerate(selected) for i, v in enumerate(part)}
                result.append(from_edges([(mapping[u], mapping[v]) for u, v in row["edges"]
                                          if u in mapping and v in mapping]))
    for filename in ("construction30_mixed_results.json", "construction30_nonlinear_results.json"):
        path = ROOT / filename
        if path.exists():
            for row in json.loads(path.read_text())["runs"]:
                result.append(from_edges(incidence_edges(row["B_lines"], row["C_lines"])))
    source = ROOT / "event30_uniform_results.json"
    if source.exists():
        row = json.loads(source.read_text()).get("best")
        if row:
            result.append(from_edges(row["edges"]))
    for filename in ("overnight30_seed.json", "overnight30_screened_seed.json"):
        source = ROOT / filename
        if source.exists():
            result.append(from_edges(json.loads(source.read_text())["edges"]))
    # Fill each seed to a maximal triangle-free support under the necessary cap.
    deduplicated = {}
    for adj in result:
        if validate(adj, require_minimum=False):
            adj = fill(adj, rng)
            if validate(adj):
                deduplicated.setdefault(encode(adj), adj)
    rng.shuffle(result := list(deduplicated))
    return result


def shared_seed_graph(record):
    """Normalize explicit context labels before interpreting a shared graph."""
    if "edges" in record and "parts" in record:
        parts = record["parts"]
        vertices = [v for part in parts for v in part]
        if list(map(len, parts)) != [10, 10, 10] or len(set(vertices)) != N:
            raise ValueError("shared seed must partition thirty vertices into three tens")
        mapping = {v: 10*k+i for k, part in enumerate(parts) for i, v in enumerate(part)}
        adj = from_edges([(mapping[u], mapping[v]) for u, v in record["edges"]])
    else:
        adj = decode(record.get("g6") or record["graph6"])
    if not validate(adj):
        raise ValueError("shared seed fails exact support validation")
    return adj


def inject_shared_seeds(state, pool, current, path=None):
    """Append once per labelled graph without replacing resumed search state."""
    path = path or ROOT / "overnight30_shared_seeds.json"
    seen = set(state.get("shared_seed_injected_g6", []))
    queued = set(state["pending_seeds"])
    queued.update(entry["g6"] for entry in pool)
    if current is not None:
        queued.add(encode(current))
    errors, added = [], 0
    if path.exists():
        try:
            payload = json.loads(path.read_text())
            rows = payload["seeds"]
            if payload.get("format") != 1 or not isinstance(rows, list):
                raise ValueError("unsupported shared seed inbox format")
            # Pending seeds are popped from the end; score the strongest first.
            rows = sorted(rows, key=lambda row: row.get("theta_numeric", row.get("theta", 0)))
            for row in rows:
                try:
                    code = encode(shared_seed_graph(row))
                    if code in seen:
                        continue
                    seen.add(code)
                    if code not in queued:
                        state["pending_seeds"].append(code)
                        queued.add(code)
                        added += 1
                except (ValueError, KeyError, TypeError, nx.NetworkXError) as exc:
                    errors.append(f"{row.get('source', 'shared seed')}: {type(exc).__name__}: {exc}")
        except (ValueError, KeyError, TypeError, OSError) as exc:
            errors.append(f"shared inbox: {type(exc).__name__}: {exc}")
    state.update({"shared_seed_injected_g6": sorted(seen),
                  "shared_seed_last_added": added, "shared_seed_last_checked_utc": utc_now(),
                  "shared_seed_errors": errors[:10]})
    return added


def diversity_key(adj):
    graph = graph_of(adj)
    nx.set_node_attributes(graph, {v: v//10 for v in graph}, "part")
    return nx.weisfeiler_lehman_graph_hash(graph, node_attr="part", iterations=3)


def coloured_isomorphic(left, right):
    first, second = graph_of(left), graph_of(right)
    for graph in (first, second):
        nx.set_node_attributes(graph, {v: v//10 for v in graph}, "part")
    return nx.is_isomorphic(first, second,
                            node_match=nx.algorithms.isomorphism.categorical_node_match("part", -1))


def screen_summary(result):
    return {key: result[key] for key in ("status", "reason", "dimension19_eligible", "coloring_search")
            if key in result}


def screened(adj, value=None, threshold=2.85, conflicts=2000):
    record = {"parts": PARTS, "edges": edges_of(adj)}
    if value is not None and value >= threshold:
        return screen_graph(record, sat_conflicts=conflicts)
    return fast_screen(record)


def search_score(entry):
    # A small explicitly excluded bridge population may escape local basins.
    return entry["theta"] - (0.06 if entry.get("screen", {}).get("status") == "excluded" else 0)


def run(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_dir / f"asymmetric_seed{args.seed}"
    lock = prefix.with_suffix(".lock").open("w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print(json.dumps({"status": "already_running", "seed": args.seed}))
        return
    begin = time.monotonic()
    deadline = begin + args.seconds
    checkpoint = prefix.with_suffix(".json")
    archive = prefix.with_suffix(".jsonl")
    controls = positive_controls()
    rng = random.Random(args.seed)
    source = checkpoint if checkpoint.exists() else args.resume_from
    previous = json.loads(source.read_text()) if source else None
    if previous is not None and previous.get("format") != 1:
        raise ValueError("incompatible checkpoint format")
    if previous is not None and previous["seed"] == args.seed:
        state = previous
        if state["seed"] != args.seed:
            raise ValueError("incompatible checkpoint")
        rng.setstate(tuple_tree(state["rng_state"]))
        pool = state["pool"]
        current = decode(state["current"]) if state.get("current") else None
        current_value = state.get("current_theta")
        state["resumed_from"] = str(source.resolve())
    else:
        state = {"format": 1, "seed": args.seed,
                 "scope": "bounded asymmetric construction heuristic, not an exclusion proof",
                 "started_utc": utc_now(), "pending_seeds": seed_graphs(args.seed),
                 "iterations": 0, "scored": 0, "invalid": 0, "duplicate": 0,
                 "solver_failures": 0, "accepted": 0, "near_hits": 0,
                 "moves": {}, "elapsed_seconds": 0.0, "best_theta": None,
                 "packages": {p: version(p) for p in ("numpy", "networkx", "cvxpy", "clarabel")}}
        pool = previous["pool"] if previous else []
        current = decode(previous["current"]) if previous and previous.get("current") else None
        current_value = previous.get("current_theta") if previous else None
        if previous:
            state["warm_started_from"] = str(source.resolve())
            state["shared_seed_injected_g6"] = previous.get("shared_seed_injected_g6", [])
    for entry in pool:
        if not validate(decode(entry["g6"])) or not math.isfinite(entry["theta"]):
            raise ValueError("invalid population entry in checkpoint")
    if current is not None and (not validate(current) or current_value is None):
        raise ValueError("invalid current graph in checkpoint")
    injected = inject_shared_seeds(state, pool, current)
    if injected or state["shared_seed_errors"]:
        print(json.dumps({"seed": args.seed, "shared_seeds_added": injected,
                          "shared_seed_errors": state["shared_seed_errors"]}), flush=True)
    state["termination"] = "running"
    state.setdefault("move_accepted", {})
    state.setdefault("move_invalid", {})
    state.setdefault("proposal_diagnostics", {})
    state.setdefault("excursions", {"remaining_evaluations": 0, "next_iteration": state["iterations"],
                                   "started": 0, "completed": 0, "evaluated": 0,
                                   "newly_scored": 0, "accepted": 0, "downhill_accepted": 0,
                                   "maximum_anchor_drop": 0.0})
    old_screen_version = state.get("screen_version", 0)
    state.update({"screen_version": 1, "fast_rejected": state.get("fast_rejected", 0),
                  "sat_excluded": state.get("sat_excluded", 0),
                  "best_excluded_theta": state.get("best_excluded_theta"),
                  "best_theta_meaning": "highest theta among graphs not excluded by the configured exact screens"})
    state["last_controls"] = controls
    state["last_started_utc"] = utc_now()
    cache = {entry["g6"]: entry["theta"] for entry in pool}
    screen_cache = {}
    if current is not None:
        cache[encode(current)] = current_value
    if archive.exists():
        for line in archive.open():
            try:
                record = json.loads(line)
                if "theta" in record:
                    cache[record["g6"]] = record["theta"]
                if "screen" in record:
                    screen_cache[record["g6"]] = record["screen"]
            except (json.JSONDecodeError, KeyError):
                # An interrupted final append is harmless: checkpoints are atomic.
                continue
    scorer = ThetaScorer(N)
    last_saved, last_improvement = begin, state["iterations"]
    total_before = state["elapsed_seconds"]

    def save():
        state.update({"pool": pool, "current": encode(current) if current else None,
                      "current_theta": current_value, "rng_state": rng.getstate(),
                      "elapsed_seconds": total_before+time.monotonic()-begin,
                      "updated_utc": utc_now(), "labelled_graphs_cached": len(cache)})
        atomic_json(checkpoint, state)

    def retain(code, value, adj, screening):
        nonlocal pool
        if any(entry["g6"] == code for entry in pool):
            return
        # WL is ONLY a diversity bucket, never an isomorphism or rejection proof.
        key = diversity_key(adj)
        same = [entry for entry in pool if entry["wl"] == key]
        for entry in same:
            if coloured_isomorphic(adj, decode(entry["g6"])):
                if value > entry["theta"]:
                    entry.update({"g6": code, "theta": value, "screen": screen_summary(screening)})
                return
        if len(same) >= 3:
            weakest = min(same, key=lambda entry: entry["theta"])
            if weakest["theta"] >= value:
                return
            pool.remove(weakest)
        pool.append({"g6": code, "theta": value, "wl": key, "screen": screen_summary(screening)})
        pool.sort(key=search_score, reverse=True)
        excluded = [entry for entry in pool if entry["screen"]["status"] == "excluded"]
        for entry in excluded[max(2, args.pool_size//8):]:
            pool.remove(entry)
        if len(pool) > args.pool_size:
            # Keep the best half and a stratified tail, rather than only clones
            # of the leading local optimum.
            keep = pool[:args.pool_size//2]
            tail = pool[args.pool_size//2:]
            keep += rng.sample(tail, args.pool_size-len(keep))
            pool = sorted(keep, key=search_score, reverse=True)

    # Old checkpoints may predate exact isomorphism deduplication of the pool.
    imported_pool, pool = pool, []
    for entry in imported_pool:
        adj = decode(entry["g6"])
        screening = screened(adj, entry["theta"], args.sat_threshold, args.sat_conflicts)
        screen_cache[entry["g6"]] = screen_summary(screening)
        retain(entry["g6"], entry["theta"], adj, screening)
    if current is not None:
        screen_cache[encode(current)] = screen_summary(screened(current, current_value, args.sat_threshold, args.sat_conflicts))
    state["pool_deduplication"] = "exact isomorphism preserving the three labelled contexts"
    if not old_screen_version:
        state["best_theta"] = None
        old_best_path = prefix.with_name(prefix.name+"_best.json")
        if old_best_path.exists():
            os.replace(old_best_path, prefix.with_name(prefix.name+"_legacy_best_unscreened.json"))
    for excluded in (False, True):
        eligible = [entry for entry in pool if (entry["screen"]["status"] == "excluded") == excluded]
        if eligible:
            leader = max(eligible, key=lambda entry: entry["theta"])
            field = "best_excluded_theta" if excluded else "best_theta"
            state[field] = max(state[field] or -math.inf, leader["theta"])
            best_path = prefix.with_name(prefix.name+("_best_excluded.json" if excluded else "_best.json"))
            if not best_path.exists():
                screening = screened(decode(leader["g6"]), leader["theta"], args.sat_threshold, args.sat_conflicts)
                atomic_json(best_path, {"scope": "excluded mutation seed" if excluded else "numerical candidate, not a certified realization",
                                        "theta_numeric": leader["theta"], "g6": leader["g6"],
                                        "parts": PARTS, "edges": edges_of(decode(leader["g6"])),
                                        "screen": screen_summary(screening), "certificate": screening.get("certificate"),
                                        "seed": args.seed, "source_checkpoint": str(source),
                                        "updated_utc": utc_now()})
    save()
    with archive.open("a", buffering=1) as stream:
        while time.monotonic() < deadline-0.3:
            state["iterations"] += 1
            exploring = False
            excursion = state["excursions"]
            if excursion["remaining_evaluations"] and state["iterations"] >= excursion["stop_iteration"]:
                excursion["remaining_evaluations"] = 0
                excursion["completed"] += 1
            if state["pending_seeds"]:
                trial = decode(state["pending_seeds"].pop())
                move = "saved_seed"
            elif current is None:
                trial, move = fill([0]*N, rng), "fresh_restart"
            else:
                if pool and not excursion["remaining_evaluations"] and state["iterations"] >= excursion["next_iteration"]:
                    lead_score = max(map(search_score, pool))
                    near_leaders = [entry for entry in pool if search_score(entry) >= lead_score-0.035]
                    selected = rng.choice(near_leaders[:16])
                    current, current_value = decode(selected["g6"]), selected["theta"]
                    excursion.update({"remaining_evaluations": 50, "next_iteration": state["iterations"]+650,
                                      "stop_iteration": state["iterations"]+240,
                                      "anchor_g6": selected["g6"], "anchor_theta": selected["theta"],
                                      "anchor_search_score": search_score(selected),
                                      "started_iteration": state["iterations"]})
                    excursion["started"] += 1
                exploring = excursion["remaining_evaluations"] > 0
                if not exploring and pool and (state["iterations"] % 120 == 0 or state["iterations"]-last_improvement > 250):
                    selected = rng.choice(pool[:max(4, len(pool)//2)])
                    current, current_value = decode(selected["g6"]), selected["theta"]
                    last_improvement = state["iterations"]
                if exploring and rng.random() < 0.6:
                    trial, move = cluster_rewire(current, rng)
                else:
                    trial, move = proposal(current, pool, rng, state["iterations"], state["proposal_diagnostics"])
            state["moves"][move] = state["moves"].get(move, 0)+1
            state["move_accepted"].setdefault(move, 0)
            state["move_invalid"].setdefault(move, 0)
            if not validate(trial):
                state["invalid"] += 1
                state["move_invalid"][move] += 1
                continue
            code = encode(trial)
            screening = screen_cache.get(code) or screened(trial)
            if screening["status"] == "excluded" and rng.random() >= 0.1:
                state["fast_rejected"] += 1
                if code not in screen_cache:
                    screen_cache[code] = screen_summary(screening)
                    stream.write(json.dumps({"g6": code, "seed": args.seed, "iteration": state["iterations"],
                                             "screen": screen_summary(screening)}, separators=(",", ":"))+"\n")
                continue
            if code in cache:
                state["duplicate"] += 1
                value = cache[code]
                result = None
            else:
                remaining = deadline-time.monotonic()
                if remaining < 0.3:
                    break
                result = scorer.solve(trial, limit=min(args.solve_seconds, remaining))
                if "error" in result:
                    state["solver_failures"] += 1
                    stream.write(json.dumps({"g6": code, "error": result["error"],
                                             "iteration": state["iterations"]}, separators=(",", ":"))+"\n")
                    continue
                value = result["theta"]
                cache[code] = value
                state["scored"] += 1
                if screening["status"] != "excluded" and value >= args.sat_threshold:
                    screening = screened(trial, value, args.sat_threshold, args.sat_conflicts)
                    if screening["status"] == "excluded":
                        state["sat_excluded"] += 1
                screen_cache[code] = screen_summary(screening)
                stream.write(json.dumps({"g6": code, "theta": value, "seed": args.seed,
                                         "iteration": state["iterations"], "move": move,
                                         "status": result["status"], "screen": screen_summary(screening),
                                         **({"excursion": True} if exploring else {})}, separators=(",", ":"))+"\n")
            # Cached legacy scores must also receive the high-score SAT screen.
            if code not in screen_cache:
                screening = screened(trial, value, args.sat_threshold, args.sat_conflicts)
                screen_cache[code] = screen_summary(screening)
            else:
                screening = screen_cache[code]
            excluded = screening["status"] == "excluded"
            retain(code, value, trial, screening)
            field = "best_excluded_theta" if excluded else "best_theta"
            if state[field] is None or value > state[field] + 1e-8:
                state[field] = value
                last_improvement = state["iterations"]
                screening = screened(trial, value, args.sat_threshold, args.sat_conflicts)
                best = {"scope": "excluded mutation seed" if excluded else "numerical candidate, not a certified realization",
                        "theta_numeric": value, "g6": code, "parts": PARTS,
                        "edges": edges_of(trial), "seed": args.seed,
                        "iteration": state["iterations"], "updated_utc": utc_now(),
                        "screen": screen_summary(screening), "certificate": screening.get("certificate")}
                if result:
                    best["solver"] = result["solver"]
                    best["diagnostics"] = result["diagnostics"]
                    best["diagonal_times_three"] = (3*np.diag(result["B"])).tolist()
                atomic_json(prefix.with_name(prefix.name+("_best_excluded.json" if excluded else "_best.json")), best)
                print(json.dumps({"seed": args.seed, field: value,
                                  "scored": state["scored"], "move": move}), flush=True)
            if value >= 3-args.hit_tolerance and result:
                state["near_hits"] += 1
                hit = {"scope": "numerical near-hit requiring exact verification",
                       "theta_numeric": value, "g6": code, "parts": PARTS,
                       "edges": edges_of(trial), "Gram_primal": result["B"].tolist(),
                       "diagnostics": result["diagnostics"], "solver": result["solver"],
                       "status": result["status"], "screen": screen_summary(screening), "created_utc": utc_now()}
                digest = hashlib.sha256(code.encode()).hexdigest()[:16]
                atomic_json(args.output_dir / f"candidate_{digest}.json", hit)
                save()
                print(json.dumps({"numerical_candidate": digest, "theta": value}), flush=True)
                if args.stop_on_hit:
                    state["termination"] = "numerical_candidate_needs_exact_check"
                    break
            # Cyclic cooling allows indefinitely resumed search without requiring
            # a predeclared total iteration count.
            phase = state["iterations"] % 200
            temperature = 0.001 + 0.016*(1-phase/200)
            if exploring:
                temperature = 0.015 + 0.035*excursion["remaining_evaluations"]/50
            adjusted = value - (0.06 if excluded else 0)
            adjusted_current = current_value - (0.06 if current is not None and screen_cache[encode(current)]["status"] == "excluded" else 0) if current_value is not None else -math.inf
            within_excursion_band = not exploring or adjusted >= excursion["anchor_search_score"]-0.10
            if within_excursion_band and (current is None or adjusted >= adjusted_current or rng.random() < math.exp(min(0, (adjusted-adjusted_current)/temperature))):
                if exploring:
                    excursion["accepted"] += 1
                    excursion["downhill_accepted"] += int(adjusted < adjusted_current-1e-8)
                    excursion["maximum_anchor_drop"] = max(excursion["maximum_anchor_drop"], excursion["anchor_search_score"]-adjusted)
                current, current_value = trial, value
                state["accepted"] += 1
                state["move_accepted"][move] += 1
            if exploring:
                excursion["evaluated"] += 1
                excursion["newly_scored"] += int(result is not None)
                excursion["remaining_evaluations"] -= 1
                if not excursion["remaining_evaluations"]:
                    excursion["completed"] += 1
                    last_improvement = state["iterations"]
            if time.monotonic()-last_saved >= 10:
                save()
                last_saved = time.monotonic()
        else:
            state["termination"] = "time_budget"
    if state["termination"] == "running":
        state["termination"] = "time_budget"
    save()
    print(json.dumps({k: state[k] for k in ("seed", "iterations", "scored", "invalid", "duplicate",
                                           "solver_failures", "best_theta", "best_excluded_theta", "fast_rejected", "sat_excluded", "near_hits", "termination")}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=float, default=300)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--resume-from", type=Path,
                        help="Import a checkpoint; same seed resumes RNG and counters, different seed imports population")
    parser.add_argument("--pool-size", type=int, default=64)
    parser.add_argument("--solve-seconds", type=float, default=5)
    parser.add_argument("--hit-tolerance", type=float, default=1e-6)
    parser.add_argument("--sat-threshold", type=float, default=2.85)
    parser.add_argument("--sat-conflicts", type=int, default=2000)
    parser.add_argument("--stop-on-hit", action="store_true")
    args = parser.parse_args()
    if args.seconds <= 0 or args.pool_size < 8 or args.solve_seconds <= 0:
        parser.error("positive time budgets and pool-size >= 8 are required")
    run(args)

"""Resumable SAT/theta search with weak C2, C3, and C5 support symmetries.

The three contexts always have ten vertices. Symmetry has multiple orbits and/or
fixed points, so this is broader than regular tri-Cayley constructions. Event
probabilities and overlap weights are free. All SDP results are numerical leads,
never nonexistence certificates. Each candidate is saved compactly as graph6.

Examples (run with the project's scientific Python environment):
  python search/overnight30_symmetry.py --controls
  python search/overnight30_symmetry.py --family c2_mixed --seed 1 \
      --seconds 300 --output-dir search/overnight30/pools

Repeat the same command to resume. A per-batch file lock prevents simultaneous
writers. SAT cuts and examined support assignments are replayed from JSONL.
Only best/near-3 records contain full edges and primal matrices. Numerical solves
run in a subprocess with a hard timeout; every batch respects its own deadline.
"""
import argparse
from collections import Counter, defaultdict
import fcntl
import gzip
import hashlib
from importlib.metadata import version
import itertools
import json
import multiprocessing as mp
import os
from pathlib import Path
import random
import time

# Every worker gets one BLAS thread even if the execution host has many cores.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

FAMILIES = {
    "c2_pairs": ((2, 2, 2, 2, 2),)*3,
    "c2_fixed": ((2, 2, 2, 1, 1, 1, 1),)*3,
    "c2_mixed": ((2, 2, 2, 2, 2), (2, 2, 2, 1, 1, 1, 1), (2, 1, 1, 1, 1, 1, 1, 1, 1)),
    "c3_fixed": ((3, 3, 3, 1),)*3,
    "c3_mixed": ((3, 3, 3, 1), (3, 3, 1, 1, 1, 1), (3, 1, 1, 1, 1, 1, 1, 1)),
    "c5_double": ((5, 5),)*3,
    "c5_fixed": ((5, 1, 1, 1, 1, 1),)*3,
    "c5_mixed": ((5, 5), (5, 1, 1, 1, 1, 1), (5, 1, 1, 1, 1, 1)),
    "c2_onepair": ((2, 1, 1, 1, 1, 1, 1, 1, 1),)*3,
    "c2_twopairs": ((2, 2, 1, 1, 1, 1, 1, 1),)*3,
    "c3_onecycle": ((3, 1, 1, 1, 1, 1, 1, 1),)*3,
}

STRUCTURAL_OBSTRUCTIONS = {"bridge_shortage", "unbridged_vertex", "degree_mass", "probability_propagation", "regular_pair", "c5_quotient", "orbit_probability_quotient"}


def generator(family):
    permutation = []
    for part, sizes in enumerate(FAMILIES[family]):
        assert sum(sizes) == 10
        offset = part*10
        for size in sizes:
            permutation.extend(offset+(i+1) % size for i in range(size))
            offset += size
    assert sorted(permutation) == list(range(30))
    return permutation


def orbit_positions(orbits):
    """graph6 stores upper-triangular adjacency bits column by column."""
    positions = []
    for orbit in orbits:
        u, v = orbit[0]
        bit = v*(v-1)//2+u
        positions.append((1+bit//6, 1 << (5-bit % 6)))
    return positions


def orbit_bits_from_graph6(graph6, positions):
    data = graph6.encode("ascii")
    if len(data) != 74 or data[0] != 93:
        raise ValueError("expected a headerless 30-vertex graph6 record")
    return "".join("1" if (data[index]-63) & mask else "0" for index, mask in positions)


def build(family):
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    from pysat.solvers import Cadical153
    permutation = generator(family)
    pairs = [(u, v) for u in range(30) for v in range(u+1, 30) if u//10 != v//10]
    edge_variables, orbits = {}, []
    for pair in pairs:
        if pair in edge_variables:
            continue
        orbit, current = [], pair
        while current not in orbit:
            orbit.append(current)
            current = tuple(sorted((permutation[current[0]], permutation[current[1]])))
        orbits.append(sorted(orbit))
        for edge in orbit:
            edge_variables[edge] = len(orbits)
    pool = IDPool(start_from=len(orbits)+1)
    clauses = []

    def edge(u, v):
        return edge_variables[tuple(sorted((u, v)))]

    for a, b, c in itertools.product(range(10), range(10, 20), range(20, 30)):
        clauses.append([-edge(a, b), -edge(a, c), -edge(b, c)])
    # An isolated K1,1 pair component has no third-context bridge: a bridge
    # would complete a triangle. Its positive common-state piece is impossible.
    # This strengthens the bridge screen without allocating auxiliary variables.
    for u, v in pairs:
        variable = edge(u, v)
        alternatives = {edge(u, w) for w in range(10*(v//10), 10*(v//10+1)) if w != v}
        alternatives.update(edge(w, v) for w in range(10*(u//10), 10*(u//10+1)) if w != u)
        if variable not in alternatives:  # otherwise the orbit clause is tautological
            clauses.append([-variable]+sorted(alternatives))
    if family == "c2_onepair":
        # Proposition 39: the six moved vertices must induce C6. Each of the
        # three invariant 2x2 blocks is one matching; triangle clauses select
        # the parity yielding C6 instead of two triangles.
        for k, l in itertools.combinations(range(3), 2):
            a, b = edge(10*k, 10*l), edge(10*k, 10*l+1)
            clauses.extend(([a, b], [-a, -b]))
    for u in range(30):
        incident = []
        for part in range(3):
            if part == u//10:
                continue
            lits = [edge(u, v) for v in range(10*part, 10*(part+1))]
            clauses.append(sorted(set(lits)))  # at least one neighbour per other context
            incident.extend(lits)
        # Repeated literals count separate incident edges in the sequential
        # counter. This is a weighted orbit-degree constraint, not orbit count.
        clauses.extend(CardEnc.atmost(incident, bound=10, vpool=pool,
                                     encoding=EncType.seqcounter).clauses)
    # Equal neighbourhoods in either other context preclude a length-three
    # path: any alleged path would complete a triangle at one endpoint.
    # Proposition 39 therefore requires distinct rows in BOTH pair supports.
    for part in range(3):
        for u, v in itertools.combinations(range(10*part, 10*(part+1)), 2):
            for other in range(3):
                if other == part:
                    continue
                pairs = {tuple(sorted((edge(u, w), edge(v, w))))
                         for w in range(10*other, 10*other+10) if edge(u, w) != edge(v, w)}
                different = []
                for a, b in sorted(pairs):
                    d = pool.id()
                    different.append(d)
                    clauses.extend(([-a, -b, -d], [a, b, -d], [a, -b, d], [-a, b, d]))
                clauses.append(different)
    if family == "c5_double":
        from overnight30_quotient import quotient_exclusion_cnf
        quotient_clauses, _ = quotient_exclusion_cnf(edge_variables, pool)
        clauses.extend(quotient_clauses)
    return Cadical153(bootstrap_with=clauses), edge_variables, orbits, permutation


def check_graph(graph, permutation):
    import networkx as nx
    assert set(graph) == set(range(30))
    assert all(u//10 != v//10 for u, v in graph.edges())
    assert not any(nx.triangles(graph).values())
    for u in graph:
        assert graph.degree(u) <= 10
        assert all(any(v//10 == part for v in graph[u])
                   for part in range(3) if part != u//10)
    assert all(graph.has_edge(permutation[u], permutation[v]) for u, v in graph.edges())


def missing_three_path(graph):
    """A missing same-context length-three path is an exact contraction obstruction."""
    adj = [set(graph[v]) for v in range(30)]
    for part in range(3):
        for u, v in itertools.combinations(range(10*part, 10*part+10), 2):
            if not any(adj[a] & adj[v] for a in adj[u]):
                return u, v
    return None


def path_requirement(solver, edge_variables, pair):
    """Existential CNF for one of the 200 length-three paths; no graph ansatz."""
    u, v = pair
    assert 0 <= u < v < 30 and u//10 == v//10
    parts = [k for k in range(3) if k != u//10]
    edge = lambda a, b: edge_variables[tuple(sorted((a, b)))]
    terms = set()
    for k, l in (parts, parts[::-1]):
        for b, c in itertools.product(range(10*k, 10*k+10), range(10*l, 10*l+10)):
            terms.add(tuple(sorted({edge(u, b), edge(b, c), edge(c, v)})))
    start = solver.nof_vars()+1
    clauses = []
    for auxiliary, term in enumerate(sorted(terms), start):
        clauses.extend([-auxiliary, literal] for literal in term)
    clauses.append(list(range(start, start+len(terms))))
    return dict(reason="same_context_pair_requires_length3_path", witness=list(pair),
                aux_start=start, aux_end=start+len(terms)-1, clauses=clauses)


def dominance_clause(certificate, orbits):
    """Exclude every support contained in the exactly theta-bounded graph."""
    if certificate["n"] != 30:
        raise ValueError("dominance certificate must have 30 vertices")
    edges = {tuple(edge) for edge in certificate["edges"]}
    clause = [i for i, orbit in enumerate(orbits, 1) if not all(edge in edges for edge in orbit)]
    if not clause:
        raise ValueError("empty dual dominance clause")
    return clause


def replay_structural_dominance(cut, orbits):
    """Reconstruct and verify a compact structural certificate using only stdlib."""
    from overnight30_screen import verify_certificate
    compact = cut["witness"]
    if compact["method"] not in STRUCTURAL_OBSTRUCTIONS:
        raise ValueError("unsupported global structural dominance method")
    data = compact["graph6"].encode("ascii")
    if len(data) != 74 or data[0] != 93 or any(x < 63 or x > 126 for x in data):
        raise ValueError("expected a headerless 30-vertex graph6 record")
    edges = []
    for v in range(1, 30):
        for u in range(v):
            bit = v*(v-1)//2+u
            if (data[1+bit//6]-63) & (1 << (5-bit % 6)):
                edges.append([u, v])
    graph = {"n": 30, "parts": [list(range(10*k, 10*k+10)) for k in range(3)],
             "edges": sorted(edges)}
    local_method = compact["method"] in {"regular_pair", "c5_quotient", "orbit_probability_quotient"}
    certificate = {"format": "overnight30-worker-obstruction-v1" if local_method else "overnight30-certificate-v1", "graph": graph,
                   "method": compact["method"], "witness": compact["witness"]}
    if local_method:
        from overnight30_quotient import verify_certificate as verify_worker_certificate
        verify_worker_certificate(certificate)
    else:
        verify_certificate(certificate)
    clause = dominance_clause(graph, orbits)
    if cut["clause"] != clause:
        raise ValueError("structural dominance clause does not match its certificate")
    return clause


def worker_structural_screen(record):
    from overnight30_quotient import find_obstruction
    from overnight30_screen import fast_screen
    local = find_obstruction(record)
    if local is not None:
        return {"status": "excluded", "reason": "exact regular-pair or automorphism-quotient obstruction",
                "certificate": local}
    return fast_screen(record)


def replay_dominance(source, cut, orbits):
    from overnight30_fast_dual import verify_certificate
    import networkx as nx
    witness = cut["witness"]
    relative = Path(witness["certificate"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("dominance certificate path must stay inside its chunk")
    raw = gzip.decompress((source/relative).read_bytes())
    if hashlib.sha256(raw).hexdigest() != witness["sha256"]:
        raise ValueError("dominance certificate hash mismatch")
    certificate = json.loads(raw)
    bound = verify_certificate(certificate)
    if str(bound) != witness["upper_bound"]:
        raise ValueError("dominance bound mismatch")
    graph = nx.Graph()
    graph.add_nodes_from(range(30))
    graph.add_edges_from(certificate["edges"])
    if nx.to_graph6_bytes(graph, header=False).decode().strip() != witness["graph6"]:
        raise ValueError("dominance certificate graph mismatch")
    clause = dominance_clause(certificate, orbits)
    if clause != cut["clause"]:
        raise ValueError("dominance clause does not match its certificate")
    return clause


def bad_coloring(graph):
    """Return an independently checked proper coloring with a class <=9.

    Such a witness excludes theta3 on30 vertices by Corollary6/Theorem5, even
    when the candidate's optimal SDP weights are not uniform. Unknown SAT
    outcomes merely leave the candidate available for numerical scoring.
    """
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    from pysat.solvers import Cadical153
    pool = IDPool(start_from=91)
    x = lambda v, color: 3*v+color+1
    clauses = []
    for v in range(30):
        clauses.append([x(v, k) for k in range(3)])
        clauses.extend([-x(v, k), -x(v, l)] for k, l in itertools.combinations(range(3), 2))
    clauses.extend([-x(u, k), -x(v, k)] for u, v in graph.edges() for k in range(3))
    clauses.extend(CardEnc.atmost([x(v, 0) for v in range(30)], bound=9,
                                 vpool=pool, encoding=EncType.seqcounter).clauses)
    with Cadical153(bootstrap_with=clauses) as solver:
        solver.conf_budget(2000)
        result = solver.solve_limited()
        if result is not True:
            return None
        model = set(solver.get_model())
        colors = [next(k for k in range(3) if x(v, k) in model) for v in range(30)]
    assert colors.count(0) <= 9
    assert all(colors[u] != colors[v] for u, v in graph.edges())
    return colors


def theta_child(connection):
    """One long-lived numerical process; parent terminates slow/hung solves."""
    import cvxpy as cp
    import numpy as np
    problems = {}
    c5_solver = None
    while True:
        request = connection.recv()
        if request is None:
            return
        n, edges = request
        try:
            reduction_fallback = None
            if n == 30:
                try:
                    from overnight30_c5 import C5ThetaSolver, support_shifts
                    if support_shifts(edges) is not None:
                        if c5_solver is None:
                            c5_solver = C5ThetaSolver()
                        reduced = c5_solver.solve(edges)
                        if "error" not in reduced:
                            connection.send(reduced)
                            continue
                        reduction_fallback = reduced["error"]
                except Exception as error:
                    reduction_fallback = f"{type(error).__name__}: {error}"
            if n not in problems:
                B = cp.Variable((n, n), symmetric=True)
                mask = cp.Parameter((n, n))
                problem = cp.Problem(cp.Maximize(cp.sum(B)),
                                     [B >> 0, cp.trace(B) == 1, cp.multiply(mask, B) == 0])
                problems[n] = B, mask, problem
            B, mask, problem = problems[n]
            zeros = np.ones((n, n))
            np.fill_diagonal(zeros, 0)
            for u, v in edges:
                zeros[u, v] = zeros[v, u] = 0
            mask.value = zeros
            attempts = []
            for solver in ("CLARABEL", "SCS"):
                try:
                    if solver == "CLARABEL":
                        problem.solve(solver=solver, max_iter=250)
                    else:
                        problem.solve(solver=solver, eps=1e-6, max_iters=20000)
                    if problem.value is None or not np.isfinite(problem.value) or B.value is None:
                        attempts.append(f"{solver}: missing finite primal result ({problem.status})")
                        continue
                    matrix = np.asarray(B.value)
                    spectrum = np.linalg.eigvalsh(matrix)
                    trace_residual = float(abs(np.trace(matrix)-1))
                    zero_residual = float(np.max(abs(zeros*matrix)))
                    if (problem.status not in ("optimal", "optimal_inaccurate")
                            or spectrum[0] < -2e-5 or trace_residual > 2e-5
                            or zero_residual > 2e-5 or problem.value > 3.0001):
                        attempts.append(f"{solver}: unreliable primal residuals/status")
                        continue
                    break
                except Exception as error:
                    attempts.append(f"{solver}: {type(error).__name__}: {error}")
            else:
                raise RuntimeError("; ".join(attempts) or str(problem.status))
            connection.send({"theta_numeric": float(problem.value), "solver": solver,
                             "solver_status": problem.status,
                             "trace_residual": trace_residual,
                             "zero_residual": zero_residual,
                             "min_eigenvalue": float(spectrum[0]),
                             "primal_rank_numeric_1e-7": int(sum(spectrum > 1e-7)),
                             "primal_diagonal": matrix.diagonal().tolist(),
                             "primal_matrix": matrix.tolist(),
                             "reduction_fallback_reason": reduction_fallback})
        except Exception as error:
            connection.send({"error": f"{type(error).__name__}: {error}"})


class ThetaProcess:
    def __init__(self):
        self.process = self.connection = None

    def stop(self):
        if self.process is not None:
            if self.process.is_alive():
                self.process.terminate()
            self.process.join(timeout=0.25)
            if self.process.is_alive():
                self.process.kill()
                self.process.join(timeout=0.25)
            self.connection.close()
        self.process = self.connection = None

    def solve(self, n, edges, timeout):
        if self.process is None or not self.process.is_alive():
            self.stop()
            context = mp.get_context("spawn")
            self.connection, child = context.Pipe()
            self.process = context.Process(target=theta_child, args=(child,), daemon=True)
            self.process.start()
            child.close()
        try:
            self.connection.send((n, edges))
        except (BrokenPipeError, EOFError, OSError):
            self.stop()
            return {"error": "numerical subprocess exited before request"}
        if self.connection.poll(timeout):
            try:
                return self.connection.recv()
            except EOFError:
                self.stop()
                return {"error": "numerical subprocess exited"}
        self.stop()
        return {"error": "numerical solve timeout"}


def controls():
    import networkx as nx
    from verify_all import tricayley_z2z6
    H33 = nx.read_graph6(Path(__file__).with_name("h33.g6"))
    S = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)),
         ((0, 1), (1, 0), (1, 1), (1, 3)))
    vertices, edges = tricayley_z2z6(S)
    index = {v: i for i, v in enumerate(vertices)}
    process = ThetaProcess()
    try:
        result = {}
        for name, n, es in (("H33", 33, list(H33.edges())),
                            ("H36", 36, [(index[u], index[v]) for u, v in edges])):
            answer = process.solve(n, es, 30)
            if "error" in answer or abs(answer["theta_numeric"]-3) > 1e-6:
                raise RuntimeError(f"positive control failed: {name}: {answer}")
            result[name] = {k: v for k, v in answer.items() if k not in ("primal_matrix", "primal_diagonal")}
        return result
    finally:
        process.stop()


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value, indent=2)+"\n")
    temporary.replace(path)


def repair_tail(path):
    """Recover one interrupted final JSONL append while retaining its bytes."""
    if not path.exists():
        return
    with path.open("r+b") as stream:
        while True:
            offset = stream.tell()
            line = stream.readline()
            if not line:
                return
            try:
                json.loads(line)
            except (json.JSONDecodeError, UnicodeDecodeError) as error:
                if stream.read(1):
                    raise RuntimeError(f"invalid nonfinal checkpoint record: {path}") from error
                path.with_name(path.name+f".interrupted-{time.time_ns()}").write_bytes(line)
                stream.truncate(offset)
                return


def records(path):
    # Normally exactly one exists. Reading archive then tail also supports
    # continuing a directory after its earlier JSONL was compressed in place.
    for source in (path.with_suffix(path.suffix+".gz"), path):
        if not source.exists():
            continue
        stream = gzip.open(source, "rt") if source.suffix == ".gz" else source.open()
        with stream:
            for line_number, line in enumerate(stream, 1):
                try:
                    yield json.loads(line)
                except json.JSONDecodeError as error:
                    raise RuntimeError(f"invalid checkpoint line {source}:{line_number}") from error


def run(args):
    start = time.monotonic()
    deadline = start+args.seconds
    import networkx as nx
    output = args.output_dir/args.family/f"seed_{args.seed}"
    output.mkdir(parents=True, exist_ok=True)
    with (output/"worker.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        solver, edge_variables, orbits, permutation = build(args.family)
        num_orbits = len(orbits)
        positions = orbit_positions(orbits)
        state_path, models_path, cuts_path = (output/name for name in ("state.json", "models.jsonl", "cuts.jsonl"))
        state = {"family": args.family, "seed": args.seed,
                 "scope": "heuristic weak-symmetry support search; probabilities and overlaps free",
                 "context_sizes": [10, 10, 10], "generator": permutation, "edge_orbits": num_orbits,
                 "color_prune": not args.no_color_prune, "models": 0, "scored": 0,
                 "rejected_coloring": 0, "rejected_connectivity": 0, "rejected_twins": 0,
                 "rejected_contraction": 0,
                 "certified_dominance_cuts": 0,
                 "rejected_pair_bridge": 0, "rejected_degree_mass": 0,
                 "rejected_regular_pair": 0, "rejected_quotient": 0,
                 "rejected_probability_propagation": 0,
                 "certified_structural_dominance_cuts": 0,
                 "duplicate_graphs": 0, "solver_errors": 0, "batches": 0,
                 "duplicate_labelled": 0, "duplicate_isomorphic": 0,
                 "numerical_representatives_admitted": 0,
                 "deduplication_note": "models count labelled support assignments; bounded exact colored isomorphism checks do not certify a complete list of nonisomorphic classes",
                 "ancestor_dirs": [],
                 "packages": {p: version(p) for p in ("numpy", "networkx", "cvxpy", "clarabel", "python-sat")}}
        resume_path = None
        if args.resume_from is not None:
            resume_path = args.resume_from/args.family/f"seed_{args.seed}"
            if not (resume_path/"state.json").exists():
                raise ValueError(f"previous chunk has no matching state: {resume_path}")
            resume_path = resume_path.resolve()
        saved = None
        if state_path.exists():
            saved = json.loads(state_path.read_text())
        elif resume_path is not None:
            saved = json.loads((resume_path/"state.json").read_text())
            saved["ancestor_dirs"] = saved.get("ancestor_dirs", [])+[str(resume_path)]
        if saved is not None:
            for key in ("family", "seed", "generator", "edge_orbits", "color_prune"):
                if saved[key] != state[key]:
                    raise ValueError(f"resume configuration differs: {key}")
            state.update(saved)
        if args.family == "c5_double":
            from overnight30_quotient import quotient_patterns, verified_catalog
            patterns = quotient_patterns()
            state["base_quotient_filter"] = {"static_theorem": "Proposition46",
                                             "verified_farkas_certificates": len(verified_catalog()),
                                             "distinct_permuted_patterns": len(patterns),
                                             "count_predicate_values": sorted({value for pattern in patterns for value in pattern})}
        sources = [Path(path) for path in state["ancestor_dirs"]]+[output]
        if len({path.resolve() for path in sources}) != len(sources):
            raise ValueError("resume chain contains the same directory twice")
        repair_tail(models_path)
        repair_tail(cuts_path)
        # Exact labelled dedup is cheap. Within each WL bucket use a bounded
        # number of exact, context-preserving isomorphism checks; skipped checks
        # only result in extra SDP work, never an unsound exclusion.
        seen, buckets = set(), defaultdict(list)
        replayed = 0
        replayed_dominance = 0
        replayed_structural = 0
        statuses = Counter()
        for source in sources:
            for record in records(source/"models.jsonl"):
                bits = record.get("orbit_bits")
                if bits is None:
                    bits = orbit_bits_from_graph6(record["graph6"], positions)
                assert len(bits) == num_orbits
                solver.add_clause([-(i+1) if bit == "1" else i+1 for i, bit in enumerate(bits)])
                seen.add(record["graph6"])
                if record.get("wl_hash") and len(buckets[record["wl_hash"]]) < 20:
                    buckets[record["wl_hash"]].append(record["graph6"])
                replayed += 1
                statuses[record["status"]] += 1
            for cut in records(source/"cuts.jsonl"):
                if cut["reason"] == "same_context_pair_requires_length3_path":
                    regenerated = path_requirement(solver, edge_variables, cut["witness"])
                    # Base CNF improvements may change the first auxiliary ID.
                    # A path cut uses only its own auxiliaries and edge-orbit
                    # variables, so deterministic rebasing preserves the cut.
                    shift = regenerated["aux_start"]-cut["aux_start"]
                    rebased = []
                    for clause in cut["clauses"]:
                        converted = []
                        for literal in clause:
                            if abs(literal) <= num_orbits:
                                converted.append(literal)
                            else:
                                assert cut["aux_start"] <= abs(literal) <= cut["aux_end"]
                                converted.append((1 if literal > 0 else -1)*(abs(literal)+shift))
                        rebased.append(converted)
                    if regenerated["clauses"] != rebased:
                        raise ValueError("auxiliary path cut replay differs from saved CNF")
                    solver.append_formula(regenerated["clauses"])
                elif cut["reason"] == "verified_theta_dual_dominance":
                    solver.add_clause(replay_dominance(source, cut, orbits))
                    replayed_dominance += 1
                elif cut["reason"] == "verified_structural_dominance":
                    solver.add_clause(replay_structural_dominance(cut, orbits))
                    replayed_structural += 1
                else:
                    assert all(1 <= literal <= num_orbits for literal in cut["clause"])
                    solver.add_clause(cut["clause"])
        state.update(models=replayed, scored=statuses["scored_numerically"],
                     certified_dominance_cuts=replayed_dominance,
                     certified_structural_dominance_cuts=replayed_structural,
                     rejected_pair_bridge=statuses["pair_bridge_obstruction"],
                     rejected_degree_mass=statuses["degree_mass_obstruction"],
                     rejected_regular_pair=statuses["regular_pair_obstruction"],
                     rejected_quotient=statuses["quotient_obstruction"],
                     rejected_probability_propagation=statuses["probability_propagation_obstruction"],
                     rejected_coloring=statuses["unbalanced_coloring"],
                     rejected_connectivity=statuses["disconnected"], rejected_twins=statuses["false_twins"],
                     rejected_contraction=statuses["same_context_merge"],
                     duplicate_labelled=statuses["duplicate_labelled"],
                     duplicate_isomorphic=statuses["duplicate_isomorphic"],
                     duplicate_graphs=statuses["duplicate_labelled"]+statuses["duplicate_isomorphic"],
                     solver_errors=statuses["solver_error"],
                     numerical_representatives_admitted=statuses["scored_numerically"]+statuses["solver_error"])
        initial_models, initial_scored = state["models"], state["scored"]
        state["batches"] += 1
        state["termination"] = "running"
        rng = random.Random(args.seed+1000003*state["batches"])
        solver.set_phases([i if rng.random() < .8 else -i for i in range(1, num_orbits+1)])
        seed_path = args.seed_graphs or Path(__file__).with_name("overnight30_symmetry_seeds.json")
        seeds = json.loads(seed_path.read_text()).get(args.family, []) if seed_path.exists() else []
        seed_bits = [orbit_bits_from_graph6(item["graph6"], positions) for item in seeds]
        pending_seeds = [bits for bits, item in zip(seed_bits, seeds) if item["graph6"] not in seen]
        state["construction_seed_graphs"] = len(seed_bits)
        state["construction_seed_file"] = str(seed_path.resolve()) if seed_bits else None
        process = ThetaProcess()
        best_path = output/"best.json"
        earlier_bests = [json.loads((source/"best.json").read_text())
                         for source in sources if (source/"best.json").exists()]
        # New exact filters may invalidate an older numerical leader. Keep its
        # evidence explicitly, but report only currently surviving best graphs.
        from overnight30_screen import screen_graph
        from overnight30_quotient import find_obstruction
        admitted_bests, excluded_bests, checked_bests = [], [], set()
        for candidate in earlier_bests:
            if candidate["graph6"] in checked_bests:
                continue
            checked_bests.add(candidate["graph6"])
            local = find_obstruction(candidate)
            screening = ({"status": "excluded", "reason": "exact regular-pair or automorphism-quotient obstruction",
                          "certificate": local} if local is not None else screen_graph(candidate, sat_conflicts=2000))
            if screening["status"] == "excluded":
                excluded_bests.append(dict(candidate, screen=screening,
                                           exclusion_status="excluded by updated exact structural screen"))
            else:
                admitted_bests.append(candidate)
        if excluded_bests:
            atomic_json(output/"legacy_best_excluded.json",
                        {"scope": "historical numerical leaders excluded by current exact checks",
                         "candidates": excluded_bests})
        state["legacy_bests_excluded"] = len(excluded_bests)
        best = max(admitted_bests, key=lambda item: item["theta_numeric"]) if admitted_bests else None
        best_bits = orbit_bits_from_graph6(best["graph6"], positions) if best else None
        if best is not None:
            atomic_json(best_path, best)
        elif best_path.exists():
            best_path.unlink()  # the full excluded record was preserved above
        near_dir = output/"near3"
        near_dir.mkdir(exist_ok=True)
        dual_dir = output/"dual_cuts"
        dual_dir.mkdir(exist_ok=True)
        batch_dominance = batch_dual_attempts = 0
        dual_degree_signatures = set()
        with models_path.open("a", buffering=1) as model_stream, cuts_path.open("a", buffering=1) as cut_stream:
            def checkpoint():
                state["batch_elapsed_seconds"] = time.monotonic()-start
                state["updated_unix"] = time.time()
                state["best_theta_numeric"] = best["theta_numeric"] if best else None
                state["distinct_labelled_supports"] = len(seen)
                state["batch_models"] = state["models"]-initial_models
                state["batch_scored"] = state["scored"]-initial_scored
                state["batch_dominance_cuts"] = batch_dominance
                state["batch_dual_attempts"] = batch_dual_attempts
                atomic_json(state_path, state)

            def learn(clause, reason, witness):
                clause = sorted(set(clause))
                if not clause:
                    raise RuntimeError(f"empty {reason} cut: check its witness")
                solver.add_clause(clause)
                cut_stream.write(json.dumps({"clause": clause, "reason": reason, "witness": witness})+"\n")

            def learn_path(pair):
                cut = path_requirement(solver, edge_variables, pair)
                solver.append_formula(cut["clauses"])
                cut_stream.write(json.dumps(cut, separators=(",", ":"))+"\n")

            checkpoint()
            try:
                while time.monotonic() < deadline:
                    if (state["models"]-initial_models) % 25 == 0:
                        density_bias = rng.choice((.2, .4, .6, .8, .95))
                        solver.set_phases([i if rng.random() < density_bias else -i
                                           for i in range(1, num_orbits+1)])
                    solver.conf_budget(2000)
                    assumptions = []
                    if pending_seeds:
                        bits = pending_seeds.pop(0)
                        assumptions = [i if bit == "1" else -i for i, bit in enumerate(bits, 1)]
                    elif (seed_bits or best_bits) and rng.random() < .75:
                        bits = best_bits if best_bits and (not seed_bits or rng.random() < .75) else rng.choice(seed_bits)
                        assumptions = [i if bit == "1" else -i for i, bit in enumerate(bits, 1)
                                       if rng.random() < .65]
                        solver.set_phases([i if bit == "1" else -i for i, bit in enumerate(bits, 1)])
                    answer = solver.solve_limited(assumptions=assumptions)
                    if answer is False and assumptions:
                        solver.conf_budget(2000)
                        answer = solver.solve_limited()
                    if answer is None:
                        checkpoint()
                        continue
                    if not answer:
                        state["termination"] = "SAT family exhausted under recorded cuts; no proof log"
                        break
                    model = set(v for v in solver.get_model() if v > 0)
                    bits = "".join("1" if i in model else "0" for i in range(1, num_orbits+1))
                    edges = [e for i, orbit in enumerate(orbits, 1) if i in model for e in orbit]
                    graph = nx.Graph()
                    graph.add_nodes_from((v, {"part": v//10}) for v in range(30))
                    graph.add_edges_from(edges)
                    check_graph(graph, permutation)
                    g6 = nx.to_graph6_bytes(graph, header=False).decode().strip()
                    record = {"family": args.family, "seed": args.seed, "index": state["models"],
                              "graph6": g6, "edges_count": len(edges),
                              "context_sizes": [10, 10, 10]}
                    deadline_reached = False
                    state["models"] += 1
                    solver.add_clause([-i if i in model else i for i in range(1, num_orbits+1)])
                    if g6 in seen:
                        record["status"] = "duplicate_labelled"
                        state["duplicate_graphs"] += 1
                        state["duplicate_labelled"] += 1
                    elif not nx.is_connected(graph):
                        component = min(nx.connected_components(graph), key=len)
                        clause = [variable for (u, v), variable in edge_variables.items()
                                  if (u in component) != (v in component)]
                        learn(clause, "connected_30_event_GHZ_required", sorted(component))
                        record["status"] = "disconnected"
                        state["rejected_connectivity"] += 1
                    elif (merge_pair := missing_three_path(graph)) is not None:
                        learn_path(merge_pair)
                        record["status"] = "same_context_merge"
                        record["contraction_pair"] = list(merge_pair)
                        state["rejected_contraction"] += 1
                    elif ((structural_screen := worker_structural_screen(record))["status"] == "excluded"
                          and structural_screen["certificate"]["method"] in STRUCTURAL_OBSTRUCTIONS):
                        certificate = structural_screen["certificate"]
                        witness = {"graph6": g6, "method": certificate["method"],
                                   "witness": certificate["witness"]}
                        cut = {"clause": dominance_clause(certificate["graph"], orbits),
                               "reason": "verified_structural_dominance", "witness": witness}
                        replay_structural_dominance(cut, orbits)
                        learn(cut["clause"], cut["reason"], witness)
                        if certificate["method"] == "degree_mass":
                            record.update(status="degree_mass_obstruction", structural_method="degree_mass")
                            state["rejected_degree_mass"] += 1
                        elif certificate["method"] == "regular_pair":
                            record.update(status="regular_pair_obstruction", structural_method="regular_pair")
                            state["rejected_regular_pair"] += 1
                        elif certificate["method"] == "probability_propagation":
                            record.update(status="probability_propagation_obstruction", structural_method="probability_propagation")
                            state["rejected_probability_propagation"] += 1
                        elif certificate["method"] in {"c5_quotient", "orbit_probability_quotient"}:
                            record.update(status="quotient_obstruction", structural_method=certificate["method"])
                            state["rejected_quotient"] += 1
                        else:
                            record.update(status="pair_bridge_obstruction", structural_method=certificate["method"])
                            state["rejected_pair_bridge"] += 1
                        state["certified_structural_dominance_cuts"] += 1
                    else:
                        colors = bad_coloring(graph) if not args.no_color_prune else None
                        if colors is not None:
                            learn([variable for (u, v), variable in edge_variables.items() if colors[u] == colors[v]],
                                  "proper_3_coloring_has_class_at_most_9", colors)
                            record["status"] = "unbalanced_coloring"
                            state["rejected_coloring"] += 1
                        elif any(set(graph[u]) == set(graph[v]) for part in range(3)
                                 for u, v in itertools.combinations(range(10*part, 10*(part+1)), 2)):
                            record["status"] = "false_twins"
                            state["rejected_twins"] += 1
                        else:
                            wl = nx.weisfeiler_lehman_graph_hash(graph, node_attr="part", iterations=3)
                            record["wl_hash"] = wl
                            duplicate = False
                            for prior_g6 in buckets[wl]:
                                prior = nx.from_graph6_bytes(prior_g6.encode())
                                nx.set_node_attributes(prior, {v: v//10 for v in prior}, "part")
                                if nx.is_isomorphic(graph, prior, node_match=lambda a, b: a["part"] == b["part"]):
                                    duplicate = True
                                    break
                            if duplicate:
                                record["status"] = "duplicate_isomorphic"
                                state["duplicate_graphs"] += 1
                                state["duplicate_isomorphic"] += 1
                            else:
                                state["numerical_representatives_admitted"] += 1
                                if len(buckets[wl]) < 20:
                                    buckets[wl].append(g6)
                                remaining = deadline-time.monotonic()
                                if remaining <= .1:
                                    result = {"error": "deadline before numerical solve"}
                                    deadline_reached = True
                                else:
                                    result = process.solve(30, edges, min(args.solve_timeout, remaining))
                                if "error" in result:
                                    record.update(status="solver_error", error=result["error"])
                                    state["solver_errors"] += 1
                                else:
                                    state["scored"] += 1
                                    record.update({k: v for k, v in result.items()
                                                   if k in ("theta_numeric", "solver", "solver_status")})
                                    record["status"] = "scored_numerically"
                                    signature = tuple(sorted(dict(graph.degree()).values()))
                                    if (args.family == "c5_double" and batch_dominance < args.dual_cut_limit
                                            and batch_dual_attempts < 3*args.dual_cut_limit
                                            and "theta_dual_candidate_matrix" in result
                                            and deadline-time.monotonic() > .1
                                            and (batch_dual_attempts < 10 or best is None
                                                 or result["theta_numeric"] >= best["theta_numeric"]-.015
                                                 or signature not in dual_degree_signatures)):
                                        batch_dual_attempts += 1
                                        dual_degree_signatures.add(signature)
                                        try:
                                            from overnight30_fast_dual import certify_from_matrix, verify_certificate
                                            certificate = certify_from_matrix(edges, result["theta_dual_candidate_matrix"],
                                                                              result["theta_dual_upper_numeric"])
                                            if certificate is not None:
                                                bound = verify_certificate(certificate)
                                                raw = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
                                                digest = hashlib.sha256(raw).hexdigest()
                                                certificate_file = dual_dir/f"{digest}.json.gz"
                                                temporary = certificate_file.with_suffix(".tmp")
                                                temporary.write_bytes(gzip.compress(raw, compresslevel=9, mtime=0))
                                                temporary.replace(certificate_file)
                                                witness = {"graph6": g6, "sha256": digest, "upper_bound": str(bound),
                                                           "certificate": str(certificate_file.relative_to(output))}
                                                learn(dominance_clause(certificate, orbits), "verified_theta_dual_dominance", witness)
                                                batch_dominance += 1
                                                state["certified_dominance_cuts"] += 1
                                                record.update(exact_theta_upper=str(bound), dominance_certificate_sha256=digest)
                                        except Exception as error:
                                            record["dual_certification_error"] = f"{type(error).__name__}: {error}"
                                    full = dict(record, edges=sorted([list(e) for e in edges]),
                                                parts=[list(range(10*k, 10*(k+1))) for k in range(3)],
                                                certificate_status="numerical only; exact certification required")
                                    full.update(result)
                                    if "exact_theta_upper" in record:
                                        full["certificate_status"] = "exact dual upper bound below three; numerical leader retained as construction guide"
                                    if best is None or result["theta_numeric"] > best["theta_numeric"]:
                                        best = full
                                        best_bits = bits
                                        atomic_json(best_path, best)
                                        print(args.family, "seed", args.seed, "best theta", best["theta_numeric"],
                                              "models", state["models"], "scored", state["scored"], flush=True)
                                    if result["theta_numeric"] >= args.near_threshold:
                                        digest = hashlib.sha256(g6.encode()).hexdigest()[:16]
                                        atomic_json(near_dir/f"{digest}.json", full)
                                        print("HIGH PRIORITY numerical theta candidate", near_dir/f"{digest}.json", flush=True)
                    seen.add(g6)
                    model_stream.write(json.dumps(record, separators=(",", ":"))+"\n")
                    checkpoint()
                    if deadline_reached:
                        state["termination"] = "batch time budget"
                        break
                else:
                    state["termination"] = "batch time budget"
            finally:
                process.stop()
                solver.delete()
                checkpoint()
        print(json.dumps({k: state[k] for k in ("family", "seed", "models", "scored", "duplicate_graphs",
                                              "rejected_coloring", "solver_errors", "best_theta_numeric", "termination")}),
              flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", choices=sorted(FAMILIES), default="c2_mixed")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--seconds", type=float, default=300)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).with_name("overnight30")/"pools")
    parser.add_argument("--resume-from", type=Path, help="previous output-dir root for this family and seed; archived .jsonl.gz is supported")
    parser.add_argument("--seed-graphs", type=Path, help="family-indexed graph6 seed file; defaults to bundled construction seeds when present")
    parser.add_argument("--dual-cut-limit", type=int, default=100, help="maximum new exact C5 dominance cuts per batch")
    parser.add_argument("--solve-timeout", type=float, default=15)
    parser.add_argument("--near-threshold", type=float, default=2.999)
    parser.add_argument("--no-color-prune", action="store_true")
    parser.add_argument("--controls", action="store_true")
    parser.add_argument("--list-families", action="store_true")
    args = parser.parse_args()
    if args.list_families:
        print("\n".join(sorted(FAMILIES)))
        return
    if args.controls:
        print(json.dumps(controls(), indent=2))
        return
    if args.seconds <= 0 or args.solve_timeout <= 0 or args.dual_cut_limit < 0:
        parser.error("time budgets must be positive and the dual-cut limit nonnegative")
    run(args)


if __name__ == "__main__":
    main()

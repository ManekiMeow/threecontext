"""Bounded C3 x C3 RTT search with two distinct quotient coordinates.

A is F3^2 plus a fixed point. B has three size-three orbits with quotient
coordinates (x,x,y); C uses (x,x,y) or (x,y,y). The two translations generate
the FULL group action. Each complete cross-context edge orbit is a SAT variable.
The CNF requires triangle-freeness, nonempty/proper cross-neighbourhoods, total
degree at most ten, and distinct rows in each directed pair of contexts.

By default independent translations of the six fibre coordinates normalize
their A-supports to {0} or {1,2}. Empty/full A-supports would give duplicate rows,
so this is a symmetry normalization, not a uniform-probability assumption.
Use --unnormalized to omit it. Additional exact screens operate on each model;
their certificates are saved independently of the SAT solver.

An exactly excluded graph H generates the positive OR of all absent edge-orbit
variables, excluding every support contained in H. Each such dominance clause
is bound to its independently verified exact certificate. Unexcluded models are
blocked only by their complete assignment. CNF, orbit definitions, assignments,
clauses, and individual certificates are retained. If SAT
reports exhaustion, a separate proof-producing solve is attempted on the final
CNF, but its DRAT proof is NOT checked here. Neither bounded failure nor an
unchecked solver exhaustion is reported as a nonexistence theorem.

Run once into an empty output directory, with one CPU and a total wall budget:
  OPENBLAS_NUM_THREADS=1 python search/overnight30_c3square_search.py \\
      --seconds 80 --output-dir /tmp/ghz30-c3square-two-kernel
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
import os
from pathlib import Path
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

PARTS = [list(range(10*k, 10*k+10)) for k in range(3)]


def atomic_json(path, data):
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(data, separators=(",", ":"))+"\n")
    temporary.replace(path)


def generators(mode):
    assert mode in ("aligned", "opposed")
    quotients = ((0, 0, 1), (0, 0, 1) if mode == "aligned" else (0, 1, 1))
    result = []
    for coordinate in range(2):
        permutation = list(range(30))
        for x, y in itertools.product(range(3), repeat=2):
            image = ((x+1) % 3, y) if coordinate == 0 else (x, (y+1) % 3)
            permutation[3*x+y] = 3*image[0]+image[1]
        for offset, labels in zip((10, 20), quotients):
            for orbit, label in enumerate(labels):
                for j in range(3):
                    permutation[offset+3*orbit+j] = offset+3*orbit+(j+(label == coordinate)) % 3
        result.append(permutation)
    assert all(result[0][result[1][v]] == result[1][result[0][v]] for v in range(30))
    assert all(p[p[p[v]]] == v for p in result for v in range(30))
    return result


def build(mode, normalize=True, orbit_bounds=False):
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    permutations = generators(mode)
    edge_variables, orbits = {}, []
    for u in range(30):
        for v in range(u+1, 30):
            if u//10 == v//10 or (u, v) in edge_variables:
                continue
            orbit, todo = {(u, v)}, [(u, v)]
            while todo:
                a, b = todo.pop()
                for permutation in permutations:
                    edge = tuple(sorted((permutation[a], permutation[b])))
                    if edge not in orbit:
                        orbit.add(edge)
                        todo.append(edge)
            orbits.append(sorted(orbit))
            edge_variables.update((edge, len(orbits)) for edge in orbit)
    assert len(edge_variables) == 300
    pool, clauses = IDPool(start_from=len(orbits)+1), []

    def edge(u, v):
        return edge_variables[tuple(sorted((u, v)))]

    for a, b, c in itertools.product(*PARTS):
        clauses.append(sorted({-edge(a, b), -edge(a, c), -edge(b, c)}))
    for u in range(30):
        incident = []
        for part in PARTS:
            if u in part:
                continue
            literals = [edge(u, v) for v in part]
            clauses.append(sorted(set(literals)))
            clauses.append(sorted({-literal for literal in literals}))
            incident.extend(literals)
        # Repeated orbit literals count separate incident edges, as in the
        # existing symmetry worker's sequential degree counter.
        clauses.extend(CardEnc.atmost(incident, bound=10, vpool=pool,
                                     encoding=EncType.seqcounter).clauses)
    for part in PARTS:
        for u, v in itertools.combinations(part, 2):
            for other in PARTS:
                if other is part:
                    continue
                differences = []
                pairs = {tuple(sorted((edge(u, w), edge(v, w)))) for w in other
                         if edge(u, w) != edge(v, w)}
                for a, b in sorted(pairs):
                    different = pool.id()
                    differences.append(different)
                    clauses.extend(([-a, -b, -different], [a, b, -different],
                                    [a, -b, different], [-a, b, different]))
                clauses.append(differences)
    if normalize:
        for offset in (10, 20):
            for orbit in range(3):
                zero, one, two = [edge(0, offset+3*orbit+j) for j in range(3)]
                # (zero,one,two) is exactly (1,0,0) or (0,1,1).
                clauses.extend(([zero, one], [-zero, -one], [-one, two], [one, -two]))
    if orbit_bounds:
        # Strict probability budgets for each transitive orbit, without
        # averaging event weights: a 9-orbit has degree pair <=(4,5), and a
        # 3-orbit has pair <=(3,6) or <=(4,5), with either degree at least two.
        for u in list(range(9))+list(range(10, 19))+list(range(20, 29)):
            incident = []
            for other in PARTS:
                if u in other:
                    continue
                literals = [edge(u, v) for v in other]
                incident.extend(literals)
                clauses.extend(CardEnc.atleast(literals, bound=2, vpool=pool,
                                               encoding=EncType.seqcounter).clauses)
                clauses.extend(CardEnc.atmost(literals, bound=5 if u < 9 else 6, vpool=pool,
                                              encoding=EncType.seqcounter).clauses)
            clauses.extend(CardEnc.atmost(incident, bound=9, vpool=pool,
                                          encoding=EncType.seqcounter).clauses)
        ab, ac = edge(9, 19), edge(9, 29)
        for offset in (10, 20):
            selected = [edge(9, offset+3*k) for k in range(3)]
            clauses.append([-ab, -ac]+selected)
            clauses.extend([-ab, -ac, -u, -v] for u, v in itertools.combinations(selected, 2))
        if mode == "opposed":
            # A mixed fixed-point pattern gives six events with only three
            # distinct degree-two neighbourhoods in the opposite singleton
            # kernel orbit. Distinct pair rows therefore force both ab and ac.
            clauses.extend(([ab], [ac], [-edge(19, 29)]))
    return clauses, pool.top, orbits, permutations


def model_graph(bits, orbits):
    import networkx as nx
    graph = nx.Graph()
    graph.add_nodes_from(range(30))
    graph.add_edges_from(edge for bit, orbit in zip(bits, orbits) if bit == "1" for edge in orbit)
    return graph


def validate_graph(graph, permutations):
    assert set(graph) == set(range(30))
    edges = {tuple(sorted(edge)) for edge in graph.edges()}
    assert all(u//10 != v//10 for u, v in edges)
    assert all(not set(graph[u]) & set(graph[v]) for u, v in edges)
    assert all(graph.degree(v) <= 10 for v in graph)
    for part in PARTS:
        for other in PARTS:
            if part is other:
                continue
            rows = [frozenset(graph[u]) & frozenset(other) for u in part]
            assert all(0 < len(row) < 10 for row in rows)
            assert len(set(rows)) == 10
    for permutation in permutations:
        assert {tuple(sorted((permutation[u], permutation[v]))) for u, v in edges} == edges


def write_cnf(path, variables, clauses):
    with path.open("w") as stream:
        stream.write(f"p cnf {variables} {len(clauses)}\n")
        for clause in clauses:
            stream.write(" ".join(map(str, clause))+" 0\n")


def enumerate_mode(mode, directory, deadline, normalize, max_models, orbit_bounds=False, dominance=True):
    import networkx as nx
    from pysat.solvers import Cadical153
    from overnight30_screen import fast_screen, verify_certificate
    directory.mkdir(parents=True, exist_ok=True)
    clauses, variables, orbits, permutations = build(mode, normalize, orbit_bounds)
    atomic_json(directory/"action.json", {"mode": mode, "parts": PARTS,
                "generators": permutations, "edge_orbits": orbits, "edge_variables": len(orbits),
                "normalized_fibres": normalize, "necessary_orbit_bounds": orbit_bounds,
                "structural_dominance": dominance})
    write_cnf(directory/"base.cnf", variables, clauses)
    counts, survivors, blocking = Counter(), [], []
    state = {"mode": mode, "edge_orbits": len(orbits), "cnf_variables": variables,
             "cnf_clauses": len(clauses), "termination": "time_budget", "models": 0,
             "scope": "bounded SAT enumeration; solver exhaustion alone is not an exclusion proof"}
    with Cadical153(bootstrap_with=clauses) as solver, (directory/"models.jsonl").open("w", buffering=1) as stream:
        while time.monotonic() < deadline and state["models"] < max_models:
            solver.conf_budget(2000)
            answer = solver.solve_limited()
            if answer is None:
                counts["sat_conflict_budget"] += 1
                continue
            if not answer:
                state["termination"] = "solver_reports_exhausted_unverified"
                break
            model = set(solver.get_model())
            bits = "".join("1" if i in model else "0" for i in range(1, len(orbits)+1))
            block = [-i if i in model else i for i in range(1, len(orbits)+1)]
            graph = model_graph(bits, orbits)
            validate_graph(graph, permutations)
            record = {"n": 30, "parts": PARTS, "edges": sorted(graph.edges())}
            screening = fast_screen(record)
            status = screening["status"]
            row = {"index": state["models"], "bits": bits,
                   "g6": nx.to_graph6_bytes(graph, header=False).decode().strip(),
                   "edges_count": len(graph.edges()), "screen_status": status,
                   "reason": screening.get("reason"),
                   "dimension19_eligible": screening.get("dimension19_eligible")}
            if status == "excluded":
                cert = screening["certificate"]
                assert verify_certificate(cert)
                name = f"certificate_{state['models']:06d}.json"
                atomic_json(directory/name, cert)
                row["certificate"] = name
                row["certificate_sha256"] = hashlib.sha256((directory/name).read_bytes()).hexdigest()
                counts[cert["method"]] += 1
                if dominance:
                    block = [i for i, bit in enumerate(bits, 1) if bit == "0"]
                    assert block and all(i not in model for i in block)
                    row["sat_cut"] = "excluded_supergraph_dominance"
                    counts["dominance_cuts"] += 1
            elif status == "not_excluded":
                counts["screen_survivor"] += 1
                survivors.append({**record, **row, "mode": mode})
            else:
                raise RuntimeError(f"unexpected screen result: {screening}")
            row.setdefault("sat_cut", "complete_assignment")
            row["sat_clause"] = block
            blocking.append(block)
            solver.add_clause(block)
            stream.write(json.dumps(row, separators=(",", ":"))+"\n")
            state["models"] += 1
        if state["models"] >= max_models:
            state["termination"] = "model_budget"
    final_clauses = clauses+blocking
    write_cnf(directory/"final.cnf", variables, final_clauses)
    if state["termination"] == "solver_reports_exhausted_unverified" and time.monotonic() < deadline:
        with Cadical153(bootstrap_with=final_clauses, with_proof=True) as proof_solver:
            proof_solver.conf_budget(10000)
            if proof_solver.solve_limited() is False:
                proof = proof_solver.get_proof()
                (directory/"exhaustion.drat").write_text("\n".join(proof)+"\n")
                state["proof_log"] = "exhaustion.drat"
                state["proof_status"] = "recorded, not independently verified"
    state["screen_counts"] = dict(counts)
    state["final_cnf_clauses"] = len(final_clauses)
    state["final_cnf_sha256"] = hashlib.sha256((directory/"final.cnf").read_bytes()).hexdigest()
    atomic_json(directory/"summary.json", state)
    print(json.dumps(state), flush=True)
    return state, survivors


def main(args):
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("use an empty output directory; this bounded experiment does not overwrite evidence")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    begin, summaries, survivors = time.monotonic(), [], []
    deadline = begin+args.seconds
    modes = ("aligned", "opposed") if args.mode == "both" else (args.mode,)
    enumeration_end = deadline-min(12, args.seconds/5) if args.theta_count else deadline
    for i, mode in enumerate(modes):
        end = time.monotonic()+max(0, enumeration_end-time.monotonic())/(len(modes)-i)
        summary, candidates = enumerate_mode(mode, args.output_dir/mode, end,
                                             not args.unnormalized, args.max_models,
                                             args.necessary_orbit_bounds, not args.no_dominance)
        summaries.append(summary)
        survivors.extend(candidates)
    numerical, controls = [], None
    if survivors and args.theta_count and time.monotonic() < deadline:
        from overnight30_symmetry import ThetaProcess
        process = ThetaProcess()
        try:
            import networkx as nx
            from verify_all import tricayley_z2z6
            h33 = nx.read_graph6(Path(__file__).with_name("h33.g6"))
            shifts = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)),
                      ((0, 1), (1, 0), (1, 1), (1, 3)))
            vertices, edges = tricayley_z2z6(shifts)
            index = {vertex: i for i, vertex in enumerate(vertices)}
            controls = {}
            for name, n, es in (("H33", 33, list(h33.edges())),
                                ("H36", 36, [(index[u], index[v]) for u, v in edges])):
                remaining = deadline-time.monotonic()
                if remaining <= 0:
                    raise RuntimeError("budget ended before numerical controls completed")
                answer = process.solve(n, es, min(4, remaining))
                if "error" in answer or abs(answer["theta_numeric"]-3) > 1e-6:
                    raise RuntimeError(f"positive control failed: {name}: {answer}")
                controls[name] = answer["theta_numeric"]
            # Edge count is only a cheap selection heuristic, never a theta bound.
            for candidate in sorted(survivors, key=lambda row: row["edges_count"], reverse=True)[:args.theta_count]:
                remaining = deadline-time.monotonic()
                if remaining < 0.2:
                    break
                answer = process.solve(30, candidate["edges"], min(4, remaining))
                result = {**candidate, **answer, "scope": "numerical candidate, not a realization certificate"}
                path = args.output_dir/f"theta_{len(numerical):02d}.json"
                atomic_json(path, result)
                numerical.append({"path": path.name, "mode": candidate["mode"],
                                  "theta_numeric": answer.get("theta_numeric"), "error": answer.get("error")})
        finally:
            process.stop()
    final = {"scope": "bounded two-kernel C3 x C3 RTT search; no nonexistence conclusion",
             "elapsed_seconds": time.monotonic()-begin, "modes": summaries,
             "surviving_models": len(survivors), "controls": controls, "numerical": numerical}
    atomic_json(args.output_dir/"summary.json", final)
    print(json.dumps(final), flush=True)


def certify_survivors(root, seconds=60):
    """Bind every surviving assignment to an exact dual via an explicit map."""
    import networkx as nx
    from overnight30_screen import exact_dual_certificate, verify_certificate
    begin = time.monotonic()
    directory = root/"aligned"
    rows = [json.loads(line) for line in (directory/"models.jsonl").read_text().splitlines()]
    action = json.loads((directory/"action.json").read_text())
    summary = json.loads((directory/"summary.json").read_text())
    target = root/"survivor_certificates"
    target.mkdir(exist_ok=True)
    representatives, bindings = [], {}
    node_match = nx.algorithms.isomorphism.categorical_node_match("part", -1)
    for row in rows:
        if row["screen_status"] != "not_excluded":
            continue
        graph = model_graph(row["bits"], action["edge_orbits"])
        nx.set_node_attributes(graph, {v: v//10 for v in graph}, "part")
        key = nx.weisfeiler_lehman_graph_hash(graph, node_attr="part", iterations=3)
        chosen = mapping = None
        for i, representative in enumerate(representatives):
            if representative["wl"] != key:
                continue
            matcher = nx.algorithms.isomorphism.GraphMatcher(graph, representative["graph"], node_match=node_match)
            if matcher.is_isomorphic():
                chosen, mapping = i, [matcher.mapping[v] for v in range(30)]
                break
        if chosen is None:
            chosen, mapping = len(representatives), list(range(30))
            representatives.append({"graph": graph, "wl": key, "index": row["index"]})
        assert sorted(mapping) == list(range(30))
        assert all(mapping[v]//10 == v//10 for v in range(30))
        assert {tuple(sorted((mapping[u], mapping[v]))) for u, v in graph.edges()} == {
            tuple(sorted(edge)) for edge in representatives[chosen]["graph"].edges()}
        bindings[row["index"]] = {"representative": chosen, "vertex_map_to_representative": mapping}
    print(json.dumps({"survivors": len(bindings), "context_preserving_isomorphism_classes": len(representatives)}), flush=True)
    records = []
    for i, representative in enumerate(representatives):
        if time.monotonic()-begin >= seconds:
            raise TimeoutError("exact survivor certification budget ended; completed certificates remain resumable")
        path = target/f"representative_{i:03d}.json"
        graph = {"n": 30, "parts": PARTS,
                 "edges": [list(edge) for edge in sorted(representative["graph"].edges())]}
        numeric = None
        if path.exists():
            certificate = json.loads(path.read_text())
            assert {tuple(edge) for edge in certificate["graph"]["edges"]} == {tuple(edge) for edge in graph["edges"]}
            assert verify_certificate(certificate)
        else:
            certificate, numeric = exact_dual_certificate(graph)
            if certificate is None:
                raise RuntimeError(f"no exact dual found for representative {i}, numerical theta {numeric}")
            assert verify_certificate(certificate)
            atomic_json(path, certificate)
        record = {"representative": i, "source_assignment_index": representative["index"],
                  "certificate": str(path.relative_to(root)),
                  "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "bound": certificate["witness"]["bound"], "theta_numeric": numeric}
        records.append(record)
        atomic_json(root/"survivor_progress.json", {"representatives": records, "bindings": bindings,
                                                    "complete": len(records) == len(representatives)})
        print(json.dumps(record), flush=True)
    assignments = []
    for row in rows:
        entry = {"assignment_index": row["index"],
                 "final_cnf_clause_number": summary["cnf_clauses"]+row["index"]+1,
                 "bits": row["bits"], "g6": row["g6"], "sat_cut": row["sat_cut"],
                 "sat_clause": row["sat_clause"]}
        if row["screen_status"] == "excluded":
            path = directory/row["certificate"]
            entry.update({"certificate": str(path.relative_to(root)),
                          "sha256": row["certificate_sha256"]})
        else:
            binding = bindings[row["index"]]
            representative = records[binding["representative"]]
            entry.update({**binding, "certificate": representative["certificate"],
                          "sha256": representative["sha256"], "bound": representative["bound"]})
        assignments.append(entry)
    files = {}
    for name in ("action.json", "base.cnf", "final.cnf", "exhaustion.drat", "exhaustion.lrat", "models.jsonl"):
        path = directory/name
        if path.exists():
            files[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {"format": "c3square-aligned-coverage-v1", "scope": "exact assignment certificates and SAT-coverage evidence; mathematical base constraints require the accompanying reduction",
                "base_cnf_clauses": summary["cnf_clauses"], "cnf_variables": summary["cnf_variables"],
                "assignment_count": len(assignments), "structural_assignments": len(rows)-len(bindings),
                "surviving_assignments": len(bindings), "representative_count": len(records),
                "maximum_survivor_theta_bound": str(max(Fraction(row["bound"]) for row in records)),
                "files_sha256": files, "representatives": records, "assignments": assignments,
                "proofs": {"cnf": "aligned/final.cnf", "drat": "aligned/exhaustion.drat", "lrat": "aligned/exhaustion.lrat"},
                "elapsed_seconds": time.monotonic()-begin}
    atomic_json(root/"coverage_manifest.json", manifest)
    print(json.dumps({key: manifest[key] for key in ("assignment_count", "surviving_assignments", "representative_count", "maximum_survivor_theta_bound", "elapsed_seconds")}), flush=True)


def verify_coverage(root, rebuild_base=False):
    """Replay hashes, exact graph/dual witnesses, isomorphisms and the LRAT proof."""
    from overnight30_screen import verify_certificate
    from event30_lrat import read_cnf, verify as verify_lrat
    manifest = json.loads((root/"coverage_manifest.json").read_text())
    assert manifest["format"] == "c3square-aligned-coverage-v1"
    for filename, digest in manifest["files_sha256"].items():
        assert hashlib.sha256((root/filename).read_bytes()).hexdigest() == digest
    action = json.loads((root/"aligned/action.json").read_text())
    assert action["mode"] == "aligned" and action["parts"] == PARTS
    assert action["normalized_fibres"] is True and action["necessary_orbit_bounds"] is True
    assert action["generators"] == generators("aligned")
    covered = set()
    for saved_orbit in action["edge_orbits"]:
        saved = {tuple(edge) for edge in saved_orbit}
        assert saved and not saved & covered
        first = min(saved)
        closure, todo = {first}, [first]
        while todo:
            u, v = todo.pop()
            for permutation in action["generators"]:
                edge = tuple(sorted((permutation[u], permutation[v])))
                if edge not in closure:
                    closure.add(edge)
                    todo.append(edge)
        assert closure == saved
        covered.update(saved)
    assert covered == {(u, v) for u in range(30) for v in range(u+1, 30) if u//10 != v//10}
    assert action["edge_variables"] == len(action["edge_orbits"]) == 54
    variables, clauses = read_cnf(root/manifest["proofs"]["cnf"])
    base_variables, base = read_cnf(root/"aligned/base.cnf")
    assert variables == base_variables == manifest["cnf_variables"]
    assert len(base) == manifest["base_cnf_clauses"]
    assert all(clauses[i] == clause for i, clause in base.items())
    if rebuild_base:
        rebuilt, rebuilt_variables, rebuilt_orbits, rebuilt_generators = build(
            "aligned", normalize=True, orbit_bounds=True)
        assert rebuilt_variables == base_variables
        assert [frozenset(clause) for clause in rebuilt] == list(base.values())
        assert [[list(edge) for edge in orbit] for orbit in rebuilt_orbits] == action["edge_orbits"]
        assert rebuilt_generators == action["generators"]
    assert len(manifest["assignments"]) == manifest["assignment_count"] == len(clauses)-len(base)
    assert {entry["final_cnf_clause_number"] for entry in manifest["assignments"]} == set(range(len(base)+1, len(clauses)+1))
    verified = {}
    for entry in manifest["assignments"]:
        path = root/entry["certificate"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
        if entry["certificate"] not in verified:
            cert = json.loads(path.read_text())
            assert cert["method"] in {"theta_dual", "disconnected", "false_twins", "same_context_merge",
                                       "unbalanced_coloring", "bridge_shortage", "unbridged_vertex", "pair_hall",
                                       "degree_mass", "probability_propagation", "disjoint_partition", "one_sided_regular"}
            assert verify_certificate(cert)
            verified[entry["certificate"]] = cert
        cert = verified[entry["certificate"]]
        assert len(entry["bits"]) == len(action["edge_orbits"])
        assert set(entry["bits"]) <= {"0", "1"}
        assert entry["final_cnf_clause_number"] == len(base)+entry["assignment_index"]+1
        edges = {tuple(edge) for bit, orbit in zip(entry["bits"], action["edge_orbits"])
                 if bit == "1" for edge in orbit}
        if "vertex_map_to_representative" in entry:
            mapping = entry["vertex_map_to_representative"]
            assert sorted(mapping) == list(range(30))
            assert all(mapping[v]//10 == v//10 for v in range(30))
            edges = {tuple(sorted((mapping[u], mapping[v]))) for u, v in edges}
            assert cert["method"] == "theta_dual"
            assert Fraction(entry["bound"]) == Fraction(cert["witness"]["bound"]) < 3
        assert edges == {tuple(edge) for edge in cert["graph"]["edges"]}
        expected = ([i for i, bit in enumerate(entry["bits"], 1) if bit == "0"]
                    if entry["sat_cut"] == "excluded_supergraph_dominance" else
                    [-i if bit == "1" else i for i, bit in enumerate(entry["bits"], 1)])
        assert expected == entry["sat_clause"]
        assert frozenset(expected) == clauses[entry["final_cnf_clause_number"]]
    proof = verify_lrat(root/manifest["proofs"]["cnf"], root/manifest["proofs"]["lrat"])
    print(json.dumps({"verified_assignments": len(manifest["assignments"]),
                      "verified_distinct_certificates": len(verified), "lrat": proof,
                      "base_cnf_regenerated": rebuild_base,
                      "maximum_survivor_theta_bound": manifest["maximum_survivor_theta_bound"]}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=float, default=80)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--certify-survivors", type=Path)
    parser.add_argument("--verify-coverage", type=Path)
    parser.add_argument("--rebuild-base", action="store_true",
                        help="also regenerate the reviewed base CNF using python-sat")
    parser.add_argument("--mode", choices=("both", "aligned", "opposed"), default="both")
    parser.add_argument("--max-models", type=int, default=100000)
    parser.add_argument("--theta-count", type=int, default=6)
    parser.add_argument("--unnormalized", action="store_true")
    parser.add_argument("--necessary-orbit-bounds", action="store_true",
                        help="add exact orbit-degree and fixed-pattern necessities from the RTT reduction")
    parser.add_argument("--no-dominance", action="store_true",
                        help="block one model at a time instead of certified supergraph dominance cuts")
    args = parser.parse_args()
    if not 0 < args.seconds <= 90 or args.max_models < 1 or args.theta_count < 0:
        parser.error("require 0<seconds<=90, max-models>=1, theta-count>=0")
    if args.certify_survivors:
        certify_survivors(args.certify_survivors, args.seconds)
    elif args.verify_coverage:
        verify_coverage(args.verify_coverage, args.rebuild_base)
    elif args.output_dir:
        main(args)
    else:
        parser.error("provide --output-dir, --certify-survivors, or --verify-coverage")

"""Exact screening and certificate mining for (10,10,10) GHZ candidates.

Fast Python API:
    fast_screen(record)
    screen_graph(record, sat_conflicts=20000, exact_dual=False, try_primal=False)

A record has ``edges`` and ``parts``, or ``graph6``/``g6`` (the latter defaults
to the three contiguous size-ten parts). All graph vertices are integers.

CLI:
    python search/overnight30_screen.py scan POOL --output-dir CERTS --seconds 60
    python search/overnight30_screen.py scan POOL --output-dir CERTS --exact-dual --min-theta 2.85
    python search/overnight30_screen.py verify CERT.json
    python search/overnight30_screen.py self-test

Verification of saved certificates uses only standard Python. SAT and SDP are
discovery tools, never trusted as numeric infeasibility certificates. The
existing rational-primal helper requires a strict feasible margin and may fail
on a boundary face; that failure is explicitly inconclusive.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import time


def normalize_graph(record):
    if "graph" in record and isinstance(record["graph"], dict):
        record = record["graph"]
    if "graph6" in record or "g6" in record:
        import networkx as nx
        encoded = record.get("graph6", record.get("g6"))
        if not isinstance(encoded, bytes):
            encoded = encoded.encode("ascii")
        graph = nx.from_graph6_bytes(encoded.strip())
        n, edges = len(graph), list(graph.edges())
    elif "B_lines" in record and "C_lines" in record:
        B, C = list(map(set, record["B_lines"])), list(map(set, record["C_lines"]))
        if len(B) != 10 or len(C) != 10:
            raise ValueError("line contexts must each have ten events")
        n, edges = 30, []
        for offset, lines in ((10, B), (20, C)):
            edges.extend((p, offset+j) for j, line in enumerate(lines) for p in line)
        edges.extend((10+j, 20+k) for j, k in itertools.product(range(10), repeat=2)
                     if B[j].isdisjoint(C[k]))
    else:
        edges = record["edges"]
        if record.get("parts"):
            n = sum(map(len, record["parts"]))
        else:
            n = int(record.get("n", 30))
    parts = record.get("parts")
    if parts is None:
        parts = [list(range(10*k, 10*k+10)) for k in range(3)] if n == 30 else []
    clean_edges = set()
    for u, v in edges:
        if not isinstance(u, int) or not isinstance(v, int) or not 0 <= u < n or not 0 <= v < n or u == v:
            raise ValueError("edges must join distinct integer vertices in range(n)")
        clean_edges.add(tuple(sorted((u, v))))
    return {"n": n, "parts": [sorted(part) for part in parts],
            "edges": [list(edge) for edge in sorted(clean_edges)]}


def graph_id(graph):
    payload = json.dumps(graph, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()[:20]


def adjacency(graph):
    adj = [set() for _ in range(graph["n"])]
    for u, v in graph["edges"]:
        adj[u].add(v)
        adj[v].add(u)
    return adj


def components(adj, vertices=None):
    unseen = set(range(len(adj))) if vertices is None else set(vertices)
    found = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        component, stack = {root}, [root]
        while stack:
            u = stack.pop()
            new = adj[u] & unseen
            unseen.difference_update(new)
            component.update(new)
            stack.extend(new)
        found.append(sorted(component))
    return found


def context_coloring(graph):
    parts = graph["parts"]
    if len(parts) != 3 or sorted(v for p in parts for v in p) != list(range(graph["n"])):
        raise ValueError("parts must partition all graph vertices into three contexts")
    return {v: k for k, part in enumerate(parts) for v in part}


def hall_deficiency(adj, left, right):
    """Return a Hall-deficient left subset, or None if a matching covers it."""
    right = set(right)
    mate = {}

    def augment(u, visited):
        for v in sorted(adj[u] & right):
            if v in visited:
                continue
            visited.add(v)
            if v not in mate or augment(mate[v], visited):
                mate[v] = u
                return True
        return False

    for u in left:
        augment(u, set())
    if len(mate) == len(left):
        return None
    reached_left = set(left)-set(mate.values())
    reached_right, stack = set(), list(reached_left)
    while stack:
        u = stack.pop()
        for v in (adj[u] & right)-reached_right:
            reached_right.add(v)
            if v in mate and mate[v] not in reached_left:
                reached_left.add(mate[v])
                stack.append(mate[v])
    assert len(reached_right) < len(reached_left)
    return sorted(reached_left)


def base_validation(graph, require30=True):
    colors = context_coloring(graph)
    if require30 and (graph["n"] != 30 or list(map(len, graph["parts"])) != [10, 10, 10]):
        raise ValueError("screening theorem applies to exactly (10,10,10)")
    adj = adjacency(graph)
    if any(colors[u] == colors[v] for u, v in graph["edges"]):
        raise ValueError("a prescribed edge lies inside a proposed context")
    if any(adj[u] & adj[v] for u, v in graph["edges"]):
        raise ValueError("candidate is not triangle-free")
    return adj, colors


def certificate(graph, method, **witness):
    return {"format": "overnight30-certificate-v1", "graph": graph,
            "method": method, "witness": witness}


def pair_bridge_data(graph, pair, adj=None):
    """Pair components and third-context vertices touching both component sides.

    This computes support geometry only; the cardinality exclusion below is
    specific to thirty events. In particular this helper also accepts H33.
    """
    if adj is None:
        adj = adjacency(graph)
    k, l = pair
    left, right = set(graph["parts"][k]), set(graph["parts"][l])
    third = graph["parts"][3-k-l]
    records = []
    for component in components(adj, left | right):
        X, Y = set(component) & left, set(component) & right
        bridges = [v for v in third if adj[v] & X and adj[v] & Y]
        records.append({"component": component, "left": sorted(X), "right": sorted(Y),
                        "bridges": sorted(bridges)})
    return records


def context_degree_mass(graph, context, adj=None):
    """Exact sum of the strict event-probability upper bounds in a context."""
    from event30_orbit_degree import relaxed_event_upper
    if adj is None:
        adj = adjacency(graph)
    others = [set(part) for k, part in enumerate(graph["parts"]) if k != context]
    total = Q(0)
    for vertex in graph["parts"][context]:
        q, r = (len(adj[vertex] & part) for part in others)
        assert 1 <= q <= 9 and 1 <= r <= 9
        total += relaxed_event_upper(q, r)
    return total


def fast_screen(record):
    """Cheap exact obstructions. 'not_excluded' never claims feasibility."""
    graph = normalize_graph(record)
    result = {"graph_id": graph_id(graph), "status": "not_excluded", "graph": graph}
    if graph["n"] != 30 or list(map(len, graph["parts"])) != [10, 10, 10]:
        result.update(status="out_of_scope", reason="30-event structural screen only")
        return result
    try:
        adj, colors = base_validation(graph)
    except ValueError as exc:
        result.update(status="invalid_candidate", reason=str(exc))
        return result
    cc = components(adj)
    if len(cc) > 1:
        cert = certificate(graph, "disconnected", component=cc[0])
        result.update(status="excluded", reason="disconnected graph has a GHZ component smaller than 30",
                      certificate=cert)
        return result
    for part in graph["parts"]:
        for u, v in itertools.combinations(part, 2):
            if adj[u] == adj[v]:
                cert = certificate(graph, "false_twins", vertices=[u, v])
                result.update(status="excluded", reason="false twins collapse to 29 events",
                              certificate=cert)
                return result
    for part in graph["parts"]:
        for u, v in itertools.combinations(part, 2):
            if not any(adj[a] & adj[v] for a in adj[u]):
                cert = certificate(graph, "same_context_merge", vertices=[u, v])
                result.update(status="excluded", reason="a same-context pair has no length-three path and merges to 29 events",
                              certificate=cert)
                return result
    pair_diagnostics = []
    for k, l in itertools.combinations(range(3), 2):
        pair = graph["parts"][k] + graph["parts"][l]
        cc = components(adj, pair)
        hall = hall_deficiency(adj, graph["parts"][k], graph["parts"][l])
        pair_diagnostics.append({"pair": [k, l], "components": len(cc), "perfect_matching": hall is None})
        for comp in cc:
            if sum(colors[v] == k for v in comp) != sum(colors[v] == l for v in comp):
                alternate = [colors[v] for v in range(30)]
                for v in comp:
                    alternate[v] = l if colors[v] == k else k
                cert = certificate(graph, "unbalanced_coloring", coloring=alternate,
                                   discovery="unbalanced Kempe component")
                result.update(status="excluded", reason="Kempe swap gives a context of at most nine events",
                              certificate=cert, pair_diagnostics=pair_diagnostics)
                return result
        bridge_data = pair_bridge_data(graph, [k, l], adj)
        pair_diagnostics[-1]["bridge_sizes"] = [len(data["bridges"]) for data in bridge_data]
        for data in bridge_data:
            if len(data["bridges"]) < len(data["left"]):
                cert = certificate(graph, "bridge_shortage", pair=[k, l],
                                   component=data["component"], bridges=data["bridges"])
                result.update(status="excluded", reason="a pair component has too few third-context bridges",
                              certificate=cert, pair_diagnostics=pair_diagnostics)
                return result
        covered = set().union(*(set(data["bridges"]) for data in bridge_data))
        missing = set(graph["parts"][3-k-l])-covered
        if missing:
            cert = certificate(graph, "unbridged_vertex", pair=[k, l], vertex=min(missing))
            result.update(status="excluded", reason="a third-context event bridges no pair component",
                          certificate=cert, pair_diagnostics=pair_diagnostics)
            return result
        if hall is not None:
            cert = certificate(graph, "pair_hall", pair=[k, l], deficient_subset=hall)
            result.update(status="excluded", reason="two-context support violates the global probability Hall bound",
                          certificate=cert, pair_diagnostics=pair_diagnostics)
            return result
        if len(cc) > 1 and "dimension19_certificate" not in result:
            result["dimension19_certificate"] = certificate(
                graph, "pair_disconnected_d19", pair=[k, l], component=cc[0])
    result["pair_diagnostics"] = pair_diagnostics
    result["dimension19_pair_connected"] = all(x["components"] == 1 for x in pair_diagnostics)
    result["dimension19_pair_perfect_matching"] = all(x["perfect_matching"] for x in pair_diagnostics)
    result["dimension19_eligible"] = "dimension19_certificate" not in result
    from event30_one_sided_regular import find_obstruction as one_sided_obstruction
    regular = one_sided_obstruction(graph)
    if regular is not None:
        cert = certificate(graph, "one_sided_regular", **regular["witness"])
        result.update(status="excluded", reason="a one-sided regular context pair has degree outside three or four",
                      certificate=cert)
        return result
    # The balanced pair-component checks above exclude cross-degree zero.
    # Triangle-freeness then excludes cross-degree ten: any neighbor in the
    # third context would have cross-degree zero toward that complete side.
    bounds = [context_degree_mass(graph, k, adj) for k in range(3)]
    result["degree_mass_bounds"] = list(map(str, bounds))
    for context, bound in enumerate(bounds):
        if bound <= 1:
            cert = certificate(graph, "degree_mass", context=context, bound=str(bound))
            result.update(status="excluded", reason="strict degree probability bounds cannot sum to one in a context",
                          certificate=cert)
            return result
    # Iterated strict neighborhood-mass bounds strengthen the one-step degree
    # budget. Lazy import avoids a cycle with the standalone exact verifier.
    from event30_probability_propagation import screen as propagate_probabilities
    propagated = propagate_probabilities(graph, max_rounds=3)
    result["probability_propagation_rounds"] = propagated["rounds"]
    result["probability_propagation_budgets"] = propagated["context_budgets"]
    if propagated["status"] == "excluded":
        proof = propagated["certificate"]
        cert = certificate(graph, "probability_propagation",
                           rounding_denominator=proof["rounding_denominator"],
                           rounds=proof["rounds"], terminal=proof["witness"])
        result.update(status="excluded", reason="iterated strict probability bounds are inconsistent",
                      certificate=cert)
        return result
    from event30_disjoint_partition import screen as partition_probabilities
    partitioned = partition_probabilities(graph)
    result["disjoint_partition_budgets"] = partitioned["context_budgets"]
    if partitioned["status"] == "excluded":
        result.update(status="excluded", reason="disjoint-neighborhood partition bounds cannot sum to one",
                      certificate=partitioned["certificate"])
    return result


def unbalanced_coloring(graph, conflict_budget=20000):
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    from pysat.solvers import Cadical153
    pool = IDPool()
    var = lambda v, k: pool.id((v, k))
    clauses = []
    for v in range(30):
        clauses.append([var(v, k) for k in range(3)])
        clauses.extend([-var(v, k), -var(v, l)] for k, l in itertools.combinations(range(3), 2))
    clauses.extend([-var(u, k), -var(v, k)] for u, v in graph["edges"] for k in range(3))
    # By renaming colors, an unbalanced coloring has color zero of size <=9.
    clauses.extend(CardEnc.atmost([var(v, 0) for v in range(30)], bound=9,
                                 vpool=pool, encoding=EncType.seqcounter).clauses)
    with Cadical153(bootstrap_with=clauses) as solver:
        solver.conf_budget(conflict_budget)
        feasible = solver.solve_limited()
        if feasible is not True:
            return None, "unknown_budget" if feasible is None else "no_unbalanced_coloring_found"
        model = {x for x in solver.get_model() if x > 0}
        return [next(k for k in range(3) if var(v, k) in model) for v in range(30)], "witness"


def exact_psd(matrix):
    """Symmetric Schur complements over Fraction; handles zero pivots exactly."""
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        return False
    if any(matrix[i][j] != matrix[j][i] for i in range(n) for j in range(i)):
        return False
    matrix = [list(map(Q, row)) for row in matrix]
    for i in range(n):
        pivot = matrix[i][i]
        row = matrix[i][i+1:]
        if pivot < 0:
            return False
        if pivot == 0:
            if any(row):
                return False
            continue
        for j in range(i+1, n):
            for k in range(j, n):
                value = matrix[j][k] - row[j-i-1]*row[k-i-1]/pivot
                matrix[j][k] = matrix[k][j] = value
    return True


def exact_dual_certificate(graph, denominators=(10000, 1000000)):
    import networkx as nx
    from dualcert import dual_upper
    H = nx.Graph()
    H.add_nodes_from(range(graph["n"]))
    H.add_edges_from(graph["edges"])
    nodes, fixed, numeric_matrix, optimum = dual_upper(H)
    if nodes != list(range(graph["n"])) or not math.isfinite(optimum):
        raise RuntimeError("unexpected dual solver output")
    n = len(nodes)
    for denominator in denominators:
        rounded = [[Q(1) if fixed[i, j] else Q(round(float(numeric_matrix[i, j])*denominator), denominator)
                    for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(i):
                rounded[i][j] = rounded[j][i]
        for slack in (1e-5, 1e-4, 1e-3, 1e-2):
            bound = Q(math.ceil((optimum+slack)*denominator), denominator)
            if bound >= 3:
                continue
            if exact_psd([[(bound if i == j else 0)-rounded[i][j] for j in range(n)] for i in range(n)]):
                free = [int(rounded[i][j]*denominator) for i in range(n) for j in range(i+1, n)
                        if not fixed[i, j]]
                cert = certificate(graph, "theta_dual", bound=str(bound), denominator=denominator,
                                   free_numerators=free)
                verify_certificate(cert)
                return cert, optimum
    return None, optimum


def rational_primal_certificate(graph):
    import networkx as nx
    from gcert import rational_certificate
    colors = context_coloring(graph)
    H = nx.Graph()
    H.add_nodes_from(range(graph["n"]))
    H.add_edges_from(graph["edges"])
    output, margin = rational_certificate(H, colors)
    if output is None:
        return None, margin
    nodes, edges, values, weights = output
    assert nodes == list(range(graph["n"]))
    cert = certificate(graph, "theta_primal", weights=weights,
                       edge_values=[[u, v, value] for (u, v), value in zip(edges, values)])
    verify_certificate(cert)
    return cert, margin


def verify_certificate(cert):
    """No NetworkX/SAT/SDP dependencies: check the complete stored witness."""
    if not __debug__:
        raise RuntimeError("certificate verification requires assertions enabled; do not use python -O")
    if cert["format"] != "overnight30-certificate-v1":
        raise ValueError("unsupported certificate format")
    graph = cert["graph"]
    # Saved graphs always have explicit edges, so this remains standard Python.
    if normalize_graph(graph) != graph:
        raise ValueError("graph is not normalized")
    n, adj = graph["n"], adjacency(graph)
    method, witness = cert["method"], cert["witness"]
    if method in ("disconnected", "false_twins", "same_context_merge", "unbalanced_coloring", "bridge_shortage", "unbridged_vertex", "degree_mass", "probability_propagation", "disjoint_partition", "one_sided_regular", "pair_hall", "pair_disconnected_d19", "pair_hall_d19"):
        adj, colors = base_validation(graph)
        if method == "disconnected":
            S = set(witness["component"])
            assert S and len(S) < 30 and S <= set(range(30))
            assert all(not (adj[v]-S) for v in S)
        elif method == "false_twins":
            u, v = witness["vertices"]
            assert u != v and colors[u] == colors[v] and adj[u] == adj[v]
        elif method == "same_context_merge":
            u, v = witness["vertices"]
            assert u != v and colors[u] == colors[v]
            assert not any(adj[a] & adj[v] for a in adj[u])
        elif method == "unbalanced_coloring":
            c = witness["coloring"]
            assert len(c) == 30 and all(k in (0, 1, 2) for k in c)
            assert all(c[u] != c[v] for u, v in graph["edges"])
            assert min(c.count(k) for k in range(3)) <= 9
        elif method == "bridge_shortage":
            k, l = witness["pair"]
            assert 0 <= k < l < 3
            matches = [data for data in pair_bridge_data(graph, [k, l], adj)
                       if data["component"] == witness["component"]]
            assert len(matches) == 1
            data = matches[0]
            assert data["bridges"] == witness["bridges"]
            assert len(data["left"]) == len(data["right"]) > 0
            assert len(data["bridges"]) < len(data["left"])
        elif method == "unbridged_vertex":
            k, l = witness["pair"]
            assert 0 <= k < l < 3
            vertex = witness["vertex"]
            assert vertex in graph["parts"][3-k-l]
            assert all(vertex not in data["bridges"] for data in pair_bridge_data(graph, [k, l], adj))
        elif method == "degree_mass":
            context = witness["context"]
            assert type(context) is int and 0 <= context < 3
            bound = context_degree_mass(graph, context, adj)
            assert str(bound) == witness["bound"] and bound <= 1
        elif method == "probability_propagation":
            from event30_probability_propagation import verify_certificate as verify_propagation
            verify_propagation({"format": "ghz30-probability-propagation-v1", "graph": graph,
                                "rounding_denominator": witness["rounding_denominator"],
                                "rounds": witness["rounds"], "witness": witness["terminal"]})
        elif method == "disjoint_partition":
            from event30_disjoint_partition import verify_certificate as verify_partition
            verify_partition(cert)
        elif method == "one_sided_regular":
            from event30_one_sided_regular import verify_certificate as verify_regular
            verify_regular({"format": "ghz30-one-sided-regular-v1", "graph": graph,
                            "witness": witness})
        elif method == "pair_disconnected_d19":
            k, l = witness["pair"]
            assert 0 <= k < l < 3
            pair = set(graph["parts"][k]+graph["parts"][l])
            S = set(witness["component"])
            assert S and S < pair
            assert all(not (adj[v] & (pair-S)) for v in S)
        else:
            k, l = witness["pair"]
            assert 0 <= k < l < 3
            S = set(witness["deficient_subset"])
            assert S and S <= set(graph["parts"][k])
            neighbors = set().union(*(adj[v] for v in S)) & set(graph["parts"][l])
            assert len(neighbors) < len(S)
        return True
    if method == "theta_dual":
        bound, denominator = Q(witness["bound"]), witness["denominator"]
        assert bound < 3 and isinstance(denominator, int) and denominator > 0
        free = iter(witness["free_numerators"])
        matrix = [[Q(1) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if j not in adj[i]:
                    matrix[i][j] = matrix[j][i] = Q(next(free), denominator)
        marker = object()
        assert next(free, marker) is marker
        assert exact_psd([[(bound if i == j else 0)-matrix[i][j] for j in range(n)] for i in range(n)])
        return True
    if method == "theta_primal":
        adj, colors = base_validation(graph, require30=False)
        weights = list(map(Q, witness["weights"]))
        assert len(weights) == n and min(weights) > 0
        matrix = [[Q(0) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            matrix[i][i] = weights[i]
        seen = set()
        for u, v, value in witness["edge_values"]:
            edge = tuple(sorted((u, v)))
            assert u != v and v in adj[u] and edge not in seen
            seen.add(edge)
            matrix[u][v] = matrix[v][u] = Q(value)
        for i in range(n):
            for k, part in enumerate(graph["parts"]):
                if k != colors[i]:
                    assert sum(matrix[i][j] for j in part) == weights[i]
        assert exact_psd(matrix)
        assert sum(map(sum, matrix)) == 3*sum(weights)
        return True
    raise ValueError("unrecognized certificate method")


def screen_graph(record, sat_conflicts=20000, exact_dual=False, try_primal=False):
    result = fast_screen(record)
    graph = result["graph"]
    if result["status"] == "excluded":
        verify_certificate(result["certificate"])
        return result
    if result["status"] == "invalid_candidate":
        return result
    if result["status"] == "not_excluded" and sat_conflicts > 0:
        coloring, status = unbalanced_coloring(graph, sat_conflicts)
        result["coloring_search"] = status
        if coloring is not None:
            result["certificate"] = certificate(graph, "unbalanced_coloring", coloring=coloring,
                                                discovery="SAT witness")
            verify_certificate(result["certificate"])
            result.update(status="excluded", reason="proper coloring has a context of at most nine events")
            return result
    if exact_dual:
        try:
            cert, numeric = exact_dual_certificate(graph)
            result["theta_dual_numeric"] = numeric
            if cert is not None:
                result.update(status="excluded", reason="exact rational theta upper bound below three", certificate=cert)
                return result
        except Exception as exc:
            result["dual_discovery_error"] = f"{type(exc).__name__}: {exc}"
    if try_primal:
        try:
            cert, margin = rational_primal_certificate(graph)
            result["primal_margin_numeric"] = float(margin) if margin is not None else None
            if cert is not None:
                result.update(status="certified_ghz", certificate=cert)
            else:
                result["primal_note"] = "inconclusive: strict-margin/positive-weight-floor helper can miss boundary faces"
        except Exception as exc:
            result["primal_discovery_error"] = f"{type(exc).__name__}: {exc}"
    return result


def graph_records(data):
    if isinstance(data, dict):
        if any(key in data for key in ("edges", "graph6", "g6")) or all(key in data for key in ("B_lines", "C_lines")):
            yield data
        else:
            for value in data.values():
                yield from graph_records(value)
    elif isinstance(data, list):
        for value in data:
            yield from graph_records(value)


def source_records(paths):
    for path in paths:
        path = Path(path)
        if path.is_dir():
            endings = (".json", ".jsonl", ".g6", ".json.gz", ".jsonl.gz", ".g6.gz")
            yield from source_records(sorted(p for p in path.rglob("*") if p.name.endswith(endings)))
            continue
        suffix = Path(path.stem).suffix if path.suffix == ".gz" else path.suffix
        opener = gzip.open if path.suffix == ".gz" else open
        try:
            with opener(path, "rt") as stream:
                if suffix in (".jsonl", ".g6"):
                    for line in stream:
                        try:
                            if suffix == ".g6":
                                if line.strip():
                                    yield {"graph6": line.strip()}
                            else:
                                yield from graph_records(json.loads(line))
                        except json.JSONDecodeError:
                            continue  # includes a worker's incomplete final append
                else:
                    yield from graph_records(json.load(stream))
        except (OSError, EOFError, json.JSONDecodeError):
            continue  # a supervisor may still be archiving or replacing a pool


def numeric_theta(record):
    for key in ("theta_numeric", "best_theta_numeric", "theta", "score"):
        value = record.get(key)
        if isinstance(value, (int, float)) and math.isfinite(value):
            return float(value)
    return None


def scan(args):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    index = args.output_dir/"screened.jsonl"
    seen = set()
    if index.exists() and not args.recheck:
        with index.open() as stream:
            for line in stream:
                try:
                    seen.add(json.loads(line)["graph_id"])
                except (json.JSONDecodeError, KeyError):
                    continue
    start, counts, processed = time.monotonic(), Counter(), 0
    for record in source_records(args.inputs):
        if time.monotonic()-start >= args.seconds or processed >= args.limit:
            break
        try:
            graph = normalize_graph(record)
            gid = graph_id(graph)
            if gid in seen:
                continue
            estimate = numeric_theta(record)
            eligible = estimate is not None and estimate >= args.min_theta
            result = screen_graph(graph, sat_conflicts=args.sat_conflicts,
                                  exact_dual=args.exact_dual and eligible,
                                  try_primal=args.try_primal and eligible and estimate >= 3-1e-5)
            # If no input score is supplied, an explicit exact-dual request
            # still acts on the record when min-theta is disabled.
            if estimate is None and args.exact_dual and args.min_theta <= 0 and result["status"] not in ("excluded", "invalid_candidate"):
                result = screen_graph(graph, sat_conflicts=0, exact_dual=True, try_primal=args.try_primal)
            for key in ("certificate", "dimension19_certificate"):
                if key in result:
                    verify_certificate(result[key])
                    path = args.output_dir/f"{gid}_{key}.json"
                    path.write_text(json.dumps(result[key], indent=2)+"\n")
                    result[key] = path.name
            result.pop("graph", None)
            result["input_theta_numeric"] = estimate
            with index.open("a") as stream:
                stream.write(json.dumps(result, separators=(",", ":"))+"\n")
            seen.add(gid)
            counts[result["status"]] += 1
            processed += 1
            if result["status"] == "certified_ghz":
                print(json.dumps(result), flush=True)
                break
        except Exception as exc:
            counts[f"error:{type(exc).__name__}"] += 1
    summary = {"processed": processed, "counts": dict(counts), "seconds": time.monotonic()-start,
               "scope": "exact witnesses exclude individual graphs; all other statuses are inconclusive"}
    (args.output_dir/"last_scan.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary), flush=True)


def self_test():
    import networkx as nx
    import tempfile
    parts = [list(range(10*k, 10*k+10)) for k in range(3)]
    # A connected bipartite star has many false twins; the witness is exact.
    star = {"parts": parts, "edges": [[0, v] for v in range(10, 30)]+[[i, 10] for i in range(1, 10)]}
    result = fast_screen(star)
    assert result["status"] == "excluded"
    verify_certificate(result["certificate"])
    # Connected path with cycling context colors has no false twins, and an
    # unbalanced proper two-coloring supplied by the independent SAT step.
    order = [10*k+i for i in range(10) for k in range(3)]
    path = normalize_graph({"parts": parts, "edges": list(zip(order, order[1:]))})
    merged = fast_screen(path)
    assert merged["certificate"]["method"] == "same_context_merge"
    verify_certificate(merged["certificate"])
    bridge = pair_bridge_data(path, [0, 1])[0]
    assert bridge["left"] == [0] and bridge["right"] == [10] and not bridge["bridges"]
    verify_certificate(certificate(path, "bridge_shortage", pair=[0, 1],
                                   component=bridge["component"], bridges=[]))
    verify_certificate(certificate(path, "unbridged_vertex", pair=[0, 1], vertex=20))
    # A balanced blow-up of C6 has degree pair (5,5) at every vertex.
    clusters = [list(range(offset, offset+5)) for offset in (0, 10, 20, 5, 15, 25)]
    edges = [(u, v) for i in range(6) for u in clusters[i] for v in clusters[(i+1) % 6]]
    degree_control = normalize_graph({"parts": parts, "edges": edges})
    assert context_degree_mass(degree_control, 0) == Q(20, 27)
    verify_certificate(certificate(degree_control, "degree_mass", context=0, bound="20/27"))
    verify_certificate(certificate(degree_control, "one_sided_regular", event_context=0,
                                   regular_context=1, regular_degree=5))
    from event30_probability_propagation import screen as propagate_probabilities
    propagated = propagate_probabilities(degree_control)["certificate"]
    verify_certificate(certificate(degree_control, "probability_propagation",
                                   rounding_denominator=propagated["rounding_denominator"],
                                   rounds=propagated["rounds"], terminal=propagated["witness"]))
    # Equality is also excluded because every event upper bound is strict.
    edges = [(u, 10) for u in range(10)]
    edges += [(u, v) for u in range(6) for v in range(20, 26)]
    edges += [(u, 26) for u in range(6, 10)]
    degree_equal = normalize_graph({"parts": parts, "edges": edges})
    assert context_degree_mass(degree_equal, 0) == 1
    verify_certificate(certificate(degree_equal, "degree_mass", context=0, bound="1"))
    coloring, status = unbalanced_coloring(path)
    assert status == "witness"
    verify_certificate(certificate(path, "unbalanced_coloring", coloring=coloring))
    # Hall witnesses exclude every dimension at thirty events. Keep checking
    # the historical dimension19 method name for saved-certificate compatibility.
    adj = adjacency(path)
    deficient = hall_deficiency(adj, parts[0], parts[2])
    assert deficient is not None
    verify_certificate(certificate(path, "pair_hall", pair=[0, 2], deficient_subset=deficient))
    verify_certificate(certificate(path, "pair_hall_d19", pair=[0, 2], deficient_subset=deficient))
    control = nx.Graph()
    control.add_nodes_from(range(30))
    control.add_edges_from(path["edges"])
    encoded = nx.to_graph6_bytes(control, header=False).decode().strip()
    assert normalize_graph({"g6": encoded}) == path
    with tempfile.TemporaryDirectory() as temporary:
        pool = Path(temporary)/"pool.jsonl.gz"
        with gzip.open(pool, "wt") as stream:
            stream.write(json.dumps({"g6": encoded})+"\n")
            stream.write(json.dumps({"graph6": encoded})+"\n{incomplete")
        records = list(source_records([Path(temporary)]))
        assert len(records) == 2
        assert all(normalize_graph(record) == path for record in records)
    # H33 is outside the 30-event screen and passes a boundary-rank primal
    # certificate even though the strict-margin discovery helper may miss it.
    H = nx.read_graph6(Path(__file__).with_name("h33.g6"))
    graph = normalize_graph({"parts": [list(range(9)), list(range(9, 21)), list(range(21, 33))],
                             "edges": list(H.edges())})
    assert fast_screen(graph)["status"] == "out_of_scope"
    for k, l in itertools.combinations(range(3), 2):
        data = pair_bridge_data(graph, [k, l])
        assert all(item["bridges"] for item in data)
        assert set().union(*(set(item["bridges"]) for item in data)) == set(graph["parts"][3-k-l])
    weights = [str(Q(1, 9) if i < 9 else Q(1, 12)) for i in range(33)]
    values = [[u, v, str(Q(1, 36) if u < 9 or v < 9 else Q(1, 24))] for u, v in H.edges()]
    verify_certificate(certificate(graph, "theta_primal", weights=weights, edge_values=values))
    # Exact dual control, independent of an optimizer: C5 has theta<3;
    # D=I+A(C5) has maximum eigenvalue3 and is not a strict certificate.
    # Instead use a bipartite path and D=2I-ss^T, with bound2.
    small = normalize_graph({"n": 4, "parts": [], "edges": [[0, 1], [1, 2], [2, 3]]})
    signs = [1, -1, 1, -1]
    free = [-signs[i]*signs[j] for i in range(4) for j in range(i+1, 4) if [i, j] not in small["edges"]]
    verify_certificate(certificate(small, "theta_dual", bound="2", denominator=1, free_numerators=free))
    print("screen controls passed: false twins, same-context merge, pair bridges, strict degree mass, one-sided regularity, probability propagation, unbalanced coloring, global Hall and legacy witness, graph6, gzip pool, H33 boundary primal, exact dual")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)
    collect = sub.add_parser("scan")
    collect.add_argument("inputs", nargs="+", type=Path)
    collect.add_argument("--output-dir", type=Path, required=True)
    collect.add_argument("--seconds", type=float, default=60)
    collect.add_argument("--limit", type=int, default=500)
    collect.add_argument("--sat-conflicts", type=int, default=20000)
    collect.add_argument("--exact-dual", action="store_true")
    collect.add_argument("--try-primal", action="store_true")
    collect.add_argument("--min-theta", type=float, default=2.85)
    collect.add_argument("--recheck", action="store_true", help="revisit existing records, e.g. to upgrade to exact dual certificates")
    check = sub.add_parser("verify")
    check.add_argument("certificate", type=Path)
    sub.add_parser("self-test")
    args = parser.parse_args()
    if args.mode == "scan":
        scan(args)
    elif args.mode == "verify":
        assert verify_certificate(json.loads(args.certificate.read_text()))
        print("exact certificate verified")
    else:
        self_test()


if __name__ == "__main__":
    main()

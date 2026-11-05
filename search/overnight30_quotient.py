"""Exact worker filters for regular supports and verified orbit quotients.

The quotient exclusion requires the actual simultaneous C5 automorphism on
six consecutive five-vertex orbits, not just an equitable partition. Its SAT
encoding defines exact block-count predicates, then excludes every context
permutation and every within-context orbit swap of Proposition46 and the
independently verified probability, saturation, moment, and finite-phase catalogs.

Certificate recognition/verification uses only standard Python. The optional
--self-test also checks the CNF with python-sat. No numerical theta is used.
"""
import argparse
from functools import lru_cache
import gzip
import hashlib
import itertools
import json
from pathlib import Path

from event30_c5_quotient import Q, check_arithmetic

METHODS = {"regular_pair", "c5_quotient", "orbit_probability_quotient", "quotient_saturation", "c5_phase_quotient", "orbit_disjoint_probability_quotient"}
PAIRS = tuple((a, b) for a in range(6) for b in range(a+1, 6) if a//2 != b//2)
ALLOWED_REGULAR_PAIRS = {(3, 3), (3, 4), (4, 3)}
CATALOG = Path(__file__).with_name("overnight30_quotient_reviews")/"2026-11-04T2300"/"summary.json"
SATURATION_CATALOG = Path(__file__).with_name("overnight30_saturation_reviews")/"2026-11-04T2300"/"summary.json"
PHASE_CATALOG = Path(__file__).with_name("overnight30_phase_reviews")/"2026-11-05T0000"/"summary.json"
GENERAL_PHASE_CATALOG = Path(__file__).with_name("overnight30_general_phase_reviews")/"2026-11-05T0000"/"summary.json"
DISJOINT_CATALOG = Path(__file__).with_name("overnight30_disjoint_reviews")/"2026-11-05T0000"/"summary.json"
EXPANSION_CATALOG = Path(__file__).with_name("overnight30_quotient_expansion")/"2026-11-05T0000"/"summary.json"
ROOT = Path(__file__).resolve().parent.parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


@lru_cache(maxsize=1)
def orbit_orders():
    return tuple(tuple(2*contexts[k]+(i ^ swaps[k]) for k in range(3) for i in range(2))
                 for contexts in itertools.permutations(range(3))
                 for swaps in itertools.product(range(2), repeat=3))


@lru_cache(maxsize=1)
def verified_catalog():
    """Load exactly the catalog's claimed certificates and verify each rational proof."""
    from event30_probability_lp import verify_certificate as verify_probability
    if not CATALOG.exists():
        return ()
    summary = json.loads(CATALOG.read_text())
    result = []
    for record in summary["results"]:
        if not record.get("exact_lp_excluded"):
            continue
        relative = Path(record["certificate"])
        require(not relative.is_absolute() and ".." not in relative.parts, "invalid quotient certificate path")
        certificate = json.loads((CATALOG.parent/relative).read_text())
        margin = verify_probability(certificate)
        require(certificate["sizes"] == [5]*6 and certificate["parts"] == [[0, 1], [2, 3], [4, 5]],
                "catalog certificate does not have six equal C5 orbits")
        require(certificate["quotient"] == record["quotient"] and str(margin) == record["margin"],
                "catalog entry differs from its exact certificate")
        result.append(certificate)
    require(len(result) == summary["certified"], "catalog certificate count differs")
    return tuple(result)


@lru_cache(maxsize=1)
def verified_saturation_catalog():
    """Verify the separate saturation catalog without changing the frozen LP data."""
    from event30_quotient_saturation import verify_certificate as verify_saturation
    if not SATURATION_CATALOG.exists():
        return ()
    summary = json.loads(SATURATION_CATALOG.read_text())
    source = Path(__file__).resolve().parent.parent/summary["source_summary"]
    require(hashlib.sha256(source.read_bytes()).hexdigest() == summary["source_summary_sha256"],
            "saturation catalog's source summary has changed")
    result = []
    for record in summary["results"]:
        relative = Path(record["certificate"])
        require(not relative.is_absolute() and ".." not in relative.parts, "invalid saturation certificate path")
        certificate = json.loads((SATURATION_CATALOG.parent/relative).read_text())
        terminal = verify_saturation(certificate)
        require(certificate["sizes"] == [5]*6 and certificate["parts"] == [[0, 1], [2, 3], [4, 5]],
                "saturation certificate does not have six equal C5 orbits")
        require(certificate["quotient"] == record["quotient"] and terminal == record["terminal"],
                "saturation catalog entry differs from its exact certificate")
        result.append(certificate)
    require(len(result) == summary["certified"], "saturation catalog certificate count differs")
    return tuple(result)


@lru_cache(maxsize=None)
def verified_phase_manifest(relative_path, digest):
    """Check every finite-phase proof once; the content hash binds its artifact."""
    relative = Path(relative_path)
    require(not relative.is_absolute() and ".." not in relative.parts, "invalid phase manifest path")
    path = ROOT/relative
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == digest, "phase manifest content hash differs")
    manifest = json.loads(raw)
    if manifest["format"] == "ghz30-c5-phase-enumeration-v1":
        from event30_c5_phases import verify as verify_phases
    else:
        require(manifest["format"] == "ghz30-c5-general-phase-enumeration-v1", "unknown phase enumeration format")
        from event30_c5_general_phases import verify as verify_phases
    report = verify_phases(path.parent)
    require(report["complete_exact_exclusion"], "incomplete phases do not exclude a quotient")
    return manifest["quotient"], report["uniform_exact_upper"], manifest["format"]


def verify_phase_reference(proof):
    formats = {"ghz30-c5-phase-catalog-entry-v1": "ghz30-c5-phase-enumeration-v1",
               "ghz30-c5-general-phase-catalog-entry-v1": "ghz30-c5-general-phase-enumeration-v1"}
    require(proof.get("format") in formats, "unknown phase reference")
    quotient, bound, kind = verified_phase_manifest(proof["manifest"], proof["manifest_sha256"])
    require(formats[proof["format"]] == kind, "phase reference names another enumeration method")
    require(proof["quotient"] == quotient and proof["uniform_exact_upper"] == bound,
            "phase reference differs from the verified complete proof")
    return quotient


@lru_cache(maxsize=1)
def verified_phase_catalog():
    if not PHASE_CATALOG.exists():
        return ()
    summary = json.loads(PHASE_CATALOG.read_text())
    result = []
    for proof in summary["results"]:
        verify_phase_reference(proof)
        result.append(proof)
    require(len(result) == summary["certified"], "phase catalog certificate count differs")
    return tuple(result)


@lru_cache(maxsize=1)
def verified_general_phase_catalog():
    if not GENERAL_PHASE_CATALOG.exists():
        return ()
    summary = json.loads(GENERAL_PHASE_CATALOG.read_text())
    result = []
    for proof in summary["results"]:
        verify_phase_reference(proof)
        result.append(proof)
    require(len(result) == summary["certified"], "general phase catalog certificate count differs")
    return tuple(result)


@lru_cache(maxsize=1)
def verified_disjoint_catalog():
    from event30_disjoint_probability import verify_quotient_certificate
    if not DISJOINT_CATALOG.exists():
        return ()
    summary = json.loads(DISJOINT_CATALOG.read_text())
    result = []
    for record in summary["results"]:
        relative = Path(record["certificate"])
        require(not relative.is_absolute() and ".." not in relative.parts, "invalid disjoint certificate path")
        proof = json.loads((DISJOINT_CATALOG.parent/relative).read_text())
        margin = verify_quotient_certificate(proof)
        require(proof["sizes"] == [5]*6 and proof["parts"] == [[0, 1], [2, 3], [4, 5]],
                "disjoint probability proof does not describe six equal C5 orbits")
        require(proof["quotient"] == record["quotient"] and str(margin) == record["margin"],
                "disjoint catalog differs from exact proof")
        result.append(proof)
    require(len(result) == summary["certified"], "disjoint catalog certificate count differs")
    return tuple(result)


@lru_cache(maxsize=1)
def verified_expansion_catalog():
    """Replay every claimed proof in the separate compressed expansion catalog."""
    from event30_disjoint_probability import verify_quotient_certificate as verify_disjoint
    from event30_probability_lp import verify_certificate as verify_probability
    from event30_quotient_saturation import verify_certificate as verify_saturation
    if not EXPANSION_CATALOG.exists():
        return ()
    folder = EXPANSION_CATALOG.parent
    summary = json.loads(EXPANSION_CATALOG.read_text())
    hashes = json.loads((folder/"artifact_sha256.json").read_text())
    names = (EXPANSION_CATALOG.name, summary["results_archive"], summary["proof_archive"])
    for name in names:
        relative = Path(name)
        require(not relative.is_absolute() and ".." not in relative.parts, "invalid expansion artifact path")
        require(hashlib.sha256((folder/relative).read_bytes()).hexdigest() == hashes[name],
                "expansion artifact content hash differs")
    with gzip.open(folder/summary["results_archive"], "rt") as stream:
        records = [json.loads(line) for line in stream]
    with gzip.open(folder/summary["proof_archive"], "rt") as stream:
        proofs = [json.loads(line) for line in stream]
    certificates = {item["name"]: item["certificate"] for item in proofs}
    require(len(certificates) == len(proofs) == summary["proof_count"], "duplicate or missing expansion proofs")
    require(len(records) == summary["reviewed"], "expansion review count differs")
    verifiers = {"probability": verify_probability, "saturation": verify_saturation,
                 "disjoint_probability": verify_disjoint}
    result = []
    for record in records:
        if not record["exact_quotient_excluded"]:
            continue
        proof = certificates[record["certificate"]]
        require(record["proof_kind"] in verifiers, "unknown expansion proof kind")
        verifiers[record["proof_kind"]](proof)
        require(proof["sizes"] == [5]*6 and proof["parts"] == [[0, 1], [2, 3], [4, 5]],
                "expansion proof does not describe six equal C5 orbits")
        require(proof["quotient"] == record["quotient"], "expansion record differs from exact proof")
        result.append(proof)
    require(len(result) == summary["certified"], "expansion certificate count differs")
    return tuple(result)


@lru_cache(maxsize=1)
def pattern_catalog():
    check_arithmetic()
    result = {}
    # Static Proposition46 has priority when a later certificate overlaps it.
    sources = [(Q, None)]+[(certificate["quotient"], certificate)
                          for certificate in verified_catalog()+verified_saturation_catalog()+verified_phase_catalog()+verified_disjoint_catalog()+verified_general_phase_catalog()+verified_expansion_catalog()]
    for matrix, certificate in sources:
        for order in orbit_orders():
            pattern = tuple(matrix[order[a]][order[b]] for a, b in PAIRS)
            inverse = tuple(order.index(i) for i in range(6))
            if pattern not in result:
                result[pattern] = (inverse, certificate)
    return result


def quotient_patterns(verified_matrices=None):
    """Return all excluded signatures; extra matrices must already be proved excluded.

    None uses the verified on-disk catalog; an explicit iterable is useful for
    independently verified future catalogs. Proposition46 is always included.
    """
    if verified_matrices is None:
        return tuple(sorted(pattern_catalog()))
    matrices = [Q]+list(verified_matrices)
    for matrix in matrices:
        require(len(matrix) == 6 and all(len(row) == 6 for row in matrix), "quotient must be six by six")
        require(all(type(matrix[a][b]) is int and 0 <= matrix[a][b] <= 5 and matrix[a][b] == matrix[b][a]
                    and (a//2 != b//2 or matrix[a][b] == 0) for a in range(6) for b in range(6)),
                "invalid equal-orbit quotient matrix")
    return tuple(sorted({tuple(matrix[order[a]][order[b]] for a, b in PAIRS)
                         for matrix in matrices for order in orbit_orders()}))


def quotient_exclusion_cnf(edge_variables, pool, verified_matrices=None):
    """Deterministic exact-count CNF; callers regenerate it as base constraints."""
    clauses, predicates, blocks = [], {}, {}
    patterns = quotient_patterns(verified_matrices)
    counts = sorted({count for pattern in patterns for count in pattern})
    for a, b in PAIRS:
        variables = [edge_variables[(5*a, 5*b+s)] for s in range(5)]
        require(len(set(variables)) == 5, "C5 quotient block must have five distinct shift variables")
        blocks[(a, b)] = variables
        for count in counts:
            equal = pool.id()
            predicates[(a, b, count)] = equal
            # For each assignment, the five-literal prefix is false exactly on
            # that assignment. The last literal therefore fixes the predicate
            # iff the number of selected shifts equals count.
            for mask in range(32):
                blocking = [-variable if mask & (1 << s) else variable
                            for s, variable in enumerate(variables)]
                clauses.append(blocking+[equal if mask.bit_count() == count else -equal])
    for pattern in patterns:
        clauses.append([-predicates[(a, b, count)] for (a, b), count in zip(PAIRS, pattern)])
    return clauses, {"patterns": len(patterns), "auxiliary_variables": len(predicates),
                     "clauses": len(clauses), "predicates": predicates, "blocks": blocks}


def graph_data(record):
    record = record.get("graph", record)
    parts = [list(range(10*k, 10*k+10)) for k in range(3)]
    require(record.get("n", 30) == 30 and record.get("parts", parts) == parts,
            "worker obstruction requires contiguous contexts (10,10,10)")
    if "edges" in record:
        edges = record["edges"]
    else:
        encoded = record.get("graph6", record.get("g6", "")).encode("ascii")
        require(len(encoded) == 74 and encoded[0] == 93 and all(63 <= x <= 126 for x in encoded),
                "invalid 30-vertex graph6 record")
        edges = []
        for v in range(1, 30):
            for u in range(v):
                bit = v*(v-1)//2+u
                if (encoded[1+bit//6]-63) & (1 << (5-bit % 6)):
                    edges.append([u, v])
    normalized = sorted({tuple(sorted(edge)) for edge in edges})
    require(all(type(u) is int and type(v) is int and 0 <= u < v < 30 and u//10 != v//10
                for u, v in normalized), "invalid support edge")
    return {"n": 30, "parts": parts, "edges": [list(edge) for edge in normalized]}


def adjacency(graph):
    adj = [set() for _ in range(30)]
    for u, v in graph["edges"]:
        adj[u].add(v)
        adj[v].add(u)
    require(all(not adj[u] & adj[v] for u, v in graph["edges"]), "support is not triangle-free")
    return adj


def regular_degrees(adj):
    result = {}
    for k, l in itertools.combinations(range(3), 2):
        degrees = {sum(v//10 == l for v in adj[u]) for u in range(10*k, 10*k+10)}
        degrees.update(sum(v//10 == k for v in adj[u]) for u in range(10*l, 10*l+10))
        result[k, l] = next(iter(degrees)) if len(degrees) == 1 else None
    return result


def quotient_counts(graph, adj):
    edges = {tuple(edge) for edge in graph["edges"]}
    action = lambda u: 5*(u//5)+(u+1) % 5
    if {tuple(sorted((action(u), action(v)))) for u, v in edges} != edges:
        return None
    counts = [[sum(v//5 == b for v in adj[5*a]) for b in range(6)] for a in range(6)]
    require(all(sum(v//5 == b for v in adj[u]) == counts[u//5][b]
                for u in range(30) for b in range(6)), "automorphism did not give equitable orbits")
    return counts


def find_obstruction(record):
    graph = graph_data(record)
    adj = adjacency(graph)
    degrees = regular_degrees(adj)
    for context in range(3):
        others = [k for k in range(3) if k != context]
        pair = tuple(degrees[tuple(sorted((context, k)))] for k in others)
        if None not in pair and pair not in ALLOWED_REGULAR_PAIRS:
            return {"format": "overnight30-worker-obstruction-v1", "graph": graph,
                    "method": "regular_pair", "witness": {"context": context,
                    "other_contexts": others, "degrees": list(pair)}}
    counts = quotient_counts(graph, adj)
    if counts is not None:
        entry = pattern_catalog().get(tuple(counts[a][b] for a, b in PAIRS))
        if entry is not None:
            order, proof = entry
            witness = {"orbit_order": list(order)}
            method = "c5_quotient"
            if proof is not None and proof["format"] == "ghz30-orbit-probability-farkas-v1":
                witness["probability_certificate"] = proof
                method = "orbit_probability_quotient"
            elif proof is not None and proof["format"] == "ghz30-quotient-saturation-v1":
                witness["saturation_certificate"] = proof
                method = "quotient_saturation"
            elif proof is not None and proof["format"] == "ghz30-disjoint-orbit-probability-farkas-v1":
                witness["disjoint_probability_certificate"] = proof
                method = "orbit_disjoint_probability_quotient"
            elif proof is not None:
                require(proof["format"] in {"ghz30-c5-phase-catalog-entry-v1", "ghz30-c5-general-phase-catalog-entry-v1"}, "unknown quotient proof format")
                witness["phase_certificate"] = proof
                method = "c5_phase_quotient"
            return {"format": "overnight30-worker-obstruction-v1", "graph": graph,
                    "method": method,
                    "witness": witness}
    return None


def verify_certificate(certificate):
    require(certificate.get("format") == "overnight30-worker-obstruction-v1", "unknown worker certificate format")
    graph = graph_data(certificate["graph"])
    adj = adjacency(graph)
    witness = certificate["witness"]
    if certificate["method"] == "regular_pair":
        context = witness["context"]
        require(type(context) is int and 0 <= context < 3, "invalid central context")
        others = [k for k in range(3) if k != context]
        require(witness["other_contexts"] == others, "incorrect regular-pair contexts")
        degrees = regular_degrees(adj)
        pair = tuple(degrees[tuple(sorted((context, k)))] for k in others)
        require(None not in pair and list(pair) == witness["degrees"], "claimed pair supports are not regular")
        require(pair not in ALLOWED_REGULAR_PAIRS, "regular degrees are not excluded by Theorem45")
    elif certificate["method"] in ("c5_quotient", "orbit_probability_quotient", "quotient_saturation", "c5_phase_quotient", "orbit_disjoint_probability_quotient"):
        if certificate["method"] == "c5_quotient":
            check_arithmetic()
            template = Q
        elif certificate["method"] == "orbit_probability_quotient":
            from event30_probability_lp import verify_certificate as verify_probability
            probability = witness["probability_certificate"]
            verify_probability(probability)
            require(probability["sizes"] == [5]*6 and probability["parts"] == [[0, 1], [2, 3], [4, 5]],
                    "probability certificate does not describe six equal C5 orbits")
            template = probability["quotient"]
        elif certificate["method"] == "quotient_saturation":
            from event30_quotient_saturation import verify_certificate as verify_saturation
            saturation = witness["saturation_certificate"]
            verify_saturation(saturation)
            require(saturation["sizes"] == [5]*6 and saturation["parts"] == [[0, 1], [2, 3], [4, 5]],
                    "saturation certificate does not describe six equal C5 orbits")
            template = saturation["quotient"]
        elif certificate["method"] == "orbit_disjoint_probability_quotient":
            from event30_disjoint_probability import verify_quotient_certificate
            disjoint = witness["disjoint_probability_certificate"]
            verify_quotient_certificate(disjoint)
            require(disjoint["sizes"] == [5]*6 and disjoint["parts"] == [[0, 1], [2, 3], [4, 5]],
                    "disjoint probability proof does not describe six equal C5 orbits")
            template = disjoint["quotient"]
        else:
            template = verify_phase_reference(witness["phase_certificate"])
        order = tuple(witness["orbit_order"])
        require(order in orbit_orders(), "orbit ordering does not preserve the context partition")
        counts = quotient_counts(graph, adj)
        require(counts is not None, "equitable counts alone do not supply the required C5 automorphism")
        require(all(counts[order[a]][order[b]] == template[a][b] for a, b in PAIRS), "quotient does not match its exact template")
    else:
        raise ValueError("unknown worker obstruction method")
    return True


def self_test():
    from pysat.formula import IDPool
    from pysat.solvers import Cadical153
    variables = {}
    for index, (a, b) in enumerate(PAIRS):
        for shift in range(5):
            for i in range(5):
                variables[(5*a+i, 5*b+(i+shift) % 5)] = 5*index+shift+1
    clauses, info = quotient_exclusion_cnf(variables, IDPool(start_from=61), verified_matrices=())
    again, _ = quotient_exclusion_cnf(variables, IDPool(start_from=61), verified_matrices=())
    require(clauses == again, "base CNF generation is nondeterministic")
    tests = 0
    with Cadical153(bootstrap_with=clauses) as solver:
        for (a, b, count), equal in info["predicates"].items():
            block = info["blocks"][(a, b)]
            for mask in range(32):
                selected = {v for s, v in enumerate(block) if mask & (1 << s)}
                assignment = [i if i in selected else -i for i in range(1, 61)]
                truth = equal if mask.bit_count() == count else -equal
                require(solver.solve(assumptions=assignment+[truth]), "correct exact-count predicate rejected")
                require(not solver.solve(assumptions=assignment+[-truth]), "incorrect exact-count predicate accepted")
                tests += 2
        for pattern in quotient_patterns(verified_matrices=()):
            def assignment(signature):
                selected = {5*i+s+1 for i, count in enumerate(signature) for s in range(count)}
                return [i if i in selected else -i for i in range(1, 61)]
            require(not solver.solve(assumptions=assignment(pattern)), "excluded quotient accepted")
            for i in range(12):
                nearby = list(pattern)
                nearby[i] += 1
                require(solver.solve(assumptions=assignment(nearby)), "a different quotient was spuriously excluded")
    def regular_fixture(q, r):
        shifts = {(0, 1): list(range(0, 2*q, 2)), (0, 2): list(range(0, 2*r, 2)),
                  (1, 2): [1, 3, 5]}
        edges = [[10*a+i, 10*b+(i+s) % 10] for (a, b), values in shifts.items()
                 for s in values for i in range(10)]
        return {"edges": edges}
    excluded = find_obstruction(regular_fixture(2, 4))
    require(excluded is not None and excluded["method"] == "regular_pair", "regular (2,4) was not excluded")
    verify_certificate(excluded)
    require(find_obstruction(regular_fixture(3, 3)) is None, "allowed regular (3,3) was excluded")
    print("Quotient CNF controls passed:", info["patterns"], "patterns,", info["auxiliary_variables"],
          "auxiliaries,", tests, "predicate checks; neighboring quotients and regular-pair controls passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    if args.self_test:
        self_test()
    elif args.check:
        result = find_obstruction(json.loads(args.check.read_text()))
        if result is not None:
            verify_certificate(result)
        print(json.dumps(result, indent=2))
    else:
        parser.error("provide --self-test or --check GRAPH.json")

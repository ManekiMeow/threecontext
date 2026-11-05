"""Complete phase coverage for C5-double quotients with a singleton spanning tree.

Vertices are (orbit, phase), orbit=0..5 and phase in Z5. For a<b, S_ab lists
the phase differences j-i of edges (a,i)--(b,j). Triangle-freeness is exactly
(S_ab+S_bc) intersect S_ac = empty for the eight context triples a<b<c.

If the singleton blocks contain a spanning tree, independent fiber shifts
uniquely make its five shifts zero, after fixing the shift of orbit0. Every
original support has one such normal form. Common multiplication by a unit
of Z5 then gives the representatives used here. Prefix enumeration rejects
only explicitly completed triangles, so its coverage is independent of SAT
or numerical optimization. Exact dual certificates verify every excluded
representative using integer arithmetic.

  python search/event30_c5_phases.py search --index 1 --seconds 60 --output-dir DIR
  python search/event30_c5_phases.py verify DIR

Verification and enumeration need only standard Python. Numerical discovery
uses one bounded subprocess. An unfinished run makes no class exclusion.
"""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import time

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
PAIRS = tuple((a, b) for a in range(6) for b in range(a+1, 6) if a//2 != b//2)
INDEX = {pair: i for i, pair in enumerate(PAIRS)}
TRIANGLES = tuple((INDEX[a, b], INDEX[b, c], INDEX[a, c])
                  for a in range(2) for b in range(2, 4) for c in range(4, 6))
UNITS = (1, 2, 3, 4)
MULTIPLY = {u: tuple(sum(1 << ((u*s) % 5) for s in range(5) if mask >> s & 1)
                    for mask in range(32)) for u in UNITS}
SUM = tuple(tuple(sum(1 << s for s in {(a+b) % 5 for a in range(5) for b in range(5)
                                     if left >> a & 1 and right >> b & 1})
                  for right in range(32)) for left in range(32))
SOURCE = Path(__file__).with_name("overnight30_quotient_reviews")/"2026-11-04T2300"/"summary.json"
ROOT = Path(__file__).resolve().parent.parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def singleton_tree(quotient):
    require(len(quotient) == 6 and all(len(row) == 6 for row in quotient), "expected six-orbit quotient")
    require(all(type(quotient[a][b]) is int and 0 <= quotient[a][b] <= 5
                and quotient[a][b] == quotient[b][a]
                and (a//2 != b//2 or quotient[a][b] == 0)
                for a in range(6) for b in range(6)), "invalid quotient")
    parent = list(range(6))
    def root(a):
        while parent[a] != a:
            a = parent[a]
        return a
    tree = []
    for a, b in PAIRS:
        if quotient[a][b] == 1 and root(a) != root(b):
            parent[root(a)] = root(b)
            tree.append((a, b))
    require(len(tree) == 5, "singleton blocks do not contain a spanning tree")
    return tuple(tree)


def canonical(masks):
    return min(tuple(MULTIPLY[u][mask] for mask in masks) for u in UNITS)


def normalize(masks, tree):
    shifts = {0: 0}
    while len(shifts) < 6:
        progress = False
        for a, b in tree:
            mask = masks[INDEX[a, b]]
            require(mask.bit_count() == 1, "tree block is not a singleton")
            phase = mask.bit_length()-1
            if a in shifts and b not in shifts:
                shifts[b] = (shifts[a]-phase) % 5
                progress = True
            elif b in shifts and a not in shifts:
                shifts[a] = (shifts[b]+phase) % 5
                progress = True
        require(progress or len(shifts) == 6, "disconnected gauge tree")
    return tuple(sum(1 << ((s+shifts[b]-shifts[a]) % 5) for s in range(5) if mask >> s & 1)
                 for (a, b), mask in zip(PAIRS, masks))


def enumerate_phases(quotient):
    tree = singleton_tree(quotient)
    fixed = {INDEX[pair] for pair in tree}
    domains = [(1,) if i in fixed else tuple(mask for mask in range(32) if mask.bit_count() == quotient[a][b])
               for i, (a, b) in enumerate(PAIRS)]
    assigned = [1 if i in fixed else None for i in range(12)]
    incident = [[triple for triple in TRIANGLES if i in triple] for i in range(12)]
    normal, prefixes = set(), 0
    def permitted(index, value):
        assigned[index] = value
        valid = all(any(assigned[i] is None for i in triple)
                    or not SUM[assigned[triple[0]]][assigned[triple[1]]] & assigned[triple[2]]
                    for triple in incident[index])
        assigned[index] = None
        return valid
    def visit():
        nonlocal prefixes
        prefixes += 1
        choices = [(i, tuple(value for value in domains[i] if permitted(i, value)))
                   for i in range(12) if assigned[i] is None]
        if not choices:
            normal.add(tuple(assigned))
            return
        index, values = min(choices, key=lambda item: (len(item[1]), item[0]))
        for value in values:
            assigned[index] = value
            visit()
        assigned[index] = None
    visit()
    representatives = sorted({canonical(masks) for masks in normal})
    orbit_sizes = Counter(len({tuple(MULTIPLY[u][mask] for mask in masks) for u in UNITS})
                          for masks in representatives)
    require(sum(size*count for size, count in orbit_sizes.items()) == len(normal), "unit-orbit coverage mismatch")
    raw = math.prod(math.comb(5, quotient[a][b]) for a, b in PAIRS)
    summary = {"gauge_tree": [list(pair) for pair in tree], "raw_phase_assignments": raw,
               "gauge_phase_assignments": math.prod(len(domain) for domain in domains),
               "normal_triangle_free_assignments": len(normal), "representatives": len(representatives),
               "raw_triangle_free_assignments": len(normal)*5**5,
               "unit_orbit_sizes": {str(k): v for k, v in sorted(orbit_sizes.items())},
               "prefixes_visited": prefixes,
               "representatives_sha256": hashlib.sha256(json.dumps(representatives, separators=(",", ":")).encode()).hexdigest()}
    return representatives, summary


def graph_edges(masks):
    return sorted((5*a+i, 5*b+(i+s) % 5) for (a, b), mask in zip(PAIRS, masks)
                  for s in range(5) if mask >> s & 1 for i in range(5))


def graph6(edges):
    edges = set(edges)
    bits = [int((u, v) in edges) for v in range(1, 30) for u in range(v)]
    bits += [0]*((-len(bits)) % 6)
    return "]"+"".join(chr(63+sum(bits[i+j] << (5-j) for j in range(6))) for i in range(0, len(bits), 6))


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value, indent=2)+"\n")
    temporary.replace(path)


def numeric_child(connection):
    import cvxpy as cp
    import numpy as np
    from overnight30_c5 import C5ThetaSolver
    reduced, full = C5ThetaSolver(), None
    while True:
        mode, edges = connection.recv()
        try:
            if mode == "reduced":
                result = reduced.solve(edges)
            else:
                if full is None:
                    D, upper = cp.Variable((30, 30), symmetric=True), cp.Variable()
                    fixed = cp.Parameter((30, 30))
                    problem = cp.Problem(cp.Minimize(upper), [cp.multiply(fixed, D) == fixed,
                                                            upper*np.eye(30)-D >> 0])
                    full = D, upper, fixed, problem
                D, upper, fixed, problem = full
                mask = np.eye(30)
                for u, v in edges:
                    mask[u, v] = mask[v, u] = 1
                fixed.value = mask
                problem.solve(solver="CLARABEL", max_iter=250)
                result = {"theta_dual_candidate_matrix": D.value.tolist(),
                          "theta_dual_upper_numeric": float(upper.value)}
            connection.send(result)
        except Exception as error:
            connection.send({"error": f"{type(error).__name__}: {error}"})


class NumericalWorker:
    def __init__(self):
        self.process = self.pipe = None
    def stop(self):
        if self.process is not None:
            self.process.terminate()
            self.process.join(.2)
            if self.process.is_alive():
                self.process.kill()
                self.process.join(.2)
            self.pipe.close()
        self.process = self.pipe = None
    def call(self, mode, edges, timeout):
        if timeout <= 0:
            return {"error": "batch deadline"}
        if self.process is None or not self.process.is_alive():
            self.stop()
            context = mp.get_context("spawn")
            self.pipe, child = context.Pipe()
            self.process = context.Process(target=numeric_child, args=(child,), daemon=True)
            self.process.start()
            child.close()
        try:
            self.pipe.send((mode, edges))
            if self.pipe.poll(timeout):
                return self.pipe.recv()
        except (BrokenPipeError, EOFError, OSError):
            pass
        self.stop()
        return {"error": "numerical process timeout or exit"}


def factor_certificate(edges, result):
    import numpy as np
    from overnight30_fast_dual import certify_from_matrix
    if "theta_dual_candidate_matrix" not in result:
        return None
    dual = np.asarray(result["theta_dual_candidate_matrix"], dtype=float)
    dual = (dual+dual.T)/2
    np.fill_diagonal(dual, 1)
    for u, v in edges:
        dual[u, v] = dual[v, u] = 1
    if not np.isfinite(dual).all():
        return None
    # This eigenvalue is only a discovery hint. The integer factor proof is
    # mandatory and independently checked before any exclusion is recorded.
    hint = float(np.linalg.eigvalsh(dual)[-1])
    return certify_from_matrix(edges, dual, hint)


def verify(directory):
    from overnight30_fast_dual import verify_certificate
    manifest = json.loads((directory/"manifest.json").read_text())
    require(manifest["format"] == "ghz30-c5-phase-enumeration-v1", "unknown phase manifest")
    source = Path(manifest["source_summary"])
    require(not source.is_absolute() and ".." not in source.parts, "source path must be repository relative")
    source_data = (ROOT/source).read_bytes()
    require(hashlib.sha256(source_data).hexdigest() == manifest["source_summary_sha256"], "source summary changed")
    index = manifest["source_index"]
    require(type(index) is int and 0 <= index < len(json.loads(source_data)["results"]), "invalid source index")
    require(json.loads(source_data)["results"][index]["quotient"] == manifest["quotient"], "source quotient differs")
    representatives, coverage = enumerate_phases(manifest["quotient"])
    require(coverage == manifest["coverage"], "coverage enumeration differs")
    expected = set(representatives)
    found, upper_bounds = set(), []
    for record in manifest["records"]:
        masks = tuple(record["phase_masks"])
        require(masks in expected and masks not in found, "unexpected or duplicate phase representative")
        edges = graph_edges(masks)
        require(graph6(edges) == record["graph6"], "phase representative graph differs")
        relative = Path(record["certificate"])
        require(not relative.is_absolute() and ".." not in relative.parts, "invalid certificate path")
        raw = gzip.decompress((directory/relative).read_bytes())
        require(hashlib.sha256(raw).hexdigest() == record["certificate_sha256"], "certificate hash differs")
        certificate = json.loads(raw)
        require(certificate["n"] == 30 and certificate["edges"] == [list(e) for e in edges], "certificate covers another graph")
        upper = verify_certificate(certificate)
        require(str(upper) == record["exact_theta_upper"], "recorded exact bound differs")
        upper_bounds.append(upper)
        found.add(masks)
    complete = found == expected
    require(manifest["complete_exact_exclusion"] == complete, "incorrect completeness claim")
    return {"complete_exact_exclusion": complete, "covered_representatives": len(found),
            "total_representatives": len(expected), "coverage": coverage,
            "uniform_exact_upper": str(max(upper_bounds)) if complete else None}


def search(args):
    from overnight30_fast_dual import verify_certificate
    start, deadline = time.monotonic(), time.monotonic()+args.seconds
    source = json.loads(args.source.read_text())
    chosen = source["results"][args.index]
    quotient = chosen["quotient"]
    representatives, coverage = enumerate_phases(quotient)
    directory = args.output_dir
    directory.mkdir(parents=True, exist_ok=True)
    (directory/"certificates").mkdir(exist_ok=True)
    manifest_path = directory/"manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        require(manifest["quotient"] == quotient and manifest["coverage"] == coverage, "resume quotient/coverage differs")
        verify(directory)
    else:
        manifest = {"format": "ghz30-c5-phase-enumeration-v1", "quotient": quotient,
                    "source_index": args.index, "source_summary": str(args.source.resolve().relative_to(ROOT)),
                    "source_summary_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
                    "coverage": coverage, "records": [], "complete_exact_exclusion": False,
                    "scope": "all C5 phases of this fixed quotient; representatives modulo fiber translations and common units"}
    completed = {tuple(row["phase_masks"]) for row in manifest["records"]}
    worker, pending = NumericalWorker(), []
    try:
        for masks in representatives:
            if masks in completed:
                continue
            if time.monotonic() >= deadline-.1:
                break
            edges = graph_edges(masks)
            numeric = worker.call("reduced", edges, min(8, deadline-time.monotonic()))
            certificate = factor_certificate(edges, numeric)
            if certificate is None and time.monotonic() < deadline-.1:
                fallback = worker.call("full_dual", edges, min(8, deadline-time.monotonic()))
                certificate = factor_certificate(edges, fallback)
            if certificate is None:
                pending.append({"phase_masks": list(masks), "graph6": graph6(edges),
                                "theta_numeric": numeric.get("theta_numeric"), "error": numeric.get("error")})
                if numeric.get("theta_numeric", 0) >= 2.999:
                    atomic_json(directory/f"near3_{len(pending)}.json", dict(numeric, phase_masks=list(masks), edges=edges))
                continue
            upper = verify_certificate(certificate)
            raw = json.dumps(certificate, sort_keys=True, separators=(",", ":")).encode()
            digest = hashlib.sha256(raw).hexdigest()
            relative = Path("certificates")/f"{digest}.json.gz"
            (directory/relative).write_bytes(gzip.compress(raw, compresslevel=9, mtime=0))
            manifest["records"].append({"phase_masks": list(masks), "graph6": graph6(edges),
                                         "theta_numeric": numeric.get("theta_numeric"),
                                         "exact_theta_upper": str(upper), "certificate": str(relative),
                                         "certificate_sha256": digest})
            completed.add(masks)
            manifest["complete_exact_exclusion"] = len(completed) == len(representatives)
            atomic_json(manifest_path, manifest)
            if len(completed) % 10 == 0:
                print("exactly covered", len(completed), "/", len(representatives), flush=True)
    finally:
        worker.stop()
        manifest["complete_exact_exclusion"] = len(completed) == len(representatives)
        manifest["last_discovery_seconds"] = time.monotonic()-start
        manifest["unresolved_attempts"] = pending
        atomic_json(manifest_path, manifest)
    report = verify(directory)
    report["discovery_seconds"] = manifest["last_discovery_seconds"]
    report["best_numeric_theta"] = max((row["theta_numeric"] for row in manifest["records"] if row["theta_numeric"] is not None), default=None)
    print(json.dumps(report), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("search")
    run.add_argument("--source", type=Path, default=SOURCE)
    run.add_argument("--index", type=int, default=1)
    run.add_argument("--seconds", type=float, default=60)
    run.add_argument("--output-dir", type=Path, required=True)
    check = sub.add_parser("verify")
    check.add_argument("directory", type=Path)
    args = parser.parse_args()
    if args.command == "verify":
        print(json.dumps(verify(args.directory), indent=2))
    else:
        require(args.seconds > 0, "positive discovery budget required")
        search(args)


if __name__ == "__main__":
    main()

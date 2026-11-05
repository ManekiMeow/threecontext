"""Complete C5 phase coverage using arbitrary proper-subset spanning trees.

Translation on a nonempty proper subset of Z5 is free: a nonzero stabilizer
would generate all of Z5. Choose the least translate on every tree edge.
There is a unique relative fiber translation producing those representatives.
This removes exactly 5**5 assignments. Common units are applied followed by
the same normalization. Completed triangle constraints are the only pruning.

The numerical worker and exact integer dual proofs are shared with the
singleton-tree search; this generalized coverage verifier is separate.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import time

from event30_c5_phases import (PAIRS, INDEX, TRIANGLES, MULTIPLY, SUM, UNITS,
                               SOURCE, ROOT, require, graph_edges, graph6,
                               atomic_json, NumericalWorker, factor_certificate)


def translate(mask, shift):
    return sum(1 << ((s+shift) % 5) for s in range(5) if mask >> s & 1)


TRANSLATION_REP = tuple(min(translate(mask, s) for s in range(5)) for mask in range(32))
TRANSLATION_SHIFT = tuple(next(s for s in range(5) if translate(mask, s) == TRANSLATION_REP[mask])
                          for mask in range(32))


def gauge_tree(quotient):
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
    result = []
    # Prefer singleton edges; their normalized phase domains are smallest.
    for a, b in sorted(PAIRS, key=lambda ab: (quotient[ab[0]][ab[1]] != 1, ab)):
        if 0 < quotient[a][b] < 5 and root(a) != root(b):
            parent[root(a)] = root(b)
            result.append((a, b))
    require(len(result) == 5, "proper subset blocks do not contain a spanning tree")
    return tuple(result)


def normalize(masks, tree):
    shifts = {0: 0}
    while len(shifts) < 6:
        progress = False
        for a, b in tree:
            mask = masks[INDEX[a, b]]
            require(0 < mask < 31, "gauge edge is not a nonempty proper subset")
            delta = TRANSLATION_SHIFT[mask]
            if a in shifts and b not in shifts:
                shifts[b] = (shifts[a]+delta) % 5
                progress = True
            elif b in shifts and a not in shifts:
                shifts[a] = (shifts[b]-delta) % 5
                progress = True
        require(progress or len(shifts) == 6, "disconnected gauge tree")
    return tuple(translate(mask, shifts[b]-shifts[a]) for (a, b), mask in zip(PAIRS, masks))


def unit_images(masks, tree):
    return {normalize(tuple(MULTIPLY[u][mask] for mask in masks), tree) for u in UNITS}


def enumerate_phases(quotient):
    tree = gauge_tree(quotient)
    fixed = {INDEX[pair] for pair in tree}
    domains = []
    for i, (a, b) in enumerate(PAIRS):
        domain = {mask for mask in range(32) if mask.bit_count() == quotient[a][b]}
        domains.append(tuple(sorted({TRANSLATION_REP[m] for m in domain} if i in fixed else domain)))
    assigned = [domain[0] if len(domain) == 1 else None for domain in domains]
    incident = [[triple for triple in TRIANGLES if i in triple] for i in range(12)]
    normal, prefixes = set(), 0
    def valid(triple):
        return any(assigned[i] is None for i in triple) or not SUM[assigned[triple[0]]][assigned[triple[1]]] & assigned[triple[2]]
    def permitted(index, value):
        assigned[index] = value
        result = all(valid(triple) for triple in incident[index])
        assigned[index] = None
        return result
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
    if all(valid(triple) for triple in TRIANGLES):
        visit()
    representatives = sorted({min(unit_images(masks, tree)) for masks in normal})
    orbit_sizes = Counter(len(unit_images(masks, tree)) for masks in representatives)
    require(sum(size*count for size, count in orbit_sizes.items()) == len(normal), "unit-orbit coverage mismatch")
    raw = math.prod(math.comb(5, quotient[a][b]) for a, b in PAIRS)
    summary = {"gauge_tree": [list(pair) for pair in tree], "raw_phase_assignments": raw,
               "gauge_phase_assignments": math.prod(len(domain) for domain in domains),
               "normal_triangle_free_assignments": len(normal), "representatives": len(representatives),
               "raw_triangle_free_assignments": len(normal)*5**5,
               "unit_orbit_sizes": {str(k): v for k, v in sorted(orbit_sizes.items())},
               "prefixes_visited": prefixes,
               "representatives_sha256": hashlib.sha256(json.dumps(representatives, separators=(",", ":")).encode()).hexdigest()}
    require(summary["gauge_phase_assignments"]*5**5 == raw, "translation gauge does not divide exactly")
    return representatives, summary


def verify(directory):
    from overnight30_fast_dual import verify_certificate
    manifest = json.loads((directory/"manifest.json").read_text())
    require(manifest["format"] == "ghz30-c5-general-phase-enumeration-v1", "unknown general phase manifest")
    source = Path(manifest["source_summary"])
    require(not source.is_absolute() and ".." not in source.parts, "source path must be repository relative")
    source_data = (ROOT/source).read_bytes()
    require(hashlib.sha256(source_data).hexdigest() == manifest["source_summary_sha256"], "source summary changed")
    index = manifest["source_index"]
    source = json.loads(source_data)
    require(type(index) is int and 0 <= index < len(source["results"]), "invalid source index")
    require(source["results"][index]["quotient"] == manifest["quotient"], "source quotient differs")
    representatives, coverage = enumerate_phases(manifest["quotient"])
    require(coverage == manifest["coverage"], "coverage enumeration differs")
    expected, found, bounds = set(representatives), set(), []
    for record in manifest["records"]:
        masks = tuple(record["phase_masks"])
        require(masks in expected and masks not in found, "unexpected or duplicate phase representative")
        edges = graph_edges(masks)
        require(graph6(edges) == record["graph6"], "phase graph differs")
        relative = Path(record["certificate"])
        require(not relative.is_absolute() and ".." not in relative.parts, "invalid certificate path")
        raw = gzip.decompress((directory/relative).read_bytes())
        require(hashlib.sha256(raw).hexdigest() == record["certificate_sha256"], "certificate hash differs")
        certificate = json.loads(raw)
        require(certificate["n"] == 30 and certificate["edges"] == [list(e) for e in edges], "certificate covers another graph")
        upper = verify_certificate(certificate)
        require(str(upper) == record["exact_theta_upper"], "recorded exact bound differs")
        bounds.append(upper)
        found.add(masks)
    complete = found == expected
    require(manifest["complete_exact_exclusion"] == complete, "incorrect completeness claim")
    return {"complete_exact_exclusion": complete, "covered_representatives": len(found),
            "total_representatives": len(expected), "coverage": coverage,
            "uniform_exact_upper": str(max(bounds)) if complete and bounds else None}


def search(args):
    from overnight30_fast_dual import verify_certificate
    start, deadline = time.monotonic(), time.monotonic()+args.seconds
    source = json.loads(args.source.read_text())
    quotient = source["results"][args.index]["quotient"]
    representatives, coverage = enumerate_phases(quotient)
    directory = args.output_dir
    directory.mkdir(parents=True, exist_ok=True)
    (directory/"certificates").mkdir(exist_ok=True)
    path = directory/"manifest.json"
    if path.exists():
        manifest = json.loads(path.read_text())
        require(manifest["quotient"] == quotient and manifest["coverage"] == coverage, "resume coverage differs")
        verify(directory)
    else:
        manifest = {"format": "ghz30-c5-general-phase-enumeration-v1", "quotient": quotient,
                    "source_index": args.index, "source_summary": str(args.source.resolve().relative_to(ROOT)),
                    "source_summary_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
                    "coverage": coverage, "records": [], "complete_exact_exclusion": False,
                    "scope": "all C5 phases of this fixed quotient; proper-subset tree translations and common units"}
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
            atomic_json(path, manifest)
    finally:
        worker.stop()
        manifest["complete_exact_exclusion"] = len(completed) == len(representatives)
        manifest["last_discovery_seconds"] = time.monotonic()-start
        manifest["unresolved_attempts"] = pending
        atomic_json(path, manifest)
    report = verify(directory)
    report["discovery_seconds"] = manifest["last_discovery_seconds"]
    print(json.dumps(report), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("search")
    run.add_argument("--source", type=Path, default=SOURCE)
    run.add_argument("--index", type=int, required=True)
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

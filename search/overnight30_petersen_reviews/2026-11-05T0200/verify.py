"""Replay the bounded Petersen-template pilot archives using exact arithmetic."""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use -O.')

import argparse
from collections import Counter
import hashlib
import io
import itertools
import json
from pathlib import Path
import sys
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
from overnight30_screen import base_validation, graph_id, verify_certificate


def independent_geometry():
    vertices = list(itertools.combinations(range(5), 2))
    neighbors = [{j for j, b in enumerate(vertices) if i == j or set(a).isdisjoint(b)}
                 for i, a in enumerate(vertices)]
    # Enumerate right subsets, independently of discovery's left-subset loop.
    rectangles = set()
    for size in range(1, 11):
        for choice in itertools.combinations(range(10), size):
            right = set(choice)
            left = {i for i in range(10) if not neighbors[i] & right}
            if left and {j for j in range(10) if all(j not in neighbors[i] for i in left)} == right:
                rectangles.add((tuple(sorted(left)), tuple(sorted(right))))
    # Discovery orders left sets by binary mask.
    rectangles = sorted(rectangles, key=lambda pair: sum(1 << i for i in pair[0]))
    assert len(rectangles) == 110
    assert Counter(tuple(map(len, pair)) for pair in rectangles) == {
        (1, 6): 10, (6, 1): 10, (2, 4): 15, (4, 2): 15, (2, 3): 30, (3, 2): 30}
    index = {pair: i for i, pair in enumerate(rectangles)}
    permutations = []
    for sigma in itertools.permutations(range(5)):
        p = [vertices.index(tuple(sorted(sigma[i] for i in vertex))) for vertex in vertices]
        for swap in (False, True):
            images = []
            for left, right in rectangles:
                pair = (tuple(sorted(p[i] for i in left)), tuple(sorted(p[i] for i in right)))
                images.append(index[pair[::-1] if swap else pair])
            assert sorted(images) == list(range(110))
            permutations.append(images)
    assert len(set(map(tuple, permutations))) == 240
    return neighbors, rectangles, permutations


def verify_payload(files):
    summary = json.loads(files['summary.json'])
    provenance = json.loads(files['provenance.json'])
    assert hashlib.sha256(files['discovery_helper.py']).hexdigest() == provenance['helper_sha256']
    neighbors, rectangles, permutations = independent_geometry()
    rows = [json.loads(line) for line in files['records.jsonl'].splitlines()]
    parts = [list(range(10*k, 10*k+10)) for k in range(3)]
    seen, canonical, methods = set(), set(), Counter()
    for row in rows:
        selected = row['selected']
        assert selected == sorted(set(selected)) and len(selected) == 10
        assert all(type(i) is int and 0 <= i < 110 for i in selected)
        edges = {(i, 10+j) for i in range(10) for j in neighbors[i]}
        for c, i in enumerate(selected, 20):
            edges.update((a, c) for a in rectangles[i][0])
            edges.update((10+b, c) for b in rectangles[i][1])
        graph = {'n': 30, 'parts': parts, 'edges': [list(edge) for edge in sorted(edges)]}
        assert graph == row['graph'] and graph_id(graph) == row['graph_id']
        adj, _ = base_validation(graph)
        assert all(1 <= len(adj[v] & set(parts[2])) <= 6 for v in range(20))
        for part in parts:
            for u, v in itertools.combinations(part, 2):
                assert any(adj[a] & adj[v] for a in adj[u])
        assert sum(min(map(len, rectangles[i])) == 1 for i in selected) <= 1
        if summary.get('perfect_matching_constraints'):
            assert len(row['matchings']) == 2
            for side, pairs in enumerate(row['matchings']):
                assert sorted(v for v, i in pairs) == list(range(10))
                assert sorted(i for v, i in pairs) == selected
                assert all(v in rectangles[i][side] for v, i in pairs)
        orbit = {tuple(sorted(p[i] for i in selected)) for p in permutations}
        assert not seen & orbit
        seen.update(orbit)
        canonical.add(min(orbit))
        if row['certificate'] is not None:
            assert row['certificate']['graph'] == graph
            assert row['status'] == 'excluded' and verify_certificate(row['certificate'])
            methods[row['certificate']['method']] += 1
    assert len(rows) == summary['representatives'] == len(canonical)
    assert len(seen) == summary['selected_subsets_blocked']
    assert dict(methods) == summary['exclusion_methods']
    assert len(rows)-sum(methods.values()) == summary['survivors']
    assert summary['numerical_solves'] == 0
    return {'records_verified': len(rows), 'exact_certificates_verified': sum(methods.values()),
            'explicit_action_orbits_verified': len(canonical),
            'blocked_selected_subsets': len(seen), 'numerical_solves': 0}, canonical


def main(data_repo):
    summary = json.loads((HERE/'summary.json').read_text())
    total, overlap, seen = 0, 0, set()
    for pilot in summary['pilots']:
        path = Path(pilot['archive_path'])
        assert not path.is_absolute() and '..' not in path.parts
        raw = (data_repo/path).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == pilot['archive_sha256']
        with tarfile.open(fileobj=io.BytesIO(raw), mode='r:gz') as archive:
            files = {member.name: archive.extractfile(member).read()
                     for member in archive.getmembers() if member.isfile()}
        report, canonical = verify_payload(files)
        assert report == pilot['verification']
        total += report['records_verified']
        overlap += len(seen & canonical)
        seen.update(canonical)
    assert total == summary['combined_generated_representatives']
    assert overlap == summary['cross_pilot_explicit_action_overlap']
    assert len(seen) == summary['combined_explicit_action_orbits']
    print(f'Verified {total} pilot records and exact individual exclusions; '
          f'{len(seen)} distinct orbits under the declared action, no family exclusion.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-repo', type=Path, default=Path('/workspace/ghz30-overnight'))
    args = parser.parse_args()
    main(args.data_repo)

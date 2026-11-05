"""Standard-library replay of exact move-audit certificates and provenance."""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use -O.')

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
from overnight30_screen import base_validation, components, verify_certificate
from analyze import analyze, COMMIT


def graph(code):
    assert len(code) == 74 and ord(code[0]) == 93
    assert all(63 <= ord(char) <= 126 for char in code)
    edges = []
    for v in range(1, 30):
        for u in range(v):
            bit = v*(v-1)//2+u
            if (ord(code[1+bit//6])-63) & (1 << (5-bit % 6)):
                edges.append([u, v])
    return {'n': 30, 'parts': [list(range(10*k, 10*k+10)) for k in range(3)], 'edges': sorted(edges)}


def edge(a, b):
    return tuple(sorted((a, b)))


def verify_provenance(record, parent):
    value = graph(record['g6'])
    adjacency, _ = base_validation(value)
    assert len(components(adjacency)) == 1
    assert all(len(neighbors) <= 10 for neighbors in adjacency)
    assert all(adjacency[v] & set(value['parts'][k])
               for v in range(30) for k in range(3) if k != v//10)
    provenance = record['provenance']
    assert provenance['parent_g6'] == parent['g6']
    original = graph(parent['g6'])
    before = set(map(tuple, original['edges']))
    after = set(map(tuple, value['edges']))
    repairs = [provenance] if record['mode'] == 'single' else provenance['repairs']
    assert len(repairs) == (1 if record['mode'] == 'single' else 2)
    used = set()
    protected = set()
    for repair in repairs:
        block = repair['certificate_block']
        assert block in parent['certificate']['witness']['blocks']
        assert not used.intersection(block['vertices'])
        used.update(block['vertices'])
        u, v = repair['sharing_pair']; b = repair['shared_neighbor']
        assert u != v and u in block['vertices'] and v in block['vertices']
        assert b in original['parts'][block['other_context']]
        assert repair['added_edge'] == [u, b]
        assert edge(u, b) not in before and edge(v, b) in before
        protected.update((edge(u, b), edge(v, b)))
    current = before.copy()
    for repair in repairs:
        u, v = repair['sharing_pair']; b = repair['shared_neighbor']
        common = [c for c in range(30) if edge(u, c) in current and edge(b, c) in current]
        removed = set(map(lambda item: edge(*item), repair['removed_conflicts']))
        assert len(removed) == len(common) and removed <= current and not removed.intersection(protected)
        assert all(sum(edge(endpoint, c) in removed for endpoint in (u, b)) == 1 for c in common)
        current.difference_update(removed)
        current.add(edge(u, b))
    assert current <= after and protected <= after
    assert provenance['edge_difference'] == len(before ^ after)
    proof = record['certificate']
    if proof:
        assert proof['graph'] == value and record['screen_status'] == 'excluded'
        assert verify_certificate(proof)
    else:
        assert record['screen_status'] == 'not_excluded'


def main(repo):
    expected = json.loads((HERE/'natural_performance.json').read_text())
    assert analyze(repo) == expected
    summary = json.loads((HERE/'natural_screen.json').read_text())
    original = {}
    for path, digest in summary['inputs_sha256'].items():
        raw = subprocess.check_output(['git', 'show', COMMIT+':'+path], cwd=repo)
        assert hashlib.sha256(raw).hexdigest() == digest
        for index, line in enumerate(gzip.decompress(raw).splitlines()):
            row = json.loads(line)
            if 'partition_repair' in row:
                original[row['g6']] = (path, index)
    assert len(original) == summary['distinct_logged_repair_labels'] == 154
    packed = (HERE/summary['certificate_archive']).read_bytes()
    assert hashlib.sha256(packed).hexdigest() == summary['certificate_archive_sha256']
    seen, counts = set(), Counter()
    for line in gzip.decompress(packed).splitlines():
        row = json.loads(line)
        assert row['g6'] not in seen and original[row['g6']] == (row['source'], row['line'])
        seen.add(row['g6'])
        proof = row['certificate']
        assert proof['graph'] == graph(row['g6']) and verify_certificate(proof)
        counts[proof['method']] += 1
        if proof['method'] == 'disjoint_partition':
            assert sorted(len(block['vertices']) for block in proof['witness']['blocks']) == [1, 4, 5]
            assert proof['witness']['budget'] == '1411/1440'
    counts['not_excluded'] = len(original)-len(seen)
    assert dict(counts) == summary['counts']
    assert Fraction(15, 32)+Fraction(2, 5)+Fraction(1, 9) == Fraction(1411, 1440) < 1

    reference = json.loads((HERE/'raw_data.json').read_text())
    raw = (repo/reference['archive_path']).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == reference['archive_sha256']
    with tarfile.open(fileobj=io.BytesIO(raw), mode='r:gz') as archive:
        report = json.load(archive.extractfile('paired_control.json'))
        records = [json.loads(line) for line in gzip.decompress(archive.extractfile('paired_control_records.jsonl.gz').read()).splitlines()]
    assert report == json.loads((HERE/'paired_control.json').read_text())
    for path, digest in report['inputs_sha256'].items():
        revision, location = (report['code_commit'], ROOT) if path.startswith('search/overnight30_') and '/chunks/' not in path else (report['data_commit'], repo)
        source = subprocess.check_output(['git', 'show', revision+':'+path], cwd=location)
        assert hashlib.sha256(source).hexdigest() == digest
    parent = json.loads(subprocess.check_output(['git', 'show', report['data_commit']+':'+report['parent']], cwd=repo))
    assert verify_certificate(parent['certificate']) and parent['certificate']['graph'] == graph(parent['g6'])
    for mode in ('single', 'paired'):
        selected = [row for row in records if row['mode'] == mode]
        assert len(selected) == report['matched_seeds_per_mode'] == 128
        assert sorted(row['seed'] for row in selected) == list(range(202611050200, 202611050328))
        observed = Counter(valid=128)
        for row in selected:
            verify_provenance(row, parent)
            observed[row['certificate']['method'] if row['certificate'] else 'not_excluded'] += 1
        assert dict(observed) == report['statistics'][mode]['counts']
        assert len({row['g6'] for row in selected}) == report['statistics'][mode]['distinct_labels']
    print('Verified natural counters, 81 exact repair exclusions, 256 matched proposals and all protected-edge provenance.')
    print('No numerical solves; screen survival is not quantum feasibility.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-repo', type=Path, default=Path('/workspace/ghz30-overnight'))
    main(parser.parse_args().data_repo)

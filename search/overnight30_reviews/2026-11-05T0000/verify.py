"""Replay exact certificates and saved isomorphisms without numerical packages."""
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

DIRECTORY = Path(__file__).resolve().parent
sys.path.insert(0, str(DIRECTORY.parents[1]))
from overnight30_screen import graph_id, normalize_graph, verify_certificate
from overnight30_fast_dual import verify_certificate as verify_fast
from overnight30_quotient import verify_certificate as verify_quotient


def verify_any(certificate):
    if certificate['format'] == 'overnight30-factor-dual-v1':
        return verify_fast(certificate)
    if certificate['format'] == 'overnight30-worker-obstruction-v1':
        return verify_quotient(certificate)
    return verify_certificate(certificate)


def certificate_graph(certificate):
    return certificate.get('graph', {'n': certificate.get('n'),
                                     'edges': certificate.get('edges')})


def rows(name):
    return list(map(json.loads, (DIRECTORY/name).read_text().splitlines()))


def edges30(code):
    code = code.removeprefix('>>graph6<<').strip()
    values = [ord(c)-63 for c in code]
    assert len(values) == 74 and values[0] == 30 and all(0 <= v <= 63 for v in values)
    bits = [bool(v & (1 << bit)) for v in values[1:] for bit in range(5, -1, -1)]
    result, offset = set(), 0
    for v in range(1, 30):
        for u in range(v):
            if bits[offset]:
                result.add((u, v))
            offset += 1
    assert not any(bits[offset:])
    return result


hashes = json.loads((DIRECTORY/'artifact_sha256.json').read_text())
for name, digest in hashes.items():
    assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest() == digest, name

certificates = {}
for path in DIRECTORY.glob('*certificate.json'):
    certificate = json.loads(path.read_text())
    verify_any(certificate)
    certificates[path.name] = certificate

for row in rows('strongest_candidates.jsonl')+rows('failure_rechecks.jsonl'):
    graph = normalize_graph({'n': 30, 'edges': sorted(edges30(row['g6']))})
    assert graph_id(graph) == row['graph_id']
    for field in ('certificate', 'dimension19_certificate'):
        if row.get(field):
            assert normalize_graph(certificate_graph(certificates[row[field]])) == graph
    if row.get('reused_certificate'):
        certificate = certificates[row['certificate']]
        assert str(verify_fast(certificate)) == row['reused_certificate']['bound']

representatives = {}
folders = [(DIRECTORY, 'selected_representatives.jsonl')]
folders.extend((DIRECTORY.parent/date, 'strongest_candidates.jsonl')
               for date in ('2026-11-04T2200', '2026-11-04T2300'))
for folder, name in folders:
    for row in map(json.loads, (folder/name).read_text().splitlines()):
        edges = edges30(row['g6'])
        assert graph_id(normalize_graph({'n': 30, 'edges': sorted(edges)})) == row['graph_id']
        representatives[(folder.name, row['graph_id'])] = edges

equivalences = rows('isomorphic_candidates.jsonl')
for row in equivalences:
    mapping = {int(k): v for k, v in row['vertex_mapping'].items()}
    assert set(mapping) == set(range(30)) and set(mapping.values()) == set(range(30))
    transformed = {tuple(sorted((mapping[u], mapping[v]))) for u, v in edges30(row['g6'])}
    assert transformed == representatives[(row['representative_review'], row['representative'])]

references = json.loads((DIRECTORY/'referenced_certificates.json').read_text())
for reference in references:
    raw = (DIRECTORY/reference['path']).read_bytes()
    assert hashlib.sha256(raw).hexdigest() == reference['sha256']
    certificate = json.loads(raw)
    verify_any(certificate)
    assert graph_id(normalize_graph(certificate_graph(certificate))) == reference['graph_id']

if '--inputs' in sys.argv:
    data = Path('/workspace/ghz30-overnight')
    commit = json.loads((DIRECTORY/'summary.json').read_text())['published_commit']
    for row in rows('inputs.jsonl'):
        path = data/row['path']
        raw = path.read_bytes() if path.exists() else b''
        if hashlib.sha256(raw).hexdigest() != row['sha256']:
            raw = subprocess.check_output(['git', 'show', f'{commit}:{row["path"]}'], cwd=data)
        assert hashlib.sha256(raw).hexdigest() == row['sha256']
        assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest() == row['git_blob']
    for row in rows('strongest_candidates.jsonl'):
        if not row.get('reused_certificate'):
            continue
        ref = row['reused_certificate']
        raw = subprocess.check_output(['git', 'show', f'{ref["commit"]}:{ref["path"]}'], cwd=data)
        decoded = gzip.decompress(raw)
        assert hashlib.sha256(decoded).hexdigest() == ref['json_sha256']
        assert json.loads(decoded) == certificates[row['certificate']]
    print('Published input hashes and reused certificate source verified.')

print(f'Verified {len(certificates)} exact certificates, {len(equivalences)} isomorphism mappings, '
      f'{len(references)} referenced proofs, and artifact hashes.')

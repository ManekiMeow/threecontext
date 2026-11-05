"""Replay exact certificates and saved isomorphisms without numerical packages."""
import gzip
import collections
import math
import hashlib
import json
from pathlib import Path
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use python -O')

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

supplement = json.loads((DIRECTORY/'active_best_review.json').read_text()) if (DIRECTORY/'active_best_review.json').exists() else None
for row in rows('strongest_candidates.jsonl')+rows('failure_rechecks.jsonl')+([supplement] if supplement else []):
    assert row['status'] == 'excluded' and row.get('certificate')
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
               for date in ('2026-11-04T2200', '2026-11-04T2300', '2026-11-05T0000', '2026-11-05T0100'))
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

history = json.loads((DIRECTORY/'snapshot_asymmetric_job98_history.json').read_text())
active = json.loads((DIRECTORY/'snapshot_asymmetric_job98_best.json').read_text())
assert active['theta_numeric'] == history['best_theta']
assert history['historical_best_theta'] >= history['best_theta']
for record in history['retired_best_records']:
    proof = record['certificate']
    verify_any(proof)
    assert normalize_graph(certificate_graph(proof)) == normalize_graph({'n':30,'edges':sorted(edges30(record['g6']))})
    assert record['theta_numeric'] <= history['best_excluded_theta']
audit = json.loads((DIRECTORY/'automorphism_audit.json').read_text())
assert audit['best']['g6'] == active['g6']
assert audit['best_excluded']['g6'] in {row['g6'] for row in history['retired_best_records']}
if supplement:
    assert supplement['g6'] == active['g6'] and supplement['theta_numeric'] == active['theta_numeric']

from audit_closure_cuts import verify as verify_closure_cuts
closure_report = verify_closure_cuts('--inputs' in sys.argv)
assert closure_report == json.loads((DIRECTORY/'closure_cut_audit.json').read_text())

summary = json.loads((DIRECTORY/'summary.json').read_text())
selected = rows('selected_representatives.jsonl')
assert len(selected) == summary['selected_new_isomorphism_representatives'] == 124
assert len(equivalences) == summary['selected_isomorphism_equivalences'] == 3
assert len(selected)+len(equivalences) == summary['selected_distinct_labeled_supports'] == 127
assert len(rows('strongest_candidates.jsonl')) == summary['exactly_reviewed_new_representatives'] == 19
assert len(certificates) == 21 and len(rows('failure_rechecks.jsonl')) == 2
assert active['g6'] == rows('strongest_candidates.jsonl')[0]['g6']

if '--inputs' in sys.argv:
    data = Path('/workspace/ghz30-overnight')
    commit = json.loads((DIRECTORY/'summary.json').read_text())['published_commit']
    sources = {row['source'] for row in selected+equivalences}
    source_bindings = {}
    labels = {False:set(), True:set()}; scored = {False:set(), True:set()}
    records = collections.Counter(); score_records = collections.Counter()
    for row in rows('inputs.jsonl'):
        path = data/row['path']
        raw = path.read_bytes() if path.exists() else b''
        if hashlib.sha256(raw).hexdigest() != row['sha256']:
            raw = subprocess.check_output(['git', 'show', f'{commit}:{row["path"]}'], cwd=data)
        assert hashlib.sha256(raw).hexdigest() == row['sha256']
        assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest() == row['git_blob']
        path_parts = row['path'].split('/')
        main_pool = row['path'].endswith('.jsonl.gz') and (row['path'].endswith('/models.jsonl.gz') or '/asymmetric_seed' in row['path'])
        if main_pool:
            job = int(next(part.split('_')[0] for part in path_parts if part[:5].isdigit() and '_' in part))
            fresh = job in summary['new_completed_jobs']
            for line in gzip.decompress(raw).splitlines():
                if not line.strip(): continue
                item=json.loads(line); records[fresh]+=1
                code=item.get('g6',item.get('graph6'))
                if code: code=code.removeprefix('>>graph6<<').strip(); labels[fresh].add(code)
                value=item.get('theta_numeric',item.get('theta'))
                if code and isinstance(value,(int,float)) and math.isfinite(value):
                    score_records[fresh]+=1; scored[fresh].add(code)
                    if row['path'] in sources:
                        source_bindings[(row['path'],item.get('index',item.get('iteration')),code)] = value
    assert len(labels[False]) == summary['baseline_distinct_labeled_supports']
    assert len(scored[False]) == summary['baseline_distinct_scored_labels']
    assert len(labels[True]) == summary['new_chunk_distinct_labeled_supports']
    assert len(scored[True]) == summary['new_chunk_distinct_scored_labels']
    assert len(labels[True]-labels[False]) == summary['new_labeled_supports_vs_baseline']
    assert len(scored[True]-scored[False]) == summary['new_scored_labels_vs_baseline']
    assert records[True] == summary['new_raw_records'] and score_records[True] == summary['new_scored_records']
    for row in selected+equivalences:
        assert source_bindings[row['source'],row['index'],row['g6']] == row['theta_numeric']
    for row in rows('strongest_candidates.jsonl'):
        if not row.get('reused_certificate'):
            continue
        ref = row['reused_certificate']
        raw = subprocess.check_output(['git', 'show', f'{ref["commit"]}:{ref["path"]}'], cwd=data)
        decoded = gzip.decompress(raw)
        assert hashlib.sha256(decoded).hexdigest() == ref['json_sha256']
        assert json.loads(decoded) == certificates[row['certificate']]
    raw = subprocess.check_output(['git','show',f"{commit}:{history['source']}"],cwd=data)
    checkpoint = json.loads(raw)
    assert all(checkpoint.get(key)==value for key,value in history.items() if key not in ('published_commit','source'))
    print('Published input hashes, full label/score coverage, selected-record bindings, history and reused proofs verified.')

print(f'Verified {len(certificates)} exact certificates, {len(equivalences)} isomorphism mappings, '
      f'{len(references)} referenced proofs, 281 natural closure cuts, and artifact hashes.')

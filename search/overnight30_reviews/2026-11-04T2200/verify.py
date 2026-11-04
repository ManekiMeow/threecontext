"""Verify this review's new exact certificates without SAT/SDP dependencies."""
import hashlib
import json
from pathlib import Path
import sys

directory = Path(__file__).resolve().parent
sys.path.insert(0, str(directory.parents[1]))
from overnight30_screen import verify_certificate

summary = json.loads((directory / 'summary.json').read_text())
snapshot = (directory / 'snapshot_00023_c5_double_best.json').read_bytes()
assert hashlib.sha256(snapshot).hexdigest() == summary['active_snapshot_sha256']
paths = sorted(directory.glob('*certificate.json'))
for path in paths:
    verify_certificate(json.loads(path.read_text()))
print(f'Verified {len(paths)} exact certificates and the immutable snapshot hash.')

from overnight30_fast_dual import verify_certificate as verify_factor_certificate
factor_bound = verify_factor_certificate(json.loads((directory/'fast_factor_c5_double.json').read_text()))
print(f'Verified the alternative integer-factor certificate: theta <= {factor_bound}.')

if '--isomorphism' in sys.argv:
    import networkx as nx
    representatives = {row['graph_id']: nx.from_graph6_bytes(row['g6'].encode())
                       for row in map(json.loads, (directory/'strongest_candidates.jsonl').read_text().splitlines())}
    duplicates = list(map(json.loads, (directory/'isomorphic_candidates.jsonl').read_text().splitlines()))
    for row in duplicates:
        graph = nx.from_graph6_bytes(row['g6'].encode())
        assert nx.is_isomorphic(graph, representatives[row['representative']])
    print(f'Verified {len(duplicates)} exact uncolored graph-isomorphism claims.')

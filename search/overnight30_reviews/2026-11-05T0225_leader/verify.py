"""Exact supplement replay for the population-crossover leader."""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use -O')
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from overnight30_screen import normalize_graph,graph_id,verify_certificate


def edges30(code):
    values=[ord(c)-63 for c in code]
    assert len(values)==74 and values[0]==30 and all(0<=v<=63 for v in values)
    bits=[bool(value & (1<<bit)) for value in values[1:] for bit in range(5,-1,-1)]
    edges=[];offset=0
    for v in range(1,30):
        for u in range(v):
            if bits[offset]:edges.append([u,v])
            offset+=1
    assert not any(bits[offset:])
    return sorted(edges)


for name,digest in json.loads((HERE/'artifact_sha256.json').read_text()).items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
raw=(HERE/'live_best_snapshot.json').read_bytes();snap=json.loads(raw)
meta=json.loads((HERE/'snapshot_source.json').read_text())
assert hashlib.sha256(raw).hexdigest()==meta['source_sha256']==meta['published_source_sha256']
assert len(raw)==meta['source_bytes'] and meta['same_as_published_best_bytes'] is True
graph=normalize_graph({'n':30,'edges':edges30(snap['g6'])})
assert graph==normalize_graph({'parts':snap['parts'],'edges':snap['edges']})
review=json.loads((HERE/'candidate_review.json').read_text())
assert review['g6']==snap['g6'] and review['theta_numeric']==snap['theta_numeric']
assert graph_id(graph)==review['graph_id']=='097316a813a1e997f5f6'
proof=json.loads((HERE/review['certificate']).read_text())
assert proof['graph']==graph and verify_certificate(proof)
assert proof['method']=='theta_dual' and proof['witness']['bound']=='29213/10000'
provenance=json.loads((HERE/'published_provenance.json').read_text())
assert len(provenance)==1
source=provenance[0];record=source['record']
assert record['g6']==snap['g6'] and record['theta']==snap['theta_numeric']
assert record['iteration']==snap['iteration']==287180
assert record['move']=='population_crossover' and 'partition_repair' not in record
assert source['job']==101 and source['published_commit']==meta['published_commit']
auto=json.loads((HERE/'automorphism_audit.json').read_text())
assert auto['g6']==snap['g6'] and auto['full_automorphism_group_order']==auto['context_preserving_automorphism_group_order']==1
if '--inputs' in sys.argv:
    data='/workspace/ghz30-overnight'
    published=subprocess.check_output(['git','show',meta['published_commit']+':'+meta['published_source']],cwd=data)
    assert published==raw
    archive=subprocess.check_output(['git','show',source['published_commit']+':'+source['source']],cwd=data)
    assert hashlib.sha256(archive).hexdigest()==source['source_sha256']
    assert json.loads(gzip.decompress(archive).splitlines()[source['source_line']])==record
    print('Published best snapshot bytes and original population-crossover scoring record verified.')
print('Verified exact theta upper29213/10000, graph binding and immutable artifact hashes.')

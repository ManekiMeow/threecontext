"""Cached exact probability-template filter for the asymmetric worker only.

The immutable48-template catalog is hashed and every source proof is checked
once per process. Target discovery uses exact arithmetic, without an LP/SDP.
Every emitted target proof is independently verified in the existing vertex
probability format. The outer global certificate is self-contained: replay
needs no template catalog, only an exact inner/outer graph binding.
"""
if not __debug__:
    raise RuntimeError('Exact probability transfer requires assertions; do not use -O.')

from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import runpy

from overnight30_screen import normalize_graph,certificate
from event30_vertex_probability_lp import geometry

ROOT=Path(__file__).resolve().parent.parent
CATALOG=ROOT/'search/overnight30_probability_transfer/2026-11-05T0200'


@lru_cache(maxsize=1)
def verified_catalog():
    """Immutable prepared templates; bounded process-local cache, no target cache."""
    hashes=json.loads((CATALOG/'artifact_sha256.json').read_text())
    raw_inputs=(CATALOG/'inputs.json').read_bytes()
    if hashlib.sha256(raw_inputs).hexdigest()!=hashes['inputs.json']:
        raise ValueError('probability-transfer catalog metadata hash mismatch')
    inputs=json.loads(raw_inputs)
    relative=Path(inputs['template_source'])
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('invalid probability-template archive path')
    raw=(ROOT/relative).read_bytes()
    digest=hashlib.sha256(raw).hexdigest()
    if digest!=inputs['template_source_sha256']:
        raise ValueError('probability-template proof archive hash mismatch')
    rows=list(map(json.loads,gzip.decompress(raw).splitlines()))
    proofs={row['name']:row['certificate'] for row in rows}
    if len(proofs)!=len(rows):
        raise ValueError('duplicate probability-template source name')
    names=inputs['template_names']
    if len(names)!=48 or len(set(names))!=48:
        raise ValueError('expected48 distinct frozen probability templates')
    algorithm=runpy.run_path(str(CATALOG/'transfer.py'))
    prepared=tuple((name,algorithm['prepare_template'](proofs[name])) for name in names)
    return digest,prepared,algorithm['transfer']


def find_transfer(record):
    """Return a verified global wrapper, or None; never call a numerical solver."""
    graph=normalize_graph(record)
    if graph['n']!=30 or list(map(len,graph['parts']))!=[10,10,10]:
        return None
    base=geometry(graph)
    digest,templates,transfer=verified_catalog()
    for name,template in templates:
        proof,details=transfer(template,graph,base)
        if proof is not None:
            # transfer() independently verifies the full target identity before
            # returning it. The outer graph is exactly that verified graph.
            assert proof['graph']==graph
            result=certificate(graph,'vertex_probability',certificate=proof)
            result['discovery']={'source_template':name,'source_archive_sha256':digest,
                                 'added_upper_rows':details['added_upper_rows']}
            return result
    return None

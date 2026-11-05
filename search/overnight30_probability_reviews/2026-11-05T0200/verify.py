"""Standard-library exact verification of the bounded vertex probability review."""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use -O.')
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'search'))
from event30_vertex_probability_lp import verify_certificate, self_test


def graph_from_g6(code):
    data=code.encode('ascii')
    assert len(data)==74 and data[0]==93 and all(63<=x<=126 for x in data)
    edges=[]
    for v in range(1,30):
        for u in range(v):
            bit=v*(v-1)//2+u
            if (data[1+bit//6]-63) & (1<<(5-bit%6)):
                edges.append([u,v])
    return {'n':30,'parts':[list(range(10*k,10*k+10)) for k in range(3)],'edges':sorted(edges)}


def read_rows(path):
    with gzip.open(path,'rt') as stream:
        return [json.loads(line) for line in stream]


def main():
    hashes=json.loads((HERE/'artifact_sha256.json').read_text())
    for name,digest in hashes.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
    summary=json.loads((HERE/'summary.json').read_text())
    inputs=[json.loads(line) for line in (HERE/'inputs.jsonl').read_text().splitlines()]
    source=ROOT/summary['selected_source']
    assert source.read_bytes()==(HERE/'inputs.jsonl').read_bytes()
    assert hashlib.sha256(source.read_bytes()).hexdigest()==summary['selected_source_sha256']
    assert len(inputs)==summary['selected_input_count']
    results=read_rows(HERE/summary['results_archive'])
    packed=read_rows(HERE/summary['proof_archive'])
    proofs={row['name']:row['certificate'] for row in packed}
    assert len(proofs)==len(packed)==summary['new_exact_exclusions']
    used,kinds=set(),Counter()
    for index,row in enumerate(results):
        assert index==row['selection_index']
        assert all(row[k]==v for k,v in inputs[index].items())
        assert not (row['already_structurally_excluded'] and row['new_exact_exclusion'])
        if row['new_exact_exclusion']:
            proof=proofs[row['certificate']]
            assert proof['graph']==graph_from_g6(row['g6'])
            verify_certificate(proof)
            assert row['contradiction_margin']==proof['contradiction_margin']
            assert row['strict_multiplier_mass']==proof['strict_multiplier_mass']
            used.add(row['certificate'])
            kinds.update(term['constraint']['kind'] for term in proof['inequalities'])
    assert used==set(proofs) and dict(kinds)==summary['constraint_kinds_used']
    assert len(results)==summary['reviewed']
    assert sum(row['already_structurally_excluded'] for row in results)==summary['existing_structural_exclusions']
    assert sum(not row['already_structurally_excluded'] for row in results)==summary['lp_attempted']
    assert len(inputs)-len(results)==summary['unreached']
    audit_raw=(ROOT/summary['prior_selected_exclusion_source']).read_bytes()
    assert hashlib.sha256(audit_raw).hexdigest()==summary['prior_selected_exclusion_source_sha256']
    earlier={row['graph_id'] for row in map(json.loads,audit_raw.splitlines()) if row['status']=='excluded'}
    certified={row['graph_id'] for row in results if row['new_exact_exclusion']}
    assert len(certified)==summary['additional_scalar_certificates']
    assert len(certified & earlier)==summary['overlap_with_prior_selected_exclusions']
    assert len(certified-earlier)==summary['newly_excluded_selected_representatives']
    assert len(certified | earlier)==summary['combined_selected_exclusions']
    assert len(inputs)-len(certified | earlier)==summary['unresolved_selected_representatives']
    assert self_test()==summary['controls']
    print(f"Verified {len(proofs)} exact per-vertex probability exclusions from {summary['lp_attempted']} LP attempts.")


if __name__=='__main__':
    main()

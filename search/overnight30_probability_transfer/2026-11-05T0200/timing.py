"""Exact replay timing, separated from immutable transfer discovery evidence."""
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import time
import transfer as transfer_module

HERE=Path(__file__).resolve().parent
ROOT=transfer_module.ROOT
MAIN=ROOT/'search/overnight30_probability_reviews/2026-11-05T0200'
graph_from_g6=runpy.run_path(str(MAIN/'verify.py'))['graph_from_g6']
from overnight30_screen import fast_screen


def main():
    if (HERE/'timing.json').exists():
        raise RuntimeError('Timing evidence already exists; preserve it.')
    inputs=json.loads((HERE/'inputs.json').read_text())
    results=list(map(json.loads,gzip.open(HERE/'results.jsonl.gz','rt')))
    catalog_started=time.perf_counter()
    raw=(ROOT/inputs['template_source']).read_bytes()
    assert hashlib.sha256(raw).hexdigest()==inputs['template_source_sha256']
    source={row['name']:row['certificate'] for row in map(json.loads,gzip.decompress(raw).splitlines())}
    prepared=[(name,transfer_module.prepare_template(source[name])) for name in inputs['template_names']]
    catalog_seconds=time.perf_counter()-catalog_started
    exact_verify=transfer_module.verify_certificate
    verify_seconds=verify_count=0
    def timed_verify(proof):
        nonlocal verify_seconds,verify_count
        start=time.perf_counter()
        result=exact_verify(proof)
        verify_seconds+=time.perf_counter()-start
        verify_count+=1
        return result
    transfer_module.verify_certificate=timed_verify
    totals={'fast_screen_seconds':0.0,'target_geometry_seconds':0.0,'transfer_inclusive_seconds':0.0}
    attempts=0
    start=time.perf_counter()
    for expected in results:
        graph=graph_from_g6(expected['g6'])
        begin=time.perf_counter(); screen=fast_screen(graph)
        totals['fast_screen_seconds']+=time.perf_counter()-begin
        assert (screen['status']=='excluded')==expected['existing_exclusion']
        if expected['existing_exclusion']:
            continue
        begin=time.perf_counter(); base=transfer_module.geometry(graph)
        totals['target_geometry_seconds']+=time.perf_counter()-begin
        found=False
        for name,template in prepared:
            begin=time.perf_counter(); proof,details=transfer_module.transfer(template,graph,base)
            totals['transfer_inclusive_seconds']+=time.perf_counter()-begin
            attempts+=1
            if proof is not None:
                assert expected['transferred'] and name==expected['source_template']
                assert details==expected['details']
                found=True
                break
        assert found==expected['transferred']
    elapsed=time.perf_counter()-start
    transfer_module.verify_certificate=exact_verify
    result={'scope':'Replay timing only; no additional targets or new certificates.',
            'processed_targets':len(results),'existing_fast_exclusions':sum(row['existing_exclusion'] for row in results),
            'template_attempts':attempts,'final_exact_proofs_verified':verify_count,
            'cold_catalog_hash_and_exact_validation_seconds':catalog_seconds,
            **totals,'final_exact_certificate_replay_seconds':verify_seconds,
            'transfer_without_final_replay_seconds':totals['transfer_inclusive_seconds']-verify_seconds,
            'full_replay_seconds':elapsed,
            'cache_scope':'Prepared source templates may be cached after exact validation; every emitted target certificate is still independently verified.'}
    (HERE/'timing.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()

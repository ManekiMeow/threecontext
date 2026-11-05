"""Nonnumerical API, exact-wrapper and checkpoint/repair-parent controls."""
from copy import deepcopy
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import random
import sys
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'search'))
import event30_probability_transfer as api
from event30_vertex_probability_lp import verify_certificate as verify_inner
from overnight30_screen import verify_certificate,normalize_graph,certificate


def reject(function):
    try:
        function()
    except (AssertionError,ValueError):
        return
    raise AssertionError('invalid integration control accepted')


def main():
    assert not any(name in sys.modules for name in ('numpy','scipy','sympy','cvxpy'))
    saved=list(map(json.loads,gzip.open(HERE/'certificates.jsonl.gz','rt')))
    wrappers=[]
    for row in saved:
        proof=api.find_transfer(row['certificate']['graph'])
        assert proof is not None and proof['method']=='vertex_probability'
        assert proof['witness']['certificate']==row['certificate']
        verify_certificate(proof)
        wrappers.append(proof)
    assert len(wrappers)==52 and api.verified_catalog.cache_info().misses==1
    assert not any(name in sys.modules for name in ('numpy','scipy','sympy','cvxpy'))
    cold_loader=api.verified_catalog
    def unavailable():
        raise AssertionError('self-contained replay touched discovery catalog')
    api.verified_catalog=unavailable
    try:
        verify_certificate(wrappers[0])
    finally:
        api.verified_catalog=cold_loader
    bad=deepcopy(wrappers[0]);bad['graph']['edges'].pop()
    reject(lambda:verify_certificate(bad))
    bad=deepcopy(wrappers[0]);bad['witness']['certificate']['equalities'][0]='1'
    reject(lambda:verify_certificate(bad))
    original_catalog=api.CATALOG
    original_inputs=(original_catalog/'inputs.json').read_bytes()
    with tempfile.TemporaryDirectory(prefix='ghz30-transfer-catalog-') as directory:
        api.CATALOG=Path(directory)
        (api.CATALOG/'inputs.json').write_bytes(original_inputs)
        (api.CATALOG/'artifact_sha256.json').write_text(json.dumps({'inputs.json':'0'*64}))
        api.verified_catalog.cache_clear()
        reject(api.verified_catalog)
        altered=json.loads(original_inputs);altered['template_source_sha256']='0'*64
        raw=json.dumps(altered).encode()
        (api.CATALOG/'inputs.json').write_bytes(raw)
        (api.CATALOG/'artifact_sha256.json').write_text(json.dumps({'inputs.json':hashlib.sha256(raw).hexdigest()}))
        api.verified_catalog.cache_clear()
        reject(api.verified_catalog)
    api.CATALOG=original_catalog
    api.verified_catalog.cache_clear()
    # Capture and replay the existing exact positive primals, then check scope.
    import overnight30_screen as screen_module
    original_verify=screen_module.verify_certificate
    primals=[]
    def capture(proof):
        result=original_verify(proof)
        if proof['method']=='theta_primal':primals.append(proof)
        return result
    screen_module.verify_certificate=capture
    try:
        from event30_degree_two_cover import positive_controls
        positive_controls()
    finally:
        screen_module.verify_certificate=original_verify
    assert len(primals)==2
    assert all(api.find_transfer(proof['graph']) is None for proof in primals)

    # Importing the worker loads its existing numerical dependencies, but these
    # controls never instantiate its scorer or execute a numerical solve.
    import overnight30_asymmetric as worker
    def forbidden(*args,**kwargs):
        raise AssertionError('numerical scoring must not run in integration controls')
    worker.ThetaScorer.solve=forbidden
    worker.positive_controls=forbidden
    observed=worker.screened(worker.from_edges(wrappers[0]['graph']['edges']))
    assert observed['certificate']['method']=='vertex_probability'
    latest_dir=ROOT/'search/overnight30_reviews/2026-11-05T0225_leader'
    latest=json.loads((latest_dir/'live_best_snapshot.json').read_text())
    latest_bound=json.loads((latest_dir/'097316a813a1e997f5f6_certificate.json').read_text())
    assert api.find_transfer({'parts':latest['parts'],'edges':latest['edges']}) is None
    historical=json.loads((ROOT/'search/overnight30_reviews/2026-11-05T0100/0711be2a3c6aff655ef8_certificate.json').read_text())
    graph=historical['graph']
    partition={'g6':worker.encode(worker.from_edges(graph['edges'])),
               'edges':graph['edges'],'parts':graph['parts'],'theta_numeric':2.9166384145141824,
               'screen':{'status':'excluded'},'certificate':historical}
    higher=dict(latest,certificate=latest_bound,screen={'status':'excluded'})
    assert higher['theta_numeric']>partition['theta_numeric']
    state={'best_records':{'best_excluded_theta':higher},'retired_best_records':[partition]}
    assert worker.partition_repair_parent(state)['g6']==partition['g6']
    state=json.loads(json.dumps(state))
    state['retired_best_records']=[]
    assert worker.partition_repair_parent(state)['g6']==partition['g6']
    replayed=0
    for function in (worker.partition_witness_repair,worker.paired_partition_witness_repair):
        first,second=random.Random(1964),random.Random(1964)
        expected=function(partition,first)
        actual=function(worker.partition_repair_parent(state),second)
        assert actual is not None and actual==expected and first.getstate()==second.getstate()
        replayed+=1
    # A recovered external partition record must survive even if a current
    # screen reclassifies that same graph with another valid exact method.
    inner={'format':'ghz30-vertex-probability-farkas-v1','graph':graph,
           'mean_bound_certificates':historical['witness']['mean_bound_certificates'],
           'inequalities':[{'constraint':dict(kind='disjoint_subset',**block),'multiplier':'1/2'}
                           for block in historical['witness']['blocks']],
           'equalities':['0','-1/2','0'],'contradiction_margin':'1/32','strict_multiplier_mass':'1'}
    assert verify_inner(inner)==Fraction(1,32)
    alternate=certificate(graph,'vertex_probability',certificate=inner)
    with tempfile.TemporaryDirectory(prefix='ghz30-transfer-history-') as directory:
        folder=Path(directory);source=folder/'previous.json';prefix=folder/'new'
        source.with_name('previous_best_excluded.json').write_text(json.dumps(partition))
        migrating={'best_theta':higher['theta_numeric'],'best_excluded_theta':None,
                   'best_records':{'best_theta':higher},'current_theta':higher['theta_numeric']}
        def rescreen(adj,value):
            chosen=alternate if worker.encode(adj)==partition['g6'] else latest_bound
            verify_certificate(chosen)
            return {'status':'excluded','certificate':chosen,'reason':'exact migration control'}
        worker.recover_best_records(migrating,source,prefix,rescreen)
        assert migrating['best_excluded_theta']==higher['theta_numeric']
        assert migrating['historical_best_theta']==higher['theta_numeric']
        assert migrating['historical_best_record']['g6']==higher['g6']
        assert migrating['current_theta']==higher['theta_numeric']
        assert migrating['partition_repair_parent']['g6']==partition['g6']
        assert migrating['partition_repair_parent']['certificate']['method']=='disjoint_partition'
    result={'frozen_target_wrappers_verified':52,'source_catalog_cache_maxsize':1,
            'self_contained_replay_without_catalog':True,'wrapper_tampers_rejected':2,
            'source_hash_tampers_rejected':2,'exact_positive_primal_scope_controls':2,
            'discovery_imports_numerical_libraries':False,'numerical_solves':0,
            'latest_guide_theta':latest['theta_numeric'],'latest_guide_transferred':False,
            'higher_nonpartition_parent_fallback':True,'checkpoint_parent_persistence':True,
            'single_and_paired_rng_replays':replayed,'external_record_reclassification_preserves_partition':True,
            'historical_graph_score_and_current_score_preserved':True}
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()

"""Bounded review of a fixed published checkpoint; never writes the data repo."""
import collections
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import networkx as nx

OUT = Path(__file__).resolve().parent
SOURCE = OUT.parents[2]
DATA = Path('/workspace/ghz30-overnight')
COMMIT = 'f37683f803842de5b3a8d7fcffcbdd781bcff492'
RUN = 'search/overnight30_runs/2026-11-04'
OLD = OUT.parent/'2026-11-04T2200'
sys.path.insert(0, str(SOURCE/'search'))
from overnight30_screen import normalize_graph, graph_id, screen_graph, verify_certificate
from overnight30_fast_dual import verify_certificate as verify_fast


def main():
    start = time.monotonic()
    previous = {r['graph_id']:r for r in map(json.loads,(OUT/'strongest_candidates.jsonl').read_text().splitlines())} if (OUT/'strongest_candidates.jsonl').exists() else {}
    previous_failures = {r['graph_id']:r for r in map(json.loads,(OUT/'failure_rechecks.jsonl').read_text().splitlines())} if (OUT/'failure_rechecks.jsonl').exists() else {}
    tree = subprocess.check_output(['git','ls-tree','-r','--format=%(objectname) %(path)',COMMIT,'--',RUN],cwd=DATA,text=True)
    blobs = dict((line.split(' ',1)[1],line.split(' ',1)[0]) for line in tree.splitlines())
    inputs = {}

    def read(path):
        raw = (DATA/path).read_bytes()
        actual = hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        if actual != blobs[path]:
            raw = subprocess.check_output(['git','show',f'{COMMIT}:{path}'],cwd=DATA)
            actual = hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        assert actual == blobs[path]
        inputs[path] = {'path':path,'git_blob':actual,'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
        return raw

    manifest = json.loads(read(RUN+'/manifest.json'))
    jobs = {j['id']:j for j in manifest['jobs']}
    completed = [j for j in jobs.values() if j['state']=='completed']
    old_summary = json.loads((OLD/'summary.json').read_text())
    old_jobs = set(old_summary['completed_jobs'])
    new_jobs = [j for j in completed if j['id'] not in old_jobs]
    old_rows = list(map(json.loads,(OLD/'strongest_candidates.jsonl').read_text().splitlines()))
    old_iso = list(map(json.loads,(OLD/'isomorphic_candidates.jsonl').read_text().splitlines()))
    old_codes = {r['g6'].removeprefix('>>graph6<<').strip() for r in old_rows+old_iso}
    old_scored = set(old_codes)
    new_codes, new_scored, leaders = set(), set(), {}
    counters, error_rows, job_reports = {}, [], []

    for job in completed:
        jid = job['id']; fresh = jid not in old_jobs
        count = collections.Counter()
        archive_paths = [p for p in blobs if p.startswith(job['directory']+'/') and p.endswith('.jsonl.gz')
                         and (p.endswith('/models.jsonl.gz') or '/asymmetric_seed' in p)]
        for path in sorted(archive_paths):
            for line in gzip.decompress(read(path)).decode().splitlines():
                if not line.strip(): continue
                r = json.loads(line); count['records'] += 1
                code = r.get('g6',r.get('graph6'))
                if code:
                    code = code.removeprefix('>>graph6<<').strip()
                    (new_codes if fresh else old_codes).add(code)
                error = r.get('error')
                if error:
                    count['solver_error_records'] += 1
                    if fresh: error_rows.append({'job':jid,'family':job['family'] or 'asymmetric','error':error,'source':path,'g6':code})
                score = r.get('theta_numeric',r.get('theta'))
                if not code or not isinstance(score,(int,float)) or not math.isfinite(score): continue
                count['scored_records'] += 1
                count['solver_status:'+str(r.get('solver_status',r.get('status','unspecified')))] += 1
                (new_scored if fresh else old_scored).add(code)
                if fresh and (code not in leaders or score > leaders[code]['theta_numeric']):
                    leaders[code] = {'g6':code,'theta_numeric':score,'job':jid,'family':job['family'] or 'asymmetric',
                                     'source':path,'index':r.get('index',r.get('iteration')),
                                     'dominance_certificate_sha256':r.get('dominance_certificate_sha256')}
        if fresh:
            counters[jid] = dict(count)
            state_paths = [p for p in blobs if p.startswith(job['directory']+'/') and
                           (p.endswith('/state.json') or (job['kind']=='asymmetric' and p.endswith(f'asymmetric_seed{job["seed"]}.json')))]
            state = json.loads(read(state_paths[0])) if state_paths else {}
            keys = ('models','scored','solver_errors','solver_failures','iterations','accepted','invalid','fast_rejected',
                    'sat_excluded','rejected_pair_bridge','rejected_degree_mass','certified_dominance_cuts',
                    'batch_models','batch_scored','best_theta_numeric','best_theta','termination')
            job_reports.append({'job':jid,'family':job['family'] or 'asymmetric','resume_job':job.get('resume_job'),
                                'returncode':job.get('returncode'),'archive_counts':dict(count),
                                'checkpoint_cumulative_or_explicit_batch_fields':{k:state[k] for k in keys if k in state}})

    novel_codes = new_codes-old_codes
    novel_scored_codes = new_scored-old_scored
    ranked = sorted((r for code,r in leaders.items() if code in novel_scored_codes),key=lambda r:r['theta_numeric'],reverse=True)
    families = {}
    for r in ranked: families.setdefault(r['family'],r)
    initial = ranked[:160]+list(families.values())
    print(json.dumps({'stage':'collection','new_jobs':len(new_jobs),'new_scored_labels':len(novel_scored_codes),
                      'top_numeric':ranked[0]['theta_numeric'] if ranked else None,'seconds':time.monotonic()-start}),flush=True)

    buckets = collections.defaultdict(list)
    for r in old_rows:
        G=nx.from_graph6_bytes(r['g6'].encode())
        buckets[nx.weisfeiler_lehman_graph_hash(G)].append((r['graph_id'],G,'2026-11-04T2200'))
    seen, representatives, equivalences = set(), [], []
    for r in initial:
        if r['g6'] in seen: continue
        seen.add(r['g6'])
        graph=normalize_graph(r); assert graph['n']==30
        r={**r,'graph_id':graph_id(graph)}
        G=nx.from_graph6_bytes(r['g6'].encode()); key=nx.weisfeiler_lehman_graph_hash(G)
        match=None
        for gid,prior,origin in buckets[key]:
            matcher=nx.algorithms.isomorphism.GraphMatcher(G,prior)
            if matcher.is_isomorphic():
                match={'representative':gid,'representative_review':origin,'vertex_mapping':matcher.mapping};break
        if match:
            equivalences.append({**r,**match});continue
        buckets[key].append((r['graph_id'],G,'2026-11-04T2300'))
        representatives.append(r)

    # Cap exact numerical discovery; retain all selected family representatives.
    family_reps={}
    for r in representatives: family_reps.setdefault(r['family'],r['graph_id'])
    chosen=[r for i,r in enumerate(representatives) if i<24 or family_reps[r['family']]==r['graph_id']]
    reviewed=[]
    for r in chosen:
        if time.monotonic()-start>110: break
        if r['graph_id'] in previous and previous[r['graph_id']].get('certificate'):
            prior=previous[r['graph_id']]
            verify_certificate(json.loads((OUT/prior['certificate']).read_text()))
            reviewed.append(prior)
            continue
        graph=normalize_graph(r); result=screen_graph(graph,sat_conflicts=20000)
        digest=r.get('dominance_certificate_sha256')
        if result['status']!='excluded' and digest:
            path=str(Path(r['source']).parent/'dual_cuts'/f'{digest}.json.gz')
            if path in blobs:
                raw=gzip.decompress(read(path)); assert hashlib.sha256(raw).hexdigest()==digest
                cert=json.loads(raw); bound=verify_fast(cert)
                assert cert['edges']==graph['edges'] and cert['n']==30
                result.update(status='excluded',reason='reused published exact factor dual',
                              reused_certificate={'path':path,'commit':COMMIT,'json_sha256':digest,'bound':str(bound)})
        if result['status']!='excluded':
            result=screen_graph(graph,sat_conflicts=0,exact_dual=True)
        for field in ('certificate','dimension19_certificate'):
            if field in result:
                verify_certificate(result[field]);name=f'{r["graph_id"]}_{field}.json'
                (OUT/name).write_text(json.dumps(result[field],indent=2)+'\n');result[field]=name
        result.pop('graph',None)
        reviewed.append({**r,**result})

    failure_audit=[]
    for r in error_rows:
        if 'deadline' in str(r['error']) or not r.get('g6'): continue
        graph=normalize_graph(r);gid=graph_id(graph)
        if gid in previous_failures:
            prior=previous_failures[gid]
            verify_certificate(json.loads((OUT/prior['certificate']).read_text()))
            failure_audit.append(prior)
            continue
        result=screen_graph(graph,sat_conflicts=20000,exact_dual=True)
        for field in ('certificate','dimension19_certificate'):
            if field in result:
                verify_certificate(result[field]);name=f'{gid}_{field}.json'
                (OUT/name).write_text(json.dumps(result[field],indent=2)+'\n');result[field]=name
        result.pop('graph',None)
        failure_audit.append({**r,**result})

    summary={'published_commit':COMMIT,'baseline_review':'2026-11-04T2200',
             'completed_jobs':[j['id'] for j in completed],'new_completed_jobs':[j['id'] for j in new_jobs],
             'baseline_distinct_labeled_supports':len(old_codes),'baseline_distinct_scored_labels':len(old_scored),
             'new_chunk_distinct_labeled_supports':len(new_codes),'new_chunk_distinct_scored_labels':len(new_scored),
             'new_labeled_supports_vs_baseline':len(novel_codes),'new_scored_labels_vs_baseline':len(novel_scored_codes),
             'new_raw_records':sum(c.get('records',0) for c in counters.values()),
             'new_scored_records':sum(c.get('scored_records',0) for c in counters.values()),
             'selected_distinct_labeled_supports':len(seen),'selected_isomorphism_equivalences':len(equivalences),
             'selected_new_isomorphism_representatives':len(representatives),'exactly_reviewed_new_representatives':len(reviewed),
             'review_statuses':dict(collections.Counter(r['status'] for r in reviewed)),
             'failure_recheck_statuses':dict(collections.Counter(r['status'] for r in failure_audit)),
             'solver_error_counts':dict(collections.Counter(str(r['error']) for r in error_rows)),
             'numerical_status_counts':dict(sum((collections.Counter({k.removeprefix('solver_status:'):v for k,v in c.items()
                                                                      if k.startswith('solver_status:')}) for c in counters.values()),collections.Counter())),
             'family_highest_new_scored_labels':families,'seconds':time.monotonic()-start,
             'scope':'uncolored exact isomorphism comparison only within selected leaders and against prior selected representatives; no entire-pool isomorphism census',
             'counter_scope':'archive counts describe new chunks; checkpoint fields are cumulative unless explicitly named batch_*'}
    for name,rows in [('strongest_candidates.jsonl',reviewed),('selected_representatives.jsonl',representatives),
                      ('isomorphic_candidates.jsonl',equivalences),('inputs.jsonl',list(inputs.values())),
                      ('job_audit.jsonl',job_reports),('solver_errors.jsonl',error_rows),('failure_rechecks.jsonl',failure_audit)]:
        (OUT/name).write_text(''.join(json.dumps(r,separators=(',',':'))+'\n' for r in rows))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if (OUT/'artifact_sha256.json').exists():
        hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.iterdir())
                if p.is_file() and p.name!='artifact_sha256.json'}
        (OUT/'artifact_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='family_highest_new_scored_labels'}),flush=True)
    print('Top reviewed:',[(r['family'],r['theta_numeric'],r['status'],r.get('reason')) for r in reviewed[:5]],flush=True)


if __name__=='__main__':
    main()

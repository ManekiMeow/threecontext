"""Materialize and verify natural degree-two closure cuts from the fixed snapshot."""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use -O')
import collections
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from overnight30_symmetry import generator, replay_structural_dominance


def orbits(family):
    permutation=generator(family)
    result=[]; seen=set()
    for u in range(30):
        for v in range(u+1,30):
            if u//10==v//10 or (u,v) in seen: continue
            orbit=[]; edge=(u,v)
            while edge not in orbit:
                orbit.append(edge)
                edge=tuple(sorted((permutation[edge[0]],permutation[edge[1]])))
            seen.update(orbit);result.append(sorted(orbit))
    return result


def materialize():
    data=Path('/workspace/ghz30-overnight')
    summary=json.loads((HERE/'summary.json').read_text())
    commit=summary['published_commit']
    manifest=json.loads(subprocess.check_output(['git','show',commit+':search/overnight30_runs/2026-11-04/manifest.json'],cwd=data))
    inputs={r['path']:r for r in map(json.loads,(HERE/'inputs.jsonl').read_text().splitlines())}
    result=[]
    for job in manifest['jobs']:
        if job['id'] not in summary['new_completed_jobs'] or job['kind']!='symmetry': continue
        path=f"{job['directory']}/{job['family']}/seed_{job['seed']}/cuts.jsonl.gz"
        raw=subprocess.check_output(['git','show',commit+':'+path],cwd=data)
        selected=[]
        for index,line in enumerate(gzip.decompress(raw).splitlines()):
            row=json.loads(line)
            if (isinstance(row,dict) and isinstance(row.get('witness'),dict)
                    and row['witness'].get('method')=='degree_two_closure'):
                selected.append({'job':job['id'],'family':job['family'],'source':path,'source_line':index,'cut':row})
        if selected:
            inputs[path]={'path':path,'git_blob':hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest(),
                          'sha256':hashlib.sha256(raw).hexdigest(),'bytes':len(raw)}
            result.extend(selected)
    (HERE/'closure_cuts.jsonl.gz').write_bytes(gzip.compress((''.join(json.dumps(r,separators=(',',':'))+'\n' for r in result)).encode(),mtime=0))
    (HERE/'inputs.jsonl').write_text(''.join(json.dumps(r,separators=(',',':'))+'\n' for r in inputs.values()))
    report=verify()
    (HERE/'closure_cut_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def verify(inputs=False):
    result=list(map(json.loads,gzip.decompress((HERE/'closure_cuts.jsonl.gz').read_bytes()).splitlines()))
    counts=collections.Counter();sizes=collections.Counter(); cache={}
    source_rows={}
    if inputs:
        summary=json.loads((HERE/'summary.json').read_text())
        for path in {row['source'] for row in result}:
            raw=subprocess.check_output(['git','show',summary['published_commit']+':'+path],cwd='/workspace/ghz30-overnight')
            source_rows[path]=list(map(json.loads,gzip.decompress(raw).splitlines()))
    for row in result:
        edge_orbits=cache.setdefault(row['family'],orbits(row['family']))
        assert replay_structural_dominance(row['cut'],edge_orbits)==row['cut']['clause']
        counts[str(row['job'])]+=1;sizes[str(len(row['cut']['witness']['witness']['vertices']))]+=1
        if inputs: assert source_rows[row['source']][row['source_line']]==row['cut']
    assert len(result)==281 and dict(counts)=={'93':50,'94':225,'96':6}
    return {'exact_compact_cuts':len(result),'by_job':dict(counts),'closed_subset_sizes':dict(sizes),
            'scope':'Natural completed-worker cuts; exact certificates and dominance clauses, not independent isomorphism classes.'}


if __name__=='__main__':
    print(json.dumps(materialize() if '--materialize' in sys.argv else verify('--inputs' in sys.argv),indent=2))

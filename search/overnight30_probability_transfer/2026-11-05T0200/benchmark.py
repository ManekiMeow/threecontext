"""Bounded exact-only transfer benchmark on frozen controls and new chunk98 labels."""
from collections import Counter,defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import time
from transfer import prepare_template,transfer,geometry,ROOT

HERE=Path(__file__).resolve().parent
MAIN=ROOT/'search/overnight30_probability_reviews/2026-11-05T0200'
graph_from_g6=runpy.run_path(str(MAIN/'verify.py'))['graph_from_g6']
from overnight30_screen import fast_screen,graph_id


def rows(path):
    with gzip.open(path,'rt') as stream:
        return list(map(json.loads,stream))


def write_rows(path,data):
    with gzip.GzipFile(filename=str(path),mode='wb',mtime=0) as stream:
        stream.write(''.join(json.dumps(row,separators=(',',':'))+'\n' for row in data).encode())


def main():
    if (HERE/'summary.json').exists():
        raise RuntimeError('Immutable benchmark already exists; use a new directory.')
    packed=rows(MAIN/'certificates.jsonl.gz')
    groups=defaultdict(list)
    for row in packed:
        key=json.dumps({k:v for k,v in row['certificate'].items() if k!='graph'},sort_keys=True)
        groups[key].append(row)
    sources=[group[0] for group in sorted(groups.values(),key=lambda group:(-len(group),group[0]['name']))]
    templates=[dict(name=row['name'],prepared=prepare_template(row['certificate'])) for row in sources]
    for template in templates:
        original=template['prepared']['proof']
        reproduced,details=transfer(template['prepared'],original['graph'])
        assert reproduced==original and details['added_upper_rows']==0
    exclusions,bindings=set(),{}
    for folder in sorted((ROOT/'search/overnight30_reviews').iterdir()):
        if not folder.is_dir():
            continue
        for name in ('selected_representatives.jsonl','strongest_candidates.jsonl','isomorphic_candidates.jsonl'):
            path=folder/name
            if path.exists():
                raw=path.read_bytes()
                bindings[str(path.relative_to(ROOT))]=hashlib.sha256(raw).hexdigest()
                exclusions.update(row['g6'] for row in map(json.loads,raw.splitlines()) if 'g6' in row)
    old=rows(MAIN/'results.jsonl.gz')
    exclusions.update(row['g6'] for row in old)
    controls=[row for row in old if not row['new_exact_exclusion'] and not row['already_structurally_excluded']][:4]
    commit='3536304a360437b6c54413602377c31438c719f0'
    source='search/overnight30_runs/2026-11-04/chunks/00098_asymmetric/asymmetric_seed20261104.jsonl.gz'
    raw=subprocess.check_output(['git','show',f'{commit}:{source}'],cwd='/workspace/ghz30-overnight')
    distinct={}
    scored=0
    for index,line in enumerate(gzip.decompress(raw).splitlines()):
        row=json.loads(line)
        if 'theta' not in row:
            continue
        scored+=1
        if row['g6'] in exclusions:
            continue
        candidate={'g6':row['g6'],'theta_numeric':row['theta'],'source_line':index,
                   'iteration':row['iteration'],'graph_id':graph_id(graph_from_g6(row['g6']))}
        distinct.setdefault(row['g6'],candidate)
    sample=sorted(distinct.values(),key=lambda row:(-row['theta_numeric'],row['iteration'],row['g6']))[:120]
    inputs={'published_commit':commit,'source':source,'source_sha256':hashlib.sha256(raw).hexdigest(),
            'source_scored_records':scored,'eligible_unreviewed_distinct_labels':len(distinct),
            'excluded_previous_review_labels':len(exclusions),'exclusion_sources':bindings,
            'template_source':str((MAIN/'certificates.jsonl.gz').relative_to(ROOT)),
            'template_source_sha256':hashlib.sha256((MAIN/'certificates.jsonl.gz').read_bytes()).hexdigest(),
            'template_names':[item['name'] for item in templates],
            'controls':controls,'new_sample':sample,
            'selection':'Top120 scores after exact graph6 deduplication and exclusion of all saved review labels; no isomorphism claim.'}
    (HERE/'inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    outputs,proofs=[],[]
    started=time.monotonic()
    for kind,candidates in [('control',controls),('new_sample',sample)]:
        limit=started+(2.0 if kind=='control' else 14.8)
        for row in candidates:
            if time.monotonic()>=limit:
                break
            begin=time.monotonic()
            graph=graph_from_g6(row['g6'])
            screen=fast_screen(graph)
            result=dict(row,sample_kind=kind,existing_exclusion=screen['status']=='excluded',
                        transferred=False,attempts=0,attempt_outcomes={})
            if screen['status']=='excluded':
                result['existing_method']=screen['certificate']['method']
            else:
                base=geometry(graph)
                outcomes=Counter()
                for template in templates:
                    if time.monotonic()>=limit:
                        break
                    proof,details=transfer(template['prepared'],graph,base)
                    result['attempts']+=1
                    if proof is None:
                        outcomes[details]+=1
                        continue
                    result.update(transferred=True,source_template=template['name'],
                                  contradiction_margin=proof['contradiction_margin'],details=details)
                    proofs.append({'name':row['graph_id'],'source_template':template['name'],
                                   'sample_kind':kind,'certificate':proof})
                    outcomes['certified']+=1
                    break
                result['attempt_outcomes']=dict(outcomes)
                result['all_templates_tried']=result['attempts']==len(templates)
            result['transfer_seconds']=time.monotonic()-begin
            outputs.append(result)
    elapsed=time.monotonic()-started
    write_rows(HERE/'results.jsonl.gz',outputs)
    write_rows(HERE/'certificates.jsonl.gz',proofs)
    summary={'scope':'Exact-only bounded certificate transfer; no numerical solver and no worker integration.',
             'requested_seconds':15,'elapsed_seconds':elapsed,'templates':len(templates),
             'template_self_controls':len(templates),'sets':{}}
    for kind,candidates in [('control',controls),('new_sample',sample)]:
        selected=[row for row in outputs if row['sample_kind']==kind]
        summary['sets'][kind]={'selected':len(candidates),'processed':len(selected),
                              'existing_structural_exclusions':sum(row['existing_exclusion'] for row in selected),
                              'transferred_exact_exclusions':sum(row['transferred'] for row in selected),
                              'template_attempts':sum(row['attempts'] for row in selected),
                              'complete_unsuccessful_attempts':sum(not row['existing_exclusion'] and not row['transferred'] and row['all_templates_tried'] for row in selected),
                              'unreached':len(candidates)-len(selected)}
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    names=['inputs.json','results.jsonl.gz','certificates.jsonl.gz','summary.json']
    (HERE/'artifact_sha256.json').write_text(json.dumps({name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()

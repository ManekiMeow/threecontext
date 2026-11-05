"""Exact proof, provenance and optional immutable-input replay for transfer review."""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use -O.')
import argparse
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
from transfer import ROOT,prepare_template,transfer,verify_certificate

HERE=Path(__file__).resolve().parent
MAIN=ROOT/'search/overnight30_probability_reviews/2026-11-05T0200'
graph_from_g6=runpy.run_path(str(MAIN/'verify.py'))['graph_from_g6']


def rows(path):
    with gzip.open(path,'rt') as stream:
        return list(map(json.loads,stream))


def main(inputs=False):
    for name,digest in json.loads((HERE/'artifact_sha256.json').read_text()).items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
    for name,digest in json.loads((HERE/'audit_artifacts_sha256.json').read_text()).items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest
    source=json.loads((HERE/'inputs.json').read_text())
    summary=json.loads((HERE/'summary.json').read_text())
    original_path=ROOT/source['template_source']
    assert hashlib.sha256(original_path.read_bytes()).hexdigest()==source['template_source_sha256']
    original={row['name']:row['certificate'] for row in rows(original_path)}
    templates={name:prepare_template(original[name]) for name in source['template_names']}
    assert len(templates)==summary['templates']==summary['template_self_controls']==48
    for template in templates.values():
        proof=template['proof']
        reproduced,details=transfer(template,proof['graph'])
        assert reproduced==proof and details['added_upper_rows']==0
    known_exclusions=set()
    for path,digest in source['exclusion_sources'].items():
        raw=(ROOT/path).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==digest
        known_exclusions.update(row['g6'] for row in map(json.loads,raw.splitlines()) if 'g6' in row)
    old=rows(MAIN/'results.jsonl.gz')
    known_exclusions.update(row['g6'] for row in old)
    assert len(known_exclusions)==source['excluded_previous_review_labels']
    assert source['controls']==[row for row in old if not row['new_exact_exclusion'] and not row['already_structurally_excluded']][:4]
    assert len(source['new_sample'])==120
    assert len({row['g6'] for row in source['new_sample']})==120
    assert all(row['g6'] not in known_exclusions for row in source['new_sample'])
    results=rows(HERE/'results.jsonl.gz')
    packed=rows(HERE/'certificates.jsonl.gz')
    proofs={row['name']:row for row in packed}
    assert len(proofs)==len(packed)
    existing_packed=rows(HERE/'existing_certificates.jsonl.gz')
    existing={row['name']:row['certificate'] for row in existing_packed}
    assert len(existing)==len(existing_packed)==6
    existing_used=set()
    used=set()
    for kind,key in [('control','controls'),('new_sample','new_sample')]:
        selected=[row for row in results if row['sample_kind']==kind]
        record=summary['sets'][kind]
        assert record['selected']==len(source[key])
        assert record['processed']==len(selected)
        assert record['unreached']==len(source[key])-len(selected)
        assert record['existing_structural_exclusions']==sum(row['existing_exclusion'] for row in selected)
        assert record['transferred_exact_exclusions']==sum(row['transferred'] for row in selected)
        assert record['template_attempts']==sum(row['attempts'] for row in selected)
        assert record['complete_unsuccessful_attempts']==sum(not row['existing_exclusion'] and not row['transferred'] and row['all_templates_tried'] for row in selected)
        for expected,row in zip(source[key],selected):
            assert all(row[k]==value for k,value in expected.items())
            assert sum(row['attempt_outcomes'].values())==row['attempts']
            if row['existing_exclusion']:
                from overnight30_screen import verify_certificate as verify_global
                certificate=existing[row['graph_id']]
                assert certificate['graph']==graph_from_g6(row['g6'])
                assert certificate['method']==row['existing_method']
                verify_global(certificate)
                existing_used.add(row['graph_id'])
            if row['transferred']:
                saved=proofs[row['graph_id']]
                assert saved['sample_kind']==kind and saved['source_template']==row['source_template']
                assert saved['certificate']['graph']==graph_from_g6(row['g6'])
                verify_certificate(saved['certificate'])
                regenerated,details=transfer(templates[row['source_template']],saved['certificate']['graph'])
                assert regenerated==saved['certificate'] and details==row['details']
                assert saved['certificate']['contradiction_margin']==row['contradiction_margin']
                used.add(row['graph_id'])
    assert used==set(proofs)
    assert existing_used==set(existing)
    if inputs:
        from overnight30_screen import graph_id
        raw=subprocess.check_output(['git','show',f"{source['published_commit']}:{source['source']}"],cwd='/workspace/ghz30-overnight')
        assert hashlib.sha256(raw).hexdigest()==source['source_sha256']
        distinct={};scored=0
        for index,line in enumerate(gzip.decompress(raw).splitlines()):
            row=json.loads(line)
            if 'theta' not in row:
                continue
            scored+=1
            if row['g6'] in known_exclusions:
                continue
            candidate={'g6':row['g6'],'theta_numeric':row['theta'],'source_line':index,
                       'iteration':row['iteration'],'graph_id':graph_id(graph_from_g6(row['g6']))}
            distinct.setdefault(row['g6'],candidate)
        assert scored==source['source_scored_records']
        assert len(distinct)==source['eligible_unreviewed_distinct_labels']
        assert source['new_sample']==sorted(distinct.values(),key=lambda row:(-row['theta_numeric'],row['iteration'],row['g6']))[:120]
    print(f"Verified {len(proofs)} exact transferred certificates, graph bindings and benchmark counts.")


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',action='store_true')
    main(parser.parse_args().inputs)

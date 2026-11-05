"""Reproduce the bounded singleton probability-LP trial; never overwrite evidence."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
from event30_vertex_probability_lp import certify, verify_certificate, self_test
from overnight30_screen import fast_screen


def write_rows(path, rows):
    with gzip.GzipFile(filename=str(path), mode='wb', mtime=0) as stream:
        stream.write(''.join(json.dumps(row, separators=(',', ':'))+'\n' for row in rows).encode())


def main(seconds):
    if (HERE/'summary.json').exists():
        raise RuntimeError('This immutable review already exists; use a new review directory.')
    source = ROOT/'search/overnight30_reviews/2026-11-05T0200/selected_representatives.jsonl'
    source_raw = source.read_bytes()
    inputs = [json.loads(line) for line in source_raw.splitlines()]
    (HERE/'inputs.jsonl').write_bytes(source_raw)
    started = time.monotonic()
    rows, proofs, kinds = [], [], Counter()
    for index, record in enumerate(inputs):
        remaining = seconds-(time.monotonic()-started)
        if remaining < .2:
            break
        begin = time.monotonic()
        screening = fast_screen(record)
        if screening['status'] == 'excluded':
            row = dict(record, selection_index=index, already_structurally_excluded=True,
                       existing_method=screening['certificate']['method'], new_exact_exclusion=False)
        else:
            proof = certify(record, denominator=4, refinements=5, time_limit=min(3.0,remaining))
            row = dict(record, selection_index=index, already_structurally_excluded=False,
                       new_exact_exclusion=proof is not None)
            if proof is not None:
                verify_certificate(proof)
                name = record['graph_id']
                proofs.append({'name':name, 'certificate':proof})
                row.update(certificate=name, contradiction_margin=proof['contradiction_margin'],
                           strict_multiplier_mass=proof['strict_multiplier_mass'])
                kinds.update(term['constraint']['kind'] for term in proof['inequalities'])
        row['seconds'] = time.monotonic()-begin
        rows.append(row)
        print(json.dumps({'index': index, 'id':record['graph_id'],
                          'excluded':row['new_exact_exclusion'],
                          'skipped':row['already_structurally_excluded'],
                          'seconds':round(row['seconds'],3)}), flush=True)
    elapsed = time.monotonic()-started
    write_rows(HERE/'results.jsonl.gz', rows)
    write_rows(HERE/'certificates.jsonl.gz', proofs)
    summary = {'scope':'Bounded exact probability relaxation on selected supports; no feasibility claims.',
               'published_cutoff':'3536304a360437b6c54413602377c31438c719f0',
               'selected_source':str(source.relative_to(ROOT)),
               'selected_source_sha256':hashlib.sha256(source_raw).hexdigest(),
               'selected_input_count':len(inputs), 'reviewed':len(rows),
               'existing_structural_exclusions':sum(row['already_structurally_excluded'] for row in rows),
               'lp_attempted':sum(not row['already_structurally_excluded'] for row in rows),
               'new_exact_exclusions':len(proofs), 'unreached':len(inputs)-len(rows),
               'requested_seconds':seconds, 'elapsed_seconds':elapsed,
               'constraint_kinds_used':dict(kinds), 'controls':self_test(),
               'results_archive':'results.jsonl.gz','proof_archive':'certificates.jsonl.gz'}
    audit_source=source.with_name('strongest_candidates.jsonl')
    audit_raw=audit_source.read_bytes()
    earlier={row['graph_id'] for row in map(json.loads,audit_raw.splitlines()) if row['status']=='excluded'}
    certified={row['graph_id'] for row in rows if row['new_exact_exclusion']}
    summary['prior_selected_exclusion_source']=str(audit_source.relative_to(ROOT))
    summary['prior_selected_exclusion_source_sha256']=hashlib.sha256(audit_raw).hexdigest()
    summary['additional_scalar_certificates']=len(certified)
    summary['overlap_with_prior_selected_exclusions']=len(certified & earlier)
    summary['newly_excluded_selected_representatives']=len(certified-earlier)
    summary['combined_selected_exclusions']=len(certified | earlier)
    summary['unresolved_selected_representatives']=len(inputs)-len(certified | earlier)
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    names=['inputs.jsonl','results.jsonl.gz','certificates.jsonl.gz','summary.json']
    (HERE/'artifact_sha256.json').write_text(json.dumps({name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                                       for name in names},indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds',type=float,default=60)
    main(parser.parse_args().seconds)

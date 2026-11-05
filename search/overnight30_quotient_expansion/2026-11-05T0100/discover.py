"""Bounded analytic review of published C5 quotient classes outside all previously reviewed classes."""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import runpy
import subprocess
import sys
import time

SEARCH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(SEARCH))
from event30_disjoint_probability import certify_quotient, verify_quotient_certificate
from event30_probability_lp import verify_certificate as verify_lp
from event30_quotient_saturation import derive, verify_certificate as verify_saturation
from overnight30_screen import fast_screen, verify_certificate as verify_graph


def pack_certificates(out, summary):
    names = sorted({row[key] for row in summary['results'] for key in ('certificate', 'representative_proof') if key in row})
    archive = 'certificates.jsonl.gz'
    with gzip.open(out/archive, 'wt') as stream:
        for name in names:
            certificate = json.loads((out/name).read_text())
            stream.write(json.dumps({'name': name, 'certificate': certificate}, separators=(',', ':'))+'\n')
    for name in names:
        (out/name).unlink()
    summary['proof_archive'] = archive
    summary['proof_count'] = len(names)
    summary['proof_kind_counts'] = {kind: sum(row.get('proof_kind') == kind for row in summary['results'])
                                    for kind in ('probability', 'saturation', 'disjoint_probability')}


def pack_results(out, summary):
    archive = 'results.jsonl.gz'
    with gzip.open(out/archive, 'wt') as stream:
        for row in summary.pop('results'):
            stream.write(json.dumps(row, separators=(',', ':'))+'\n')
    summary['results_archive'] = archive


def attach_guide_certificate(out, summary, data_repo):
    guide = summary['best_surviving_guide']
    if guide is None or 'dominance_certificate_sha256' not in guide:
        return
    from overnight30_fast_dual import verify_certificate as verify_dual
    from overnight30_quotient import graph_data
    digest = guide['dominance_certificate_sha256']
    source = data_repo/Path(guide['source']).parent/'dual_cuts'/(digest+'.json.gz')
    packed = source.read_bytes()
    raw = gzip.decompress(packed)
    assert hashlib.sha256(raw).hexdigest() == digest
    certificate = json.loads(raw)
    assert str(verify_dual(certificate)) == guide['exact_theta_upper']
    assert certificate['edges'] == graph_data(guide)['edges']
    name = 'best_guide_individual_dual.json.gz'
    (out/name).write_bytes(packed)
    guide['individual_certificate'] = name
    guide['role'] = 'mutation guide; support individually excluded, quotient unresolved by this analytic review'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-repo', type=Path, default=Path('/workspace/ghz30-overnight'))
    parser.add_argument('--commit', default='5a2b45478d066d5180d1364259e1c7b1d233fc33')
    parser.add_argument('--seconds', type=float, default=90)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    if (out/'summary.json').exists():
        raise RuntimeError('Refusing to overwrite a completed review; choose a new output directory.')
    prior_dir = SEARCH/'overnight30_quotient_reviews'/'2026-11-04T2300'
    helpers = runpy.run_path(str(prior_dir/'verify.py'))
    graph_quotient, canonical = helpers['graph_quotient'], helpers['canonical']
    prior = json.loads((prior_dir/'summary.json').read_text())
    already = {canonical(item['quotient']) for item in prior['results']}
    baseline_dir = SEARCH/'overnight30_quotient_expansion'/'2026-11-05T0000'
    baseline = json.loads((baseline_dir/'summary.json').read_text())
    with gzip.open(baseline_dir/baseline['results_archive'], 'rt') as stream:
        already.update(canonical(json.loads(line)['quotient']) for line in stream)
    mean_path = SEARCH/'overnight30_disjoint_reviews'/'2026-11-05T0000'/'disjoint_05_mean_bound.json'
    mean = json.loads(mean_path.read_text())
    files = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', args.commit],
                                    cwd=args.data_repo, text=True).splitlines()
    files = [p for p in files if '_c5_double/' in p and p.endswith(('models.jsonl', 'models.jsonl.gz'))]
    pool, cache, scored = {}, {}, 0
    for relative in files:
        path = args.data_repo/relative
        with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as stream:
            for line in stream:
                record = json.loads(line)
                if record.get('status') != 'scored_numerically':
                    continue
                scored += 1
                quotient = graph_quotient(record['graph6'])
                raw = tuple(sum(quotient, []))
                if raw not in cache:
                    cache[raw] = canonical(quotient)
                key = cache[raw]
                if key not in pool or record['theta_numeric'] > pool[key]['theta_numeric']:
                    pool[key] = {'theta_numeric': record['theta_numeric'], 'graph6': record['graph6'],
                                 'source': relative, 'index': record['index'], 'quotient': quotient}
                    for field in ('exact_theta_upper', 'dominance_certificate_sha256'):
                        if field in record:
                            pool[key][field] = record[field]
    ordered = sorted(pool.values(), key=lambda item: -item['theta_numeric'])
    with gzip.open(out/'published_pool.jsonl.gz', 'wt') as stream:
        for record in ordered:
            stream.write(json.dumps(record,separators=(',', ':'))+'\n')
    print(json.dumps({'stage': 'collected', 'scored_records': scored, 'quotient_classes': len(pool),
                      'published_commit': args.commit}), flush=True)
    rows, count, started, candidate = [], 0, time.monotonic(), None
    for item in ordered:
        if canonical(item['quotient']) in already:
            continue
        if time.monotonic()-started >= args.seconds:
            break
        count += 1
        result = {**item, 'exact_quotient_excluded': False}
        saturation = derive(item['quotient'])
        certificate, kind = None, None
        if saturation['status'] == 'excluded':
            certificate, kind = saturation['certificate'], 'saturation'
            verify_saturation(certificate)
        else:
            certificate = certify_quotient(item['quotient'], mean, refinements=5)
            if certificate is not None:
                eligible = {i for i, _ in certificate['disjoint_contexts']}
                relies_on_moment = any(term['constraint']['kind'] == 'upper'
                                      and term['constraint']['orbit'] in eligible
                                      for term in certificate['inequalities'])
                if relies_on_moment:
                    kind = 'disjoint_probability'
                    verify_quotient_certificate(certificate)
                else:
                    certificate = {k: v for k, v in certificate.items()
                                   if k not in ('mean_bound_certificate', 'disjoint_contexts')}
                    certificate['format'] = 'ghz30-orbit-probability-farkas-v1'
                    kind = 'probability'
                    verify_lp(certificate)
        if certificate is not None:
            digest = hashlib.sha256(json.dumps(certificate, sort_keys=True).encode()).hexdigest()[:20]
            name = digest+'_certificate.json'
            (out/name).write_text(json.dumps(certificate, indent=2)+'\n')
            result.update(exact_quotient_excluded=True, proof_kind=kind, certificate=name)
        else:
            screened = fast_screen(item)
            result['representative_screen'] = screened['status']
            if screened['status'] == 'excluded':
                graph_certificate = screened['certificate']
                verify_graph(graph_certificate)
                digest = hashlib.sha256(json.dumps(graph_certificate, sort_keys=True).encode()).hexdigest()[:20]
                name = digest+'_graph_certificate.json'
                (out/name).write_text(json.dumps(graph_certificate, indent=2)+'\n')
                result.update(representative_proof=name, representative_method=graph_certificate['method'])
            elif candidate is None:
                candidate = dict(result)
                print(json.dumps({'stage': 'best_surviving_guide', **candidate}), flush=True)
        rows.append(result)
        if count % 250 == 0:
            print(json.dumps({'stage': 'progress', 'reviewed': count,
                              'exact_quotient_exclusions': sum(r['exact_quotient_excluded'] for r in rows),
                              'elapsed': time.monotonic()-started}), flush=True)
    summary = {'published_commit': args.commit, 'source_files': files, 'scored_records': scored,
               'published_quotient_classes': len(pool), 'prior_reviewed_classes_skipped': len(already),
               'reviewed': len(rows), 'certified': sum(row['exact_quotient_excluded'] for row in rows),
               'inconclusive_quotients': sum(not row['exact_quotient_excluded'] for row in rows),
               'discovery_seconds': time.monotonic()-started, 'requested_seconds': args.seconds,
               'best_surviving_guide': candidate,
               'scope': 'Bounded review of published scored C5 quotients, excluding the prior150selected plus3560reviewed classes. '
                        'No claim of complete C5-family coverage; numerical scores and failed certificate discovery '
                        'do not establish quantum feasibility or nonexistence.', 'results': rows}
    pack_certificates(out, summary)
    attach_guide_certificate(out, summary, args.data_repo)
    pack_results(out, summary)
    summary.update(previous_selected_catalog=str(prior_dir.relative_to(SEARCH.parent)/'summary.json'),
                   previous_expansion_catalog=str(baseline_dir.relative_to(SEARCH.parent)/'summary.json'),
                   published_pool_archive='published_pool.jsonl.gz',
                   published_pool_exactly_covered_by_reviews=len(rows)+len(already)==len(pool),
                   combined_certified_before_new_phase_proofs=150+baseline['certified']+summary['certified'],
                   combined_residual_before_new_phase_proofs=baseline['inconclusive_quotients']+summary['inconclusive_quotients'],
                   source_sha256={name: hashlib.sha256((args.data_repo/name).read_bytes()).hexdigest() for name in files})
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('results', 'best_surviving_guide')}), flush=True)


if __name__ == '__main__':
    main()

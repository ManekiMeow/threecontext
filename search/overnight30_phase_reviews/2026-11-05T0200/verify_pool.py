"""Verify incremental published-pool counts and seed bindings from pinned Git blobs."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys

if not __debug__:
    raise RuntimeError('Exact source verification requires assertions; do not use -O')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
from overnight30_fast_dual import verify_certificate as verify_dual
from overnight30_quotient import graph_data


def main(data_repo):
    report = json.loads((HERE/'pool_coverage.json').read_text())
    helpers = runpy.run_path(str(ROOT/'search/overnight30_quotient_reviews/2026-11-04T2300/verify.py'))
    canonical, graph_quotient = helpers['canonical'], helpers['graph_quotient']
    prior = ROOT/report['previous_pool']
    assert hashlib.sha256(prior.read_bytes()).hexdigest() == report['previous_pool_sha256']
    with gzip.open(prior, 'rt') as stream:
        pool = {canonical(row['quotient']): row for row in map(json.loads, stream)}
    assert len(pool) == report['previous_pool_classes'] == 5589
    old = set(pool)
    raw = subprocess.check_output(['git', 'show', report['published_commit']+':'+report['source']], cwd=data_repo)
    assert hashlib.sha256(raw).hexdigest() == report['source_sha256']
    scored, seen, cache = 0, set(), {}
    for line in gzip.decompress(raw).splitlines():
        row = json.loads(line)
        if row.get('status') != 'scored_numerically':
            continue
        scored += 1
        quotient = graph_quotient(row['graph6'])
        raw_key = tuple(sum(quotient, []))
        if raw_key not in cache:
            cache[raw_key] = canonical(quotient)
        key = cache[raw_key]
        seen.add(key)
        if key not in pool or row['theta_numeric'] > pool[key]['theta_numeric']:
            pool[key] = dict(row, quotient=quotient, source=report['source'])
    assert scored == report['scored_records'] == 2540
    assert len(seen) == report['quotient_classes'] == 149
    assert len(seen-old) == report['new_quotient_classes'] == 34
    assert len(pool) == report['updated_pool_classes'] == 5623
    selected = json.loads((ROOT/'search/overnight30_quotient_reviews/2026-11-04T2300/summary.json').read_text())
    excluded = {canonical(row['quotient']) for row in selected['results']}
    for hour in ('0000', '0100'):
        with gzip.open(ROOT/f'search/overnight30_quotient_expansion/2026-11-05T{hour}/results.jsonl.gz', 'rt') as stream:
            excluded |= {canonical(row['quotient']) for row in map(json.loads, stream) if row['exact_quotient_excluded']}
    for hour in ('0100', '0200'):
        summary = json.loads((ROOT/f'search/overnight30_phase_reviews/2026-11-05T{hour}/summary.json').read_text())
        excluded |= {canonical(row['quotient']) for row in summary['results']}
    assert len(excluded) == report['known_excluded_quotient_classes'] == 5387
    assert excluded <= old
    remaining = sorted((row for key, row in pool.items() if key not in excluded), key=lambda row: -row['theta_numeric'])
    assert len(remaining) == report['updated_pool_remaining'] == 236
    integration = json.loads((HERE/'integration_checks.json').read_text())
    guide = remaining[0]
    assert guide['theta_numeric'] == integration['next_guide_theta']
    assert guide['source'] == integration['next_guide_source']
    assert guide['index'] == integration['next_guide_index']
    source = subprocess.check_output(['git', 'show', report['published_commit']+':'+guide['source']], cwd=data_repo)
    row = json.loads(gzip.decompress(source).splitlines()[integration['next_guide_line_number_zero_based']])
    assert row['index'] == guide['index'] and row['graph6'] == guide['graph6'] and row['theta_numeric'] == guide['theta_numeric']
    proof = json.loads(gzip.decompress((HERE/'best_guide_individual_dual.json.gz').read_bytes()))
    assert proof['edges'] == graph_data(guide)['edges']
    assert str(verify_dual(proof)) == integration['next_guide_exact_upper']
    print('Verified 2540 new scored records, 34 additional quotient classes, updated 5623-class pool, and exact mutation-guide binding.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-repo', type=Path, default=Path('/workspace/ghz30-overnight'))
    main(parser.parse_args().data_repo)

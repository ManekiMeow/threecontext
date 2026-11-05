"""Exact replay and complete published-pool quotient coverage for the 01:00 review."""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use -O.')

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys

HERE = Path(__file__).resolve().parent
SEARCH = HERE.parents[1]
ROOT = SEARCH.parent
sys.path.insert(0, str(SEARCH))
from event30_disjoint_probability import verify_quotient_certificate as verify_disjoint
from event30_probability_lp import verify_certificate as verify_probability
from event30_quotient_saturation import verify_certificate as verify_saturation
from overnight30_fast_dual import verify_certificate as verify_dual
from overnight30_quotient import graph_data
from overnight30_screen import verify_certificate as verify_graph


def read_rows(path):
    with gzip.open(path, 'rt') as stream:
        return [json.loads(line) for line in stream]


def main(inputs=False, data_repo=Path('/workspace/ghz30-overnight')):
    hashes = json.loads((HERE/'artifact_sha256.json').read_text())
    for name, digest in hashes.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest
    summary = json.loads((HERE/'summary.json').read_text())
    results = read_rows(HERE/summary['results_archive'])
    packed = read_rows(HERE/summary['proof_archive'])
    proofs = {row['name']: row['certificate'] for row in packed}
    assert len(proofs) == len(packed) == summary['proof_count']
    first = ROOT/summary['previous_selected_catalog']
    helpers = runpy.run_path(str(first.with_name('verify.py')))
    graph_quotient, canonical = helpers['graph_quotient'], helpers['canonical']
    old = {canonical(row['quotient']) for row in json.loads(first.read_text())['results']}
    second = ROOT/summary['previous_expansion_catalog']
    prior = json.loads(second.read_text())
    prior_results = read_rows(second.parent/prior['results_archive'])
    for row in prior_results:
        key = canonical(row['quotient'])
        assert key not in old
        old.add(key)
    assert len(old) == summary['prior_reviewed_classes_skipped'] == 3710
    seen, used, kinds = set(), set(), {}
    verifiers = {'probability': verify_probability, 'saturation': verify_saturation,
                 'disjoint_probability': verify_disjoint}
    for row in results:
        assert graph_quotient(row['graph6']) == row['quotient']
        key = canonical(row['quotient'])
        assert key not in old and key not in seen
        seen.add(key)
        if row['exact_quotient_excluded']:
            name, kind = row['certificate'], row['proof_kind']
            proof = proofs[name]
            assert proof['quotient'] == row['quotient']
            verifiers[kind](proof)
            used.add(name)
            kinds[kind] = kinds.get(kind, 0)+1
        if 'representative_proof' in row:
            name = row['representative_proof']
            proof = proofs[name]
            assert proof['graph'] == graph_data(row)
            verify_graph(proof)
            used.add(name)
    assert used == set(proofs)
    assert len(results) == summary['reviewed'] == 1879
    assert sum(kinds.values()) == summary['certified'] == 1843
    assert kinds == summary['proof_kind_counts']
    assert len(results)-sum(kinds.values()) == summary['inconclusive_quotients'] == 36
    pool = read_rows(HERE/summary['published_pool_archive'])
    pool_keys = set()
    bindings = {}
    for row in pool:
        assert graph_quotient(row['graph6']) == row['quotient']
        key = canonical(row['quotient'])
        assert key not in pool_keys
        pool_keys.add(key)
        bindings[row['source'], row['index']] = row
    assert len(pool_keys) == summary['published_quotient_classes'] == 5589
    assert pool_keys == old | seen
    assert summary['published_pool_exactly_covered_by_reviews'] is True
    assert summary['combined_certified_before_new_phase_proofs'] == 5371
    assert summary['combined_residual_before_new_phase_proofs'] == 218
    guide = summary['best_surviving_guide']
    if guide is not None and 'individual_certificate' in guide:
        raw = gzip.decompress((HERE/guide['individual_certificate']).read_bytes())
        assert hashlib.sha256(raw).hexdigest() == guide['dominance_certificate_sha256']
        proof = json.loads(raw)
        assert str(verify_dual(proof)) == guide['exact_theta_upper']
        assert proof['edges'] == graph_data(guide)['edges']
    if inputs:
        scored, actual_keys, found, cache = 0, set(), set(), {}
        for relative in summary['source_files']:
            raw = subprocess.check_output(['git', 'show', summary['published_commit']+':'+relative], cwd=data_repo)
            assert hashlib.sha256(raw).hexdigest() == summary['source_sha256'][relative]
            text = gzip.decompress(raw).decode() if relative.endswith('.gz') else raw.decode()
            for line in text.splitlines():
                record = json.loads(line)
                if record.get('status') != 'scored_numerically':
                    continue
                scored += 1
                quotient = graph_quotient(record['graph6'])
                raw_key = tuple(sum(quotient, []))
                if raw_key not in cache:
                    cache[raw_key] = canonical(quotient)
                actual_keys.add(cache[raw_key])
                binding = bindings.get((relative, record['index']))
                if binding is not None:
                    assert binding['graph6'] == record['graph6']
                    assert binding['theta_numeric'] == record['theta_numeric']
                    found.add((relative, record['index']))
        assert scored == summary['scored_records'] == 34721
        assert actual_keys == pool_keys and found == set(bindings)
        print('Published source hashes, all 34,721 scored records, and complete 5,589-class census verified.')
    print('Verified 1,843 additional exact exclusions; all 5,589 published quotient classes have now been reviewed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', action='store_true')
    parser.add_argument('--data-repo', type=Path, default=Path('/workspace/ghz30-overnight'))
    args = parser.parse_args()
    main(args.inputs, args.data_repo)

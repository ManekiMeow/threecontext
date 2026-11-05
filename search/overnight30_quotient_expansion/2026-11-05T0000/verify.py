"""Exact replay of the bounded expanded published-quotient review."""
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import sys
import time

HERE = Path(__file__).resolve().parent
SEARCH = HERE.parents[1]
sys.path.insert(0, str(SEARCH))
from event30_disjoint_probability import verify_quotient_certificate as verify_disjoint
from event30_probability_lp import verify_certificate as verify_probability
from event30_quotient_saturation import verify_certificate as verify_saturation
from overnight30_fast_dual import verify_certificate as verify_dual
from overnight30_quotient import graph_data
from overnight30_screen import verify_certificate as verify_graph


def load_catalog():
    summary = json.loads((HERE/'summary.json').read_text())
    with gzip.open(HERE/summary['results_archive'], 'rt') as stream:
        results = [json.loads(line) for line in stream]
    with gzip.open(HERE/summary['proof_archive'], 'rt') as stream:
        records = [json.loads(line) for line in stream]
    certificates = {record['name']: record['certificate'] for record in records}
    assert len(certificates) == len(records) == summary['proof_count']
    return summary, results, certificates


def main():
    started = time.monotonic()
    hashes = json.loads((HERE/'artifact_sha256.json').read_text())
    for name, digest in hashes.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest
    summary, results, certificates = load_catalog()
    prior_dir = SEARCH/'overnight30_quotient_reviews'/'2026-11-04T2300'
    helpers = runpy.run_path(str(prior_dir/'verify.py'))
    graph_quotient, canonical = helpers['graph_quotient'], helpers['canonical']
    old = {canonical(row['quotient']) for row in json.loads((prior_dir/'summary.json').read_text())['results']}
    seen, kinds, used = set(), {}, set()
    verify = {'probability': verify_probability, 'saturation': verify_saturation,
              'disjoint_probability': verify_disjoint}
    for row in results:
        assert graph_quotient(row['graph6']) == row['quotient']
        key = canonical(row['quotient'])
        assert key not in old and key not in seen
        seen.add(key)
        if row['exact_quotient_excluded']:
            name, kind = row['certificate'], row['proof_kind']
            certificate = certificates[name]
            assert certificate['quotient'] == row['quotient']
            verify[kind](certificate)
            used.add(name)
            kinds[kind] = kinds.get(kind, 0)+1
        if 'representative_proof' in row:
            name = row['representative_proof']
            certificate = certificates[name]
            assert certificate['graph'] == graph_data(row)
            verify_graph(certificate)
            used.add(name)
    assert used == set(certificates)
    assert len(results) == summary['reviewed'] == 3560
    assert sum(kinds.values()) == summary['certified'] == 3378
    assert kinds == summary['proof_kind_counts']
    assert len(results)-sum(kinds.values()) == summary['inconclusive_quotients'] == 182
    guide = summary['best_surviving_guide']
    assert not guide['exact_quotient_excluded'] and canonical(guide['quotient']) in seen
    raw = gzip.decompress((HERE/guide['individual_certificate']).read_bytes())
    assert hashlib.sha256(raw).hexdigest() == guide['dominance_certificate_sha256']
    individual = json.loads(raw)
    assert str(verify_dual(individual)) == guide['exact_theta_upper']
    assert individual['edges'] == graph_data(guide)['edges']
    print('Verified 3378 further exact quotient exclusions, 3560 actual C5 representatives, and guide individual dual.')
    print('Exact replay seconds:', round(time.monotonic()-started, 3))


if __name__ == '__main__':
    main()

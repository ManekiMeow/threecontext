"""Verify five exact saturation proofs and the actual C5 quotient fixtures."""
import hashlib
import json
from pathlib import Path
import runpy
import sys

HERE = Path(__file__).resolve().parent
SEARCH = HERE.parents[1]
sys.path.insert(0, str(SEARCH))
from event30_c5_quotient import Q
from event30_quotient_saturation import verify_certificate


def main():
    hashes = json.loads((HERE / 'artifact_sha256.json').read_text())
    for name, digest in hashes.items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
    summary = json.loads((HERE / 'summary.json').read_text())
    prior_path = SEARCH.parent / summary['source_summary']
    assert hashlib.sha256(prior_path.read_bytes()).hexdigest() == summary['source_summary_sha256']
    prior = json.loads(prior_path.read_text())
    helpers = runpy.run_path(str(prior_path.parent / 'verify.py'))
    graph_quotient, canonical = helpers['graph_quotient'], helpers['canonical']
    old = {canonical(item['quotient']) for item in prior['results'] if item['exact_lp_excluded']}
    new = set()
    for item in summary['results']:
        certificate = json.loads((HERE / item['certificate']).read_text())
        assert certificate['quotient'] == item['quotient'] == graph_quotient(item['graph6'])
        assert verify_certificate(certificate) == item['terminal']
        source = prior['results'][item['source_index']]
        assert source['quotient'] == item['quotient'] and source['graph6'] == item['graph6']
        assert not source['exact_lp_excluded']
        key = canonical(item['quotient'])
        assert key not in new
        new.add(key)
        assert item['already_covered_by_proposition46'] == (key == canonical(Q))
    assert len(new) == summary['certified'] == summary['new_beyond_lp'] == 5
    assert len(new-{canonical(Q)}) == summary['new_beyond_lp_and_proposition46'] == 4
    assert len(old | new) == summary['combined_certified'] == 131
    assert prior['reviewed']-len(old | new) == summary['combined_unresolved'] == 19
    print('Verified five exact saturation certificates and their actual C5 quotient fixtures.')
    print('Combined with the prior LP review: 131 of150 selected quotient classes excluded;19 unresolved.')


if __name__ == '__main__':
    main()

"""Verify exact strict-LP controls and the scope of the inconclusive trial."""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use -O.')

import gzip
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
SEARCH = HERE.parents[1]
ROOT = SEARCH.parent
sys.path.insert(0, str(SEARCH))
from event30_strict_probability_lp import self_test


def read(path):
    with gzip.open(path, 'rt') as stream:
        return [json.loads(line) for line in stream]


def main():
    for name, digest in json.loads((HERE/'artifact_sha256.json').read_text()).items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest
    summary = json.loads((HERE/'summary.json').read_text())
    results = read(HERE/summary['results_archive'])
    proofs = read(HERE/summary['proof_archive'])
    assert not proofs and summary['certified'] == summary['zero_margin'] == 0
    assert len(results) == summary['reviewed'] == 218
    expected = {}
    for date in ('2026-11-05T0000', '2026-11-05T0100'):
        relative = Path('search/overnight30_quotient_expansion')/date
        rows = read(ROOT/relative/'results.jsonl.gz')
        for i, row in enumerate(rows):
            if not row['exact_quotient_excluded']:
                expected[str(relative), i] = row
    seen = set()
    for row in results:
        key = row['source_catalog'], row['source_result_index']
        assert key in expected and key not in seen
        seen.add(key)
        assert row['quotient'] == expected[key]['quotient']
        assert row['graph6'] == expected[key]['graph6']
        assert not row['exact_strict_excluded']
    assert seen == set(expected)
    report = self_test()
    print('Verified 218-class trial scope; no new exclusion or feasibility claim.')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

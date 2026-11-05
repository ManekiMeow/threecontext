"""Exact replay of five further complete published C5 phase classes."""
import gzip
import hashlib
import json
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use python -O')

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
from event30_c5_general_phases import verify


def verify_targets():
    source = json.loads((HERE/'targets.json').read_text())
    path = ROOT/source['source_archive']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == source['source_archive_sha256']
    with gzip.open(path, 'rt') as stream:
        original = [json.loads(line) for line in stream]
    expected = [dict(row, expansion_index=i) for i, row in enumerate(original)
                if not row['exact_quotient_excluded']]
    assert source['results'] == expected and len(expected) == 182


def main():
    verify_targets()
    targets = json.loads((HERE/'targets.json').read_text())['results']
    summary = json.loads((HERE/'summary.json').read_text())
    normalized = representatives = 0
    seen = set()
    for row in summary['results']:
        relative = Path(row['manifest'])
        assert not relative.is_absolute() and '..' not in relative.parts
        path = ROOT/relative
        raw = path.read_bytes()
        assert hashlib.sha256(raw).hexdigest() == row['manifest_sha256']
        manifest = json.loads(raw)
        assert manifest['quotient'] == row['quotient']
        assert manifest['source_index'] == row['source_index'] not in seen
        assert row['source_expansion_index'] == targets[row['source_index']]['expansion_index']
        seen.add(row['source_index'])
        report = verify(path.parent)
        assert report['complete_exact_exclusion'] is True
        assert report['uniform_exact_upper'] == row['uniform_exact_upper']
        normalized += report['coverage']['normal_triangle_free_assignments']
        representatives += report['covered_representatives']
    assert len(seen) == summary['certified'] == 5
    assert normalized == summary['normalized_assignments'] == 9062
    assert representatives == summary['certified_representatives'] == 2269
    print('Verified five complete phase exclusions: 9062 normalized assignments, 2269 exact graph certificates.')


if __name__ == '__main__':
    main()

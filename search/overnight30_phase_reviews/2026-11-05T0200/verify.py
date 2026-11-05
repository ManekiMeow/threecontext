"""Exact replay of the 02:00 residual C5 phase exclusions and source selection."""
import gzip
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import runpy
import sys

if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use python -O')

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
from event30_c5_general_phases import verify


def verify_bound_comparison():
    report = json.loads((HERE/'upper_bound_comparison.json').read_text())
    records = []
    for relative, digest in report['source_manifest_sha256'].items():
        path = Path(relative)
        assert not path.is_absolute() and '..' not in path.parts
        raw = (ROOT/path).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == digest
        records += [dict(row, manifest=relative) for row in json.loads(raw)['records']
                    if row.get('theta_numeric') is not None]
    assert len(records) == report['saved_records'] == 6532
    lookup = {(row['manifest'], row['graph6']): row for row in records}
    for row in report['top']+[report['largest_all_exact_bound_record'], report['largest_0200_exact_bound_record']]:
        saved = lookup[row['manifest'], row['graph6']]
        assert saved['theta_numeric'] == row['theta_numeric']
        assert saved['exact_theta_upper'] == row['exact_theta_upper']
    assert report['top'][0]['theta_numeric'] == max(row['theta_numeric'] for row in records)
    assert Fraction(report['largest_all_exact_bound_record']['exact_theta_upper']) == max(
        Fraction(row['exact_theta_upper']) for row in records)
    assert Fraction(report['largest_0200_exact_bound_record']['exact_theta_upper']) == max(
        Fraction(row['exact_theta_upper']) for row in records if '2026-11-05T0200' in row['manifest'])


def verify_targets():
    targets = json.loads((HERE/'targets.json').read_text())
    canonical = runpy.run_path(str(ROOT/'search/overnight30_quotient_reviews/2026-11-04T2300/verify.py'))['canonical']
    for relative, digest in targets['source_sha256'].items():
        path = Path(relative)
        assert not path.is_absolute() and '..' not in path.parts
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    prior_phase = ROOT/'search/overnight30_phase_reviews/2026-11-05T0100/summary.json'
    excluded = {canonical(row['quotient']) for row in json.loads(prior_phase.read_text())['results']}
    expected = []
    for hour in ('0000', '0100'):
        path = ROOT/f'search/overnight30_quotient_expansion/2026-11-05T{hour}/results.jsonl.gz'
        with gzip.open(path, 'rt') as stream:
            for index, line in enumerate(stream):
                row = json.loads(line)
                if not row['exact_quotient_excluded'] and canonical(row['quotient']) not in excluded:
                    expected.append(dict(row, source_expansion_archive=str(path.relative_to(ROOT)),
                                         source_expansion_index=index))
    expected.sort(key=lambda row: (-row['theta_numeric'], row['source_expansion_archive'], row['source_expansion_index']))
    assert targets['results'] == expected and len(expected) == 213
    assert len({canonical(row['quotient']) for row in expected}) == 213
    closures = []
    for index, row in enumerate(expected):
        quotient, remaining = row['quotient'], set(range(6))
        while True:
            removed = {v for v in remaining if not any(
                sum(quotient[v][w] for w in range(2*k, 2*k+2)) <= 2
                and all(w in remaining for w in range(2*k, 2*k+2) if quotient[v][w])
                for k in range(3) if k != v//2)}
            if not removed:
                break
            remaining -= removed
        if any({2*k, 2*k+1} <= remaining for k in range(3)):
            closures.append({'index': index, 'closed_orbits': sorted(remaining)})
    prescreen = json.loads((HERE/'closure_prescreen.json').read_text())
    assert prescreen['reviewed'] == 213 and prescreen['excluded'] == len(closures) == 0
    assert prescreen['results'] == closures
    return expected


def main():
    targets = verify_targets()
    verify_bound_comparison()
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
        assert row['quotient'] == targets[row['source_index']]['quotient']
        seen.add(row['source_index'])
        report = verify(path.parent)
        assert report['complete_exact_exclusion'] is True
        assert report['uniform_exact_upper'] == row['uniform_exact_upper']
        normalized += report['coverage']['normal_triangle_free_assignments']
        representatives += report['covered_representatives']
    assert len(seen) == summary['certified']
    assert normalized == summary['normalized_assignments']
    assert representatives == summary['certified_representatives']
    assert summary['fixed_pool_classes'] == 5589
    assert summary['fixed_pool_previously_excluded'] == 5376
    assert summary['fixed_pool_now_excluded'] == 5376+len(seen)
    assert summary['fixed_pool_remaining'] == 213-len(seen)
    print(f'Verified {len(seen)} complete phase exclusions: {normalized} normalized assignments, {representatives} exact graph certificates.')


if __name__ == '__main__':
    main()

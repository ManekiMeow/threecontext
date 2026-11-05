"""Reuse the independently developed fixed-order residue-set phase audit.

The discovery path uses MRV bitmask enumeration. This audit instead uses the
frozen, separately written fixed-order explicit-set implementation from the
preceding general-phase review, then verifies every newly bound exact dual.
"""
import argparse
import hashlib
import json
from pathlib import Path
import runpy

if not __debug__:
    raise RuntimeError('Exact audit requires assertions; do not use python -O')

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
AUDIT = ROOT/'search/overnight30_general_phase_reviews/2026-11-05T0000/independent_coverage.py'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    audit = runpy.run_path(str(AUDIT))['audit']
    summary = json.loads((HERE/'summary.json').read_text())
    results = [audit((ROOT/row['manifest']).parent) for row in summary['results']]
    report = {'method': 'independent fixed-order residue sets, unit composition, translations, exact duals',
              'audit_implementation': str(AUDIT.relative_to(ROOT)),
              'audit_implementation_sha256': hashlib.sha256(AUDIT.read_bytes()).hexdigest(),
              'classes': results,
              'total_certificates': sum(row['exact_certificates_verified'] for row in results)}
    text = json.dumps(report, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)


if __name__ == '__main__':
    main()

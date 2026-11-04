"""Replay the saved probability-propagation exclusions using only stdlib."""
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from overnight30_screen import verify_certificate


def main():
    hashes = json.loads((HERE / 'artifact_sha256.json').read_text())
    for name, digest in hashes.items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
    summary = json.loads((HERE / 'summary.json').read_text())
    entries = summary['new_fast_exclusions']
    assert len(entries) == 16
    for entry in entries:
        certificate = json.loads((HERE / entry['certificate']).read_text())
        assert certificate['method'] == 'probability_propagation'
        assert len(certificate['witness']['rounds']) <= 3
        verify_certificate(certificate)
    print('Verified16 exact propagation certificates and artifact hashes.')


if __name__ == '__main__':
    main()

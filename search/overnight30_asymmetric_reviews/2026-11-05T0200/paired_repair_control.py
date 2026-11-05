"""Nonnumerical matched-seed comparison of single- and paired-block repairs.

The original worker and screen are loaded from the recorded source commit.
No SDP, SAT discovery, or optimizer is invoked. Screen survival is inconclusive.
"""
import argparse
from collections import Counter
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
CODE = 'fd05e9b3d5b16f806795db0e5d98ad1e6acb1895'
DATA = '3536304a360437b6c54413602377c31438c719f0'
PARENT = 'search/overnight30_runs/2026-11-04/chunks/00098_asymmetric/asymmetric_seed20261104_best_excluded.json'


def paired_repair(record, rng, w):
    """Repair two disjoint nonsingleton witness blocks before maximal filling."""
    if not record or (record.get('certificate') or {}).get('method') != 'disjoint_partition':
        return None
    cert = record['certificate']
    w['verify_certificate'](cert)
    parent = w['decode'](record['g6'])
    if (cert['graph']['parts'] != w['PARTS']
            or w['from_edges'](cert['graph']['edges']) != parent
            or not w['validate'](parent)):
        raise ValueError('partition repair certificate is not bound to its guide')
    blocks = [block for block in cert['witness']['blocks'] if len(block['vertices']) >= 2]
    if len(blocks) < 2:
        return None
    for _ in range(8):
        selected = rng.sample(blocks, 2)
        targets = []
        for block in selected:
            vertices, other = block['vertices'], block['other_context']
            choices = [(u, v, b) for u in vertices for v in vertices if u != v
                       for b in w['PARTS'][other]
                       if parent[v] >> b & 1 and not parent[u] >> b & 1]
            if not choices:
                break
            targets.append((*rng.choice(choices), block))
        if len(targets) != 2:
            continue
        protected = {tuple(sorted(edge)) for u, v, b, _ in targets for edge in ((u, b), (v, b))}
        trial, repairs = parent[:], []
        for u, v, b, block in targets:
            cleared, blocked = [], False
            for c in range(w['N']):
                if (trial[u] & trial[b]) >> c & 1:
                    choices = [edge for edge in ((u, c), (b, c)) if tuple(sorted(edge)) not in protected]
                    if not choices:
                        blocked = True
                        break
                    edge = rng.choice(choices)
                    w['remove'](trial, *edge)
                    cleared.append(list(edge))
            if blocked or not w['addable'](trial, u, b):
                break
            w['add'](trial, u, b)
            repairs.append({'certificate_block': block, 'added_edge': [u, b],
                            'sharing_pair': [u, v], 'shared_neighbor': b,
                            'removed_conflicts': cleared})
        if len(repairs) != 2:
            continue
        w['fill'](trial, rng)
        if not all(trial[u] >> v & 1 for u, v in protected) or not w['validate'](trial):
            continue
        return trial, {'parent_g6': record['g6'], 'certificate_method': 'disjoint_partition',
                       'certificate_context': cert['witness']['context'], 'repairs': repairs,
                       'edge_difference': len(set(w['edges_of'](parent)) ^ set(w['edges_of'](trial)))}
    return None


def frozen_environment(repo):
    bindings = {}
    namespaces = {}
    for name in ('overnight30_screen.py', 'overnight30_asymmetric.py'):
        path = 'search/'+name
        raw = subprocess.check_output(['git', 'show', CODE+':'+path], cwd=ROOT)
        bindings[path] = hashlib.sha256(raw).hexdigest()
        namespace = {'__name__': 'frozen_'+name[:-3], '__file__': str(ROOT/path)}
        exec(compile(raw, str(ROOT/path), 'exec'), namespace)
        namespaces[name] = namespace
    w = namespaces['overnight30_asymmetric.py']
    w['fast_screen'] = namespaces['overnight30_screen.py']['fast_screen']
    raw = subprocess.check_output(['git', 'show', DATA+':'+PARENT], cwd=repo)
    bindings[PARENT] = hashlib.sha256(raw).hexdigest()
    return w, json.loads(raw), bindings


def run(repo, count):
    w, parent, bindings = frozen_environment(repo)
    before = deepcopy(parent)
    records, statistics = [], {}
    for mode in ('single', 'paired'):
        started = time.monotonic()
        counters = Counter()
        labels = set()
        for seed in range(count):
            rng = random.Random(202611050200+seed)
            result = w['partition_witness_repair'](parent, rng) if mode == 'single' else paired_repair(parent, rng, w)
            if result is None:
                counters['fallback'] += 1
                continue
            trial, provenance = result
            assert w['validate'](trial) and parent == before
            again_rng = random.Random(202611050200+seed)
            again = w['partition_witness_repair'](parent, again_rng) if mode == 'single' else paired_repair(parent, again_rng, w)
            assert again == result and rng.getstate() == again_rng.getstate()
            code = w['encode'](trial)
            labels.add(code)
            screened = w['fast_screen']({'g6': code})
            cert = screened.get('certificate')
            if cert:
                w['verify_certificate'](cert)
                counters[cert['method']] += 1
            else:
                counters['not_excluded'] += 1
            counters['valid'] += 1
            records.append({'mode': mode, 'seed': 202611050200+seed, 'g6': code,
                            'provenance': provenance, 'screen_status': screened['status'],
                            'certificate': cert})
        statistics[mode] = {'counts': dict(counters), 'distinct_labels': len(labels),
                            'elapsed_seconds': time.monotonic()-started}
    return {'code_commit': CODE, 'data_commit': DATA, 'parent': PARENT,
            'matched_seeds_per_mode': count, 'inputs_sha256': bindings,
            'statistics': statistics,
            'scope': 'Combinatorial proposal comparison only; no numerical scoring, feasibility or graph-isomorphism claim.'}, records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-repo', type=Path, default=Path('/workspace/ghz30-overnight'))
    parser.add_argument('--seeds', type=int, default=128)
    parser.add_argument('--output-dir', type=Path, default=Path('/tmp/ghz30-paired0200-replay'))
    args = parser.parse_args()
    report, records = run(args.data_repo, args.seeds)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/'paired_control.json').write_text(json.dumps(report, indent=2)+'\n')
    raw = ''.join(json.dumps(row, separators=(',', ':'))+'\n' for row in records).encode()
    (args.output_dir/'paired_control_records.jsonl.gz').write_bytes(gzip.compress(raw, mtime=0))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

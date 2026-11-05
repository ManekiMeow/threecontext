"""Reproduce the pinned natural asymmetric-move performance audit.

This reads immutable Git blobs. It launches no searches or numerical solvers.
Novel labels are compared with earlier completed asymmetric main-worker pools;
neither global main-pool novelty nor graph isomorphism is inferred.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
COMMIT = '3536304a360437b6c54413602377c31438c719f0'
RUN = 'search/overnight30_runs/2026-11-04'
REVIEW_JOBS = (82, 85, 88, 92, 95, 98)
REPAIR_JOBS = (88, 92, 95, 98)


def analyze(repo):
    inputs = {}

    def read(path):
        raw = subprocess.check_output(['git', 'show', f'{COMMIT}:{path}'], cwd=repo)
        inputs[path] = hashlib.sha256(raw).hexdigest()
        return raw

    manifest = json.loads(read(RUN+'/manifest.json'))
    jobs = [job for job in manifest['jobs']
            if job['kind'] == 'asymmetric' and job['state'] == 'completed']
    seen_scored, seen_repair_labels = set(), set()
    previous = None
    reports, all_repair_rows = [], []
    for job in jobs:
        prefix = job['directory']+'/asymmetric_seed20261104'
        state = json.loads(read(prefix+'.json'))
        rows = [json.loads(line) for line in gzip.decompress(read(prefix+'.jsonl.gz')).splitlines()]
        records = [row for row in rows if 'theta' in row]
        fresh = Counter()
        for row in records:
            fresh[row.get('move', 'missing')] += row['g6'] not in seen_scored
            seen_scored.add(row['g6'])
        if job['id'] not in REVIEW_JOBS:
            previous = state
            continue
        assert previous is not None
        by_move = defaultdict(list)
        for row in records:
            by_move[row.get('move', 'missing')].append(row)
        move_rows = []
        for move in sorted(state['moves']):
            proposed = state['moves'][move]-previous['moves'].get(move, 0)
            if not proposed:
                continue
            scored = by_move[move]
            surviving = [row for row in scored if row['screen']['status'] != 'excluded']
            best = max(surviving, key=lambda row: row['theta'], default=None)
            move_rows.append({
                'move': move, 'proposed': proposed,
                'invalid': state.get('move_invalid', {}).get(move, 0)-previous.get('move_invalid', {}).get(move, 0),
                'accepted': state.get('move_accepted', {}).get(move, 0)-previous.get('move_accepted', {}).get(move, 0),
                'scored_records': len(scored),
                'new_scored_labels_against_asymmetric_history': fresh[move],
                'configured_screen_surviving_scores': len(surviving),
                'best_surviving_score': best['theta'] if best else None,
                'best_surviving_g6': best['g6'] if best else None,
                'scores_during_excursions': sum(row.get('excursion', False) for row in scored),
            })
        fields = ('iterations', 'scored', 'invalid', 'duplicate', 'accepted',
                  'fast_rejected', 'sat_excluded', 'solver_failures', 'elapsed_seconds')
        delta = {key: state[key]-previous[key] for key in fields}
        assert delta['scored'] == len(records)
        assert delta['iterations'] == sum(row['proposed'] for row in move_rows)
        excursion_fields = ('started', 'completed', 'evaluated', 'newly_scored', 'accepted', 'downhill_accepted')
        excursions = {key: state['excursions'][key]-previous['excursions'][key] for key in excursion_fields}
        assert excursions['newly_scored'] == sum(row.get('excursion', False) for row in records)
        repair_rows = [row for row in rows if 'partition_repair' in row]
        if job['id'] in REPAIR_JOBS:
            all_repair_rows.extend(repair_rows)
        before = previous.get('partition_repair', {})
        repair_counters = {key: state.get('partition_repair', {}).get(key, 0)-before.get(key, 0)
                           for key in ('ordinary_proposals', 'attempted', 'proposed', 'fallbacks')}
        best = state['best_records']['best_theta']
        provenance = next((row for row in records if row['g6'] == best['g6']), None)
        reports.append({
            'job': job['id'], 'source': prefix, 'started_utc': job['started_utc'],
            'finished_utc': job['finished_utc'], 'code_commit': job['code_commit'],
            'source_sha256': job['source_sha256'],
            'record_count': len(rows), 'counter_deltas': delta, 'moves': move_rows,
            'excursion_deltas': excursions, 'repair_deltas': repair_counters,
            'repair_logged_records': len(repair_rows),
            'repair_logged_statuses': dict(Counter(row['screen']['status'] for row in repair_rows)),
            'active_best_theta': state['best_theta'], 'active_best_g6': best['g6'],
            'active_best_origin_move_in_this_chunk': provenance.get('move') if provenance else None,
            'historical_best_theta': state['historical_best_theta'],
            'shared_seeds_added': state.get('shared_seed_last_added'),
        })
        previous = state
    assert [row['job'] for row in reports] == list(REVIEW_JOBS)
    repair = {key: sum(row['repair_deltas'][key] for row in reports)
              for key in ('ordinary_proposals', 'attempted', 'proposed', 'fallbacks')}
    repair['logged_records'] = len(all_repair_rows)
    repair['distinct_logged_labels'] = len({row['g6'] for row in all_repair_rows})
    repair['logged_screen_statuses'] = dict(Counter(row['screen']['status'] for row in all_repair_rows))
    repair_moves = [move for row in reports for move in row['moves'] if move['move'] == 'partition_witness_repair']
    for field in ('accepted', 'invalid', 'scored_records', 'new_scored_labels_against_asymmetric_history',
                  'configured_screen_surviving_scores'):
        repair[field] = sum(row[field] for row in repair_moves)
    repair['highest_direct_score'] = max(row['best_surviving_score'] for row in repair_moves)
    summary = {
        'pinned_data_commit': COMMIT, 'scope': 'Natural completed-worker observations, not a randomized move comparison.',
        'label_novelty_scope': 'Exact labels against earlier completed asymmetric main-worker scored pools only; no isomorphism or global-main-pool claim.',
        'jobs': reports, 'repair_totals': repair,
        'cost_limitation': 'Per-move solver times and complete parent ancestry were not logged. Worker elapsed-time deltas cannot isolate a move cost or causal contribution.',
        'recommendation': 'Keep the 128-proposal cadence and excursion settings; direct repair scores have not exceeded the injected guide. The paired-block variant is assessed separately by the matched combinatorial control.',
        'inputs_sha256': inputs,
    }
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-repo', type=Path, default=Path('/workspace/ghz30-overnight'))
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = analyze(args.data_repo)
    destination = HERE/'natural_performance.json'
    if args.write:
        destination.write_text(json.dumps(result, indent=2)+'\n')
    else:
        assert json.loads(destination.read_text()) == result
    print(json.dumps(result['repair_totals'], indent=2))
    print('Verified pinned counters, move attribution and asymmetric-history label novelty.')


if __name__ == '__main__':
    main()

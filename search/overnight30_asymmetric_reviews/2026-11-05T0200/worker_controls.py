"""Nonnumerical controls for the paired repair and persisted selector."""
import ast
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import random
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
import overnight30_asymmetric as worker
from event30_disjoint_partition import screen as partition_screen
from paired_repair_control import frozen_environment, paired_repair


def run():
    frozen, parent, bindings = frozen_environment(Path('/workspace/ghz30-overnight'))
    singleton_path = ROOT/'search/overnight30_reviews/2026-11-05T0100/4c1f39a97ef002e02c6f_certificate.json'
    proof = json.loads(singleton_path.read_text())
    singleton = {'g6': worker.encode(worker.from_edges(proof['graph']['edges'])), 'certificate': proof}
    original_inputs = deepcopy((parent, singleton))
    controls = 0
    for record in (parent, singleton):
        for seed in range(32):
            rng = random.Random(202611050210+seed)
            before = json.loads(json.dumps(rng.getstate()))
            result = worker.paired_partition_witness_repair(record, rng)
            assert result is not None
            trial, provenance = result
            assert worker.validate(trial) and provenance['repair_mode'] == 'paired'
            assert len(provenance['repairs']) == 2
            first, second = provenance['repairs']
            assert set(first['certificate_block']['vertices']).isdisjoint(second['certificate_block']['vertices'])
            for repair in provenance['repairs']:
                u, v = repair['sharing_pair']; b = repair['shared_neighbor']
                assert trial[u] >> b & 1 and trial[v] >> b & 1
            replay_rng = random.Random()
            replay_rng.setstate(worker.tuple_tree(before))
            assert worker.paired_partition_witness_repair(record, replay_rng) == result
            assert replay_rng.getstate() == rng.getstate()
            # The approved standalone prototype and production helper agree.
            prototype_rng = random.Random(202611050210+seed)
            prototype = paired_repair(record, prototype_rng, frozen)
            assert prototype[0] == trial
            expected = dict(provenance); expected.pop('repair_mode')
            assert prototype[1] == expected and prototype_rng.getstate() == rng.getstate()
            controls += 1
    assert (parent, singleton) == original_inputs
    for seed in range(8):
        assert worker.partition_witness_repair(parent, random.Random(seed)) == frozen['partition_witness_repair'](parent, random.Random(seed))
    bad = deepcopy(parent); bad['g6'] = singleton['g6']
    try:
        worker.paired_partition_witness_repair(bad, random.Random(1))
    except ValueError:
        pass
    else:
        raise AssertionError('mismatched certificate accepted')

    cycle = [v for i in range(10) for v in (i, 10+i, 20+i)]
    edges = [[cycle[i], cycle[(i+1) % 30]] for i in range(30)]
    graph = {'n': 30, 'parts': worker.PARTS, 'edges': [list(e) for e in worker.edges_of(worker.from_edges(edges))]}
    single_block_proof = partition_screen(graph)['certificate']
    assert len(single_block_proof['witness']['blocks']) == 1
    fallback_parent = {'g6': worker.encode(worker.from_edges(edges)), 'certificate': single_block_proof}
    assert worker.paired_partition_witness_repair(fallback_parent, random.Random(1)) is None

    source = (ROOT/'search/overnight30_asymmetric.py').read_text()
    tree = ast.parse(source)
    defaults = [node for node in ast.walk(tree) if isinstance(node, ast.For)
                and 'paired_attempted' in ast.unparse(node.iter)]
    selectors = [node for node in ast.walk(tree) if isinstance(node, ast.If)
                 and ast.unparse(node.test) == 'not exploring'
                 and 'repair_state' in ast.unparse(node)]
    assert len(defaults) == len(selectors) == 1
    compiled_defaults = compile(ast.fix_missing_locations(ast.Module(body=defaults, type_ignores=[])), '<actual counter defaults>', 'exec')
    compiled_selector = compile(ast.fix_missing_locations(ast.Module(body=selectors, type_ignores=[])), '<actual proposal selector>', 'exec')
    selector_checks = 0
    for attempted, expected_mode in ((258, 'single'), (259, 'paired')):
        old = {'partition_repair': {'ordinary_proposals': 127, 'attempted': attempted,
                                   'proposed': attempted, 'fallbacks': 0},
               'best_records': {'best_excluded_theta': parent}}
        namespace = dict(worker.__dict__, state=deepcopy(old), exploring=False,
                         repaired=None, rng=random.Random(24))
        exec(compiled_defaults, namespace)
        initial_state = json.loads(json.dumps(namespace['state']))
        initial_rng = json.loads(json.dumps(namespace['rng'].getstate()))
        exec(compiled_selector, namespace)
        assert namespace['state']['partition_repair']['ordinary_proposals'] == 128
        assert namespace['state']['partition_repair']['attempted'] == attempted+1
        assert namespace['repaired'] is not None
        assert namespace['repaired'][1].get('repair_mode', 'single') == expected_mode
        replay = dict(worker.__dict__, state=initial_state, exploring=False, repaired=None, rng=random.Random())
        replay['rng'].setstate(worker.tuple_tree(initial_rng))
        exec(compiled_selector, replay)
        assert namespace['state'] == replay['state'] and namespace['repaired'] == replay['repaired']
        assert namespace['rng'].getstate() == replay['rng'].getstate()
        selector_checks += 1
    namespace = dict(worker.__dict__, state={'partition_repair': {'ordinary_proposals': 127, 'attempted': 1,
                      'proposed': 1, 'fallbacks': 0}, 'best_records': {'best_excluded_theta': fallback_parent}},
                     exploring=False, repaired=None, rng=random.Random(3))
    exec(compiled_defaults, namespace); exec(compiled_selector, namespace)
    counters = namespace['state']['partition_repair']
    assert counters['paired_attempted'] == counters['paired_fallbacks'] == 1
    assert counters['paired_proposed'] == 0
    assert namespace['repaired'] is not None and 'repair_mode' not in namespace['repaired'][1]
    return {'paired_proposals': controls, 'checkpoint_rng_replays': controls,
            'prototype_matches': controls, 'unchanged_single_controls': 8,
            'real_singleton_partition_control': [5, 4, 1],
            'input_records_unmodified': True, 'mismatched_certificate_rejected': True,
            'actual_selector_json_replays': selector_checks,
            'actual_paired_to_single_fallback': True,
            'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
            'scope': 'Exact combinatorial and checkpoint controls; no numerical solves.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='Explicitly regenerate the saved control report.')
    args = parser.parse_args()
    report = run()
    if args.write:
        (HERE/'worker_controls.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))

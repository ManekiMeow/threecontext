"""Exact context probability budgets from disjoint-neighborhood partitions.

Each block U of a context partition may choose either other context, provided
the neighborhoods of its vertices there are pairwise disjoint. A verified
moment certificate gives c(U)<=|U|*u_|U|; positive definiteness of its trivial
dual makes this strict. Singleton blocks use the strict30-event bound1/9.
Summing a partition excludes a context budget<1, or a budget1 with at least
one strict block. No symmetry or uniform probabilities are assumed.

Discovery minimizes the budget over all partitions of each10-event context.
All verified catalog costs have common denominator1440, so the1024-mask DP
uses integers. Eligibility is computed by low-bit removal and neighborhood
unions. Catalog certificates are loaded and checked once per process. Saved
exclusion witnesses embed the exact bounds they use and replay independently
of the discovery catalog or dynamic program.
"""
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path


SCALE = 1440
FOLDER = Path(__file__).with_name('event30_disjoint_partition_certificates')


@lru_cache(maxsize=1)
def catalog():
    from event30_disjoint_probability import verify_bound_certificate, bound_is_strict
    costs, proofs = {1: 160}, {}
    for m in range(2, 11):
        proof = json.loads((FOLDER/f'mean_bound_{m:02}.json').read_text())
        assert proof['subset_size'] == m
        total = m*verify_bound_certificate(proof)
        assert bound_is_strict(proof)
        scaled = SCALE*total
        assert scaled.denominator == 1
        costs[m], proofs[m] = scaled.numerator, proof
    return costs, proofs


def optimal_partition(graph, context, adj, costs):
    part = graph['parts'][context]
    eligible = {1 << i: None for i in range(10)}
    for other in range(3):
        if other == context:
            continue
        positions = {v: i for i, v in enumerate(graph['parts'][other])}
        masks = [sum(1 << positions[v] for v in adj[u] if v in positions) for u in part]
        unions, disjoint = [0]*1024, [False]*1024
        disjoint[0] = True
        for mask in range(1, 1024):
            low = mask & -mask
            i, rest = low.bit_length()-1, mask ^ low
            unions[mask] = unions[rest] | masks[i]
            disjoint[mask] = disjoint[rest] and not (unions[rest] & masks[i])
            if disjoint[mask] and mask.bit_count() >= 2 and mask not in eligible:
                eligible[mask] = other
    by_vertex = [[] for _ in range(10)]
    block_cost = {}
    for subset in eligible:
        block_cost[subset] = costs[subset.bit_count()]
        for i in range(10):
            if subset & (1 << i):
                by_vertex[i].append(subset)
    dp, choice = [10*SCALE]*1024, [0]*1024
    dp[0] = 0
    for mask in range(1, 1024):
        i = (mask & -mask).bit_length()-1
        for subset in by_vertex[i]:
            if subset & mask == subset:
                value = block_cost[subset]+dp[mask ^ subset]
                if value < dp[mask]:
                    dp[mask], choice[mask] = value, subset
    blocks, mask = [], 1023
    while mask:
        subset = choice[mask]
        assert subset
        blocks.append({'vertices': [part[i] for i in range(10) if subset & (1 << i)],
                       'other_context': eligible[subset]})
        mask ^= subset
    return dp[1023], blocks


def verify_certificate(certificate):
    if not __debug__:
        raise RuntimeError('partition verification requires assertions; do not use python -O')
    from event30_disjoint_probability import verify_bound_certificate, bound_is_strict
    from overnight30_screen import normalize_graph, base_validation
    assert certificate['format'] == 'overnight30-certificate-v1'
    assert certificate['method'] == 'disjoint_partition'
    graph, witness = certificate['graph'], certificate['witness']
    assert normalize_graph(graph) == graph
    adj, _ = base_validation(graph)
    context = witness['context']
    assert type(context) is int and 0 <= context < 3
    bounds = {}
    for proof in witness['mean_bound_certificates']:
        m = proof['subset_size']
        assert type(m) is int and 2 <= m <= 10 and m not in bounds
        bounds[m] = (m*verify_bound_certificate(proof), bound_is_strict(proof))
    covered, used = set(), set()
    budget, strict = F(0), False
    for block in witness['blocks']:
        vertices, other = block['vertices'], block['other_context']
        assert vertices and all(type(v) is int for v in vertices)
        assert vertices == sorted(set(vertices))
        assert set(vertices) <= set(graph['parts'][context]) and not covered.intersection(vertices)
        covered.update(vertices)
        m = len(vertices)
        if m == 1:
            assert other is None
            budget += F(1, 9)
            strict = True
        else:
            assert type(other) is int and 0 <= other < 3 and other != context
            seen = set()
            for vertex in vertices:
                neighbors = adj[vertex] & set(graph['parts'][other])
                assert not seen.intersection(neighbors)
                seen.update(neighbors)
            total, block_strict = bounds[m]
            budget += total
            strict |= block_strict
            used.add(m)
    assert covered == set(graph['parts'][context]) and used == set(bounds)
    assert str(budget) == witness['budget']
    assert budget < 1 or (budget == 1 and strict)
    return True


def screen(record):
    from overnight30_screen import normalize_graph, base_validation
    graph = normalize_graph(record)
    adj, _ = base_validation(graph)
    costs, proofs = catalog()
    partitions = [optimal_partition(graph, k, adj, costs) for k in range(3)]
    result = {'status': 'not_excluded',
              'context_budgets': [str(F(value, SCALE)) for value, _ in partitions]}
    for context, (value, blocks) in enumerate(partitions):
        if value <= SCALE:
            sizes = sorted({len(block['vertices']) for block in blocks if len(block['vertices']) > 1})
            certificate = {'format': 'overnight30-certificate-v1', 'graph': graph,
                           'method': 'disjoint_partition',
                           'witness': {'context': context, 'blocks': blocks,
                                       'mean_bound_certificates': [deepcopy(proofs[m]) for m in sizes],
                                       'budget': str(F(value, SCALE))}}
            verify_certificate(certificate)
            result.update(status='excluded', certificate=certificate)
            break
    return result


def self_test():
    from overnight30_screen import normalize_graph, fast_screen, verify_certificate as verify_global
    order = [10*k+i for i in range(10) for k in range(3)]
    graph = normalize_graph({'parts': [list(range(10*k, 10*k+10)) for k in range(3)],
                             'edges': list(zip(order, order[1:]+order[:1]))})
    result = screen(graph)
    assert result['status'] == 'excluded' and result['context_budgets'] == ['27/40']*3
    proof = result['certificate']
    assert len(proof['witness']['blocks']) == 1
    _, bounds = catalog()
    equality = deepcopy(proof)
    equality['witness'] = {'context': 0, 'blocks': [
        {'vertices': list(range(5)), 'other_context': 1},
        {'vertices': list(range(5, 8)), 'other_context': 1},
        {'vertices': [8, 9], 'other_context': 1}],
        'mean_bound_certificates': [deepcopy(bounds[m]) for m in (2, 3, 5)], 'budget': '1'}
    assert verify_certificate(equality)
    invalid = []
    # Valid Bessel bound: c<=1/m, witnessed by a rank-one trivial dual.
    # Its total budget1 is NOT strict and must not exclude a whole context.
    nonstrict = deepcopy(proof)
    nonstrict['witness']['budget'] = '1'
    nonstrict['witness']['mean_bound_certificates'] = [{
        'format': 'ghz-disjoint-neighborhood-mean-bound-v1', 'subset_size': 10,
        'mean_probability_upper': '1/10', 'standard_dual': [['0']*3 for _ in range(3)],
        'trivial_dual': [['1/10', '-1', '0', '0'], ['-1', '10', '0', '0'],
                         ['0', '0', '0', '0'], ['0', '0', '0', '0']]}]
    from event30_disjoint_probability import verify_bound_certificate, bound_is_strict
    moment = nonstrict['witness']['mean_bound_certificates'][0]
    assert verify_bound_certificate(moment) == F(1, 10) and not bound_is_strict(moment)
    invalid.append(nonstrict)
    changed = deepcopy(proof); changed['witness']['budget'] = '1/2'; invalid.append(changed)
    changed = deepcopy(proof); changed['witness']['blocks'][0]['vertices'].pop(); invalid.append(changed)
    changed = deepcopy(proof); changed['witness']['blocks'][0]['other_context'] = 0; invalid.append(changed)
    changed = deepcopy(proof)
    changed['graph']['edges'] = sorted(changed['graph']['edges']+[[0, 11]])
    invalid.append(changed)
    changed = deepcopy(equality)
    changed['witness']['mean_bound_certificates'][0]['trivial_dual'][0][0] = '0'
    invalid.append(changed)
    for changed in invalid:
        try:
            verify_certificate(changed)
        except (AssertionError, KeyError):
            pass
        else:
            raise AssertionError('modified partition certificate accepted')
    fixtures = sorted(FOLDER.glob('control_*.json'))
    for path in fixtures:
        saved = json.loads(path.read_text())
        assert verify_global(saved)
        assert saved['witness']['budget'] == '15/16'
        assert fast_screen(saved['graph'])['certificate'] == saved
    return {'catalog_bounds_verified': 9, 'optimal_cycle_budget': '27/40',
            'strict_equality_partition': '5+3+2=1', 'tamper_controls_rejected': len(invalid),
            'valid_nonstrict_equality_rejected': True,
            'real_support_integration_controls': len(fixtures)}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('self-test', 'screen', 'verify'))
    parser.add_argument('path', type=Path, nargs='?')
    args = parser.parse_args()
    if args.command == 'self-test':
        output = self_test()
    elif args.command == 'screen':
        output = screen(json.loads(args.path.read_text()))
    else:
        output = {'verified': verify_certificate(json.loads(args.path.read_text()))}
    print(json.dumps(output, indent=2))

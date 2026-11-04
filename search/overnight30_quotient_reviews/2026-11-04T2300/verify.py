"""Verify saved quotient certificates and all 150 representative graphs exactly."""
from pathlib import Path
import itertools
import json
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
from event30_probability_lp import verify_certificate


def graph_quotient(encoded):
    data = encoded.encode('ascii')
    assert len(data) == 74 and data[0] == 93 and all(63 <= x <= 126 for x in data)
    adjacency = [set() for _ in range(30)]
    for v in range(1, 30):
        for u in range(v):
            bit = v*(v-1)//2+u
            if (data[1+bit//6]-63) & (1 << (5-bit % 6)):
                assert u//10 != v//10
                adjacency[u].add(v)
                adjacency[v].add(u)
    permutation = [5*(u//5)+(u+1) % 5 for u in range(30)]
    for u in range(30):
        assert {permutation[v] for v in adjacency[u]} == adjacency[permutation[u]]
        assert all(not adjacency[u].intersection(adjacency[v]) for v in adjacency[u])
    quotient = []
    for i in range(6):
        rows = [[sum(v//5 == j for v in adjacency[u]) for j in range(6)]
                for u in range(5*i, 5*i+5)]
        assert all(row == rows[0] for row in rows)
        quotient.append(rows[0])
    return quotient


def canonical(quotient):
    values = []
    for contexts in itertools.permutations(range(3)):
        for flips in itertools.product(range(2), repeat=3):
            p = [2*k+(a ^ flips[k]) for k in contexts for a in range(2)]
            values.append(tuple(quotient[p[i]][p[j]] for i in range(6) for j in range(6)))
    return min(values)


def main():
    summary = json.loads((ROOT/'summary.json').read_text())
    seen, certified, new = set(), 0, 0
    for item in summary['results']:
        quotient = graph_quotient(item['graph6'])
        assert quotient == item['quotient']
        key = canonical(quotient)
        assert key not in seen
        seen.add(key)
        if item['exact_lp_excluded']:
            certificate = json.loads((ROOT/item['certificate']).read_text())
            assert certificate['quotient'] == quotient
            assert certificate['sizes'] == [5]*6 and certificate['parts'] == [[0, 1], [2, 3], [4, 5]]
            assert str(verify_certificate(certificate)) == item['margin']
            certified += 1
            new += item['old_screen'] is None
    assert len(seen) == summary['reviewed'] == 150
    assert certified == summary['certified'] == 126
    assert new == summary['new_beyond_fast_screen'] == 108
    print(f'Verified {certified} exact quotient exclusions, with {len(seen)} distinct quotient representatives and their actual C5 automorphisms.')
    print('Uncertified quotient representatives remain unresolved by this probability test.')


if __name__ == '__main__':
    main()

"""Independent residue-set audit for proper-subset tree gauges.

This fixed-order enumeration does not import either phase discovery helper.
It checks normalized unit composition, all relative translations of one
representative per class, and every exact graph-bound certificate.
"""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError('Exact audit requires assertions; do not use -O')

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
from overnight30_fast_dual import verify_certificate

PAIRS = [(a, b) for a in range(6) for b in range(a+1, 6) if a//2 != b//2]
TRIPLES = list(itertools.product(range(2), range(2, 4), range(4, 6)))


def code(sets):
    return tuple(sum(1 << x for x in values) for values in sets)


def translated(values, shift):
    return frozenset((x+shift) % 5 for x in values)


def least_translate(values):
    translates = [translated(values, shift) for shift in range(5)]
    if 0 < len(values) < 5:
        assert len(set(translates)) == 5
    return min(translates, key=lambda value: sum(1 << x for x in value))


def normalize(sets, tree):
    edges = dict(zip(PAIRS, sets))
    shifts = {0: 0}
    while len(shifts) < 6:
        previous = len(shifts)
        for a, b in tree:
            target = least_translate(edges[a, b])
            possible = [shift for shift in range(5) if translated(edges[a, b], shift) == target]
            assert len(possible) == 1
            delta, = possible
            if a in shifts and b not in shifts:
                shifts[b] = (shifts[a]+delta) % 5
            elif b in shifts and a not in shifts:
                shifts[a] = (shifts[b]-delta) % 5
        assert len(shifts) > previous
    return tuple(translated(edges[a, b], shifts[b]-shifts[a]) for a, b in PAIRS)


def multiply(sets, unit, tree):
    return normalize(tuple(frozenset(unit*x % 5 for x in values) for values in sets), tree)


def enumerate_direct(quotient, tree):
    tree = [tuple(pair) for pair in tree]
    assert len(set(tree)) == len(tree) == 5
    assert all(pair in PAIRS and 0 < quotient[pair[0]][pair[1]] < 5 for pair in tree)
    connected = {0}
    for _ in range(5):
        connected |= {v for pair in tree if connected.intersection(pair) for v in pair}
    assert len(connected) == 6
    domains = {}
    for pair in PAIRS:
        values = set(map(frozenset, itertools.combinations(range(5), quotient[pair[0]][pair[1]])))
        if pair in tree:
            values = {least_translate(value) for value in values}
        domains[pair] = sorted(values, key=lambda value: sum(1 << x for x in value))
    assigned, normal, prefixes = {}, {}, 0

    def visit(depth):
        nonlocal prefixes
        prefixes += 1
        for a, b, c in TRIPLES:
            if all(pair in assigned for pair in ((a, b), (b, c), (a, c))):
                if any((x+y) % 5 in assigned[a, c] for x in assigned[a, b] for y in assigned[b, c]):
                    return
        if depth == 12:
            sets = tuple(assigned[pair] for pair in PAIRS)
            normal[code(sets)] = sets
            return
        pair = PAIRS[depth]
        for values in domains[pair]:
            assigned[pair] = values
            visit(depth+1)
        del assigned[pair]

    visit(0)
    representatives = {min(code(multiply(sets, unit, tree)) for unit in range(1, 5))
                       for sets in normal.values()}
    return normal, representatives, prefixes


def audit(folder):
    raw = (folder/'manifest.json').read_bytes()
    manifest = json.loads(raw)
    relative = Path(manifest['source_summary'])
    assert not relative.is_absolute() and '..' not in relative.parts
    source = (ROOT/relative).read_bytes()
    assert hashlib.sha256(source).hexdigest() == manifest['source_summary_sha256']
    assert json.loads(source)['results'][manifest['source_index']]['quotient'] == manifest['quotient']
    assert manifest['complete_exact_exclusion'] is True
    tree = [tuple(pair) for pair in manifest['coverage']['gauge_tree']]
    normal, representatives, prefixes = enumerate_direct(manifest['quotient'], tree)
    assert len(normal) == manifest['coverage']['normal_triangle_free_assignments']
    assert len(representatives) == manifest['coverage']['representatives']
    assert manifest['coverage']['raw_triangle_free_assignments'] == len(normal)*5**5
    recorded = {tuple(row['phase_masks']) for row in manifest['records']}
    assert len(recorded) == len(manifest['records']) and recorded == representatives
    compositions = 0
    for sets in normal.values():
        assert normalize(sets, tree) == sets
        for u, v in itertools.product(range(1, 5), repeat=2):
            assert multiply(multiply(sets, u, tree), v, tree) == multiply(sets, u*v % 5, tree)
            compositions += 1
    sizes = Counter(len({code(multiply(normal[rep], u, tree)) for u in range(1, 5)})
                    for rep in representatives)
    assert {str(k): v for k, v in sizes.items()} == manifest['coverage']['unit_orbit_sizes']
    assert sum(k*v for k, v in sizes.items()) == len(normal)
    initial = normal[min(representatives)]
    translations = set()
    for tail in itertools.product(range(5), repeat=5):
        shifts = (0,)+tail
        moved = tuple(translated(values, shifts[b]-shifts[a]) for (a, b), values in zip(PAIRS, initial))
        translations.add(code(moved))
        assert normalize(moved, tree) == initial
    assert len(translations) == 5**5
    bounds = []
    for row in manifest['records']:
        sets = normal[tuple(row['phase_masks'])]
        edges = sorted((5*a+i, 5*b+j) for (a, b), values in zip(PAIRS, sets)
                       for i in range(5) for j in range(5) if (j-i) % 5 in values)
        g6 = [ord(c)-63 for c in row['graph6']]
        assert len(g6) == 74 and g6[0] == 30 and all(0 <= x < 64 for x in g6)
        bits = [(value >> shift) & 1 for value in g6[1:] for shift in range(5, -1, -1)]
        possible = [(u, v) for v in range(1, 30) for u in range(v)]
        assert not any(bits[len(possible):])
        assert sorted(edge for edge, bit in zip(possible, bits) if bit) == edges
        relative = Path(row['certificate'])
        assert not relative.is_absolute() and '..' not in relative.parts
        data = gzip.decompress((folder/relative).read_bytes())
        assert hashlib.sha256(data).hexdigest() == row['certificate_sha256']
        certificate = json.loads(data)
        assert certificate['n'] == 30 and certificate['edges'] == [list(edge) for edge in edges]
        bound = verify_certificate(certificate)
        assert bound == Fraction(row['exact_theta_upper']) < 3
        bounds.append(bound)
    return {'source_index': manifest['source_index'], 'manifest_sha256': hashlib.sha256(raw).hexdigest(),
            'fixed_order_prefixes': prefixes, 'normal_triangle_free_assignments': len(normal),
            'unit_representatives': len(representatives), 'unit_compositions_checked': compositions,
            'relative_translations_checked': len(translations), 'exact_certificates_verified': len(bounds),
            'uniform_exact_upper': str(max(bounds))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records = [audit(HERE/f'index_{index:03}') for index in (3, 11, 103, 147)]
    report = {'method': 'independent fixed-order residue sets; normalized unit group action and exact duals',
              'classes': records, 'total_certificates': sum(r['exact_certificates_verified'] for r in records)}
    text = json.dumps(report, indent=2)+'\n'
    if args.output:
        args.output.write_text(text)
    print(text)

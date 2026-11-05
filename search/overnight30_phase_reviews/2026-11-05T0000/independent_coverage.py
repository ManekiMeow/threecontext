"""Independent fixed-order phase enumeration and exact certificate audit.

Uses explicit residue sets and direct triangle equations, independently of
event30_c5_phases' dynamic ordering, bit-mask addition tables, and gauge code.
The optional --output writes a replayable report; numerical tools are unused.
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
    raise RuntimeError("Exact audit requires Python assertions; do not use -O")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'search'))
from overnight30_fast_dual import verify_certificate

PAIRS = [(a, b) for a in range(6) for b in range(a+1, 6) if a//2 != b//2]
TRIPLES = list(itertools.product(range(2), range(2, 4), range(4, 6)))


def graph_edges(sets):
    return sorted((5*a+i, 5*b+j) for (a, b), values in zip(PAIRS, sets)
                  for i in range(5) for j in range(5) if (j-i) % 5 in values)


def decode_graph6(code):
    values = [ord(c)-63 for c in code]
    assert len(values) == 74 and values[0] == 30 and all(0 <= x < 64 for x in values)
    bits = [(value >> shift) & 1 for value in values[1:] for shift in range(5, -1, -1)]
    pairs = [(u, v) for v in range(1, 30) for u in range(v)]
    assert not any(bits[len(pairs):])
    return sorted(edge for edge, bit in zip(pairs, bits) if bit)


def masks(sets):
    return tuple(sum(1 << value for value in values) for values in sets)


def units(sets):
    return {masks([{u*x % 5 for x in values} for values in sets]) for u in range(1, 5)}


def enumerate_direct(quotient, tree):
    assert len(quotient) == 6 and all(len(row) == 6 for row in quotient)
    assert all(type(quotient[a][b]) is int and 0 <= quotient[a][b] <= 5
               and quotient[a][b] == quotient[b][a]
               and (a//2 != b//2 or quotient[a][b] == 0)
               for a in range(6) for b in range(6))
    tree = [tuple(pair) for pair in tree]
    assert len(tree) == len(set(tree)) == 5
    assert all(pair in PAIRS and quotient[pair[0]][pair[1]] == 1 for pair in tree)
    connected = {0}
    while True:
        expanded = connected | {v for edge in tree if connected.intersection(edge) for v in edge}
        if expanded == connected:
            break
        connected = expanded
    assert connected == set(range(6))  # Five connected edges on six vertices form a tree.
    assigned = {pair: frozenset([0]) for pair in tree}
    free = [pair for pair in PAIRS if pair not in assigned]
    domains = {pair: list(map(frozenset, itertools.combinations(range(5), quotient[pair[0]][pair[1]])))
               for pair in free}
    normal, prefixes = {}, 0

    def visit(depth):
        nonlocal prefixes
        prefixes += 1
        # Check all completed triangle equations directly in F5.
        for a, b, c in TRIPLES:
            if all(pair in assigned for pair in ((a, b), (b, c), (a, c))):
                if any((x+y) % 5 in assigned[a, c]
                       for x in assigned[a, b] for y in assigned[b, c]):
                    return
        if depth == len(free):
            sets = tuple(assigned[pair] for pair in PAIRS)
            normal[masks(sets)] = sets
            return
        pair = free[depth]
        for values in domains[pair]:
            assigned[pair] = values
            visit(depth+1)
        del assigned[pair]

    visit(0)
    representatives = {min(units(sets)) for sets in normal.values()}
    return normal, representatives, prefixes


def check_gauge(sets, tree):
    """Check all relative translations independently of the discovery normalizer."""
    orbit, expected = set(), masks(sets)
    for tail in itertools.product(range(5), repeat=5):
        shift = (0,)+tail
        moved = {pair: {(value+shift[pair[1]]-shift[pair[0]]) % 5 for value in values}
                 for pair, values in zip(PAIRS, sets)}
        orbit.add(masks([moved[pair] for pair in PAIRS]))
        correction = {0: 0}
        while len(correction) < 6:
            before = len(correction)
            for a, b in tree:
                value, = moved[a, b]
                if a in correction and b not in correction:
                    correction[b] = (correction[a]-value) % 5
                elif b in correction and a not in correction:
                    correction[a] = (correction[b]+value) % 5
            assert len(correction) > before
        restored = [{(value+correction[b]-correction[a]) % 5 for value in moved[a, b]}
                    for a, b in PAIRS]
        assert masks(restored) == expected
    assert len(orbit) == 5**5
    return len(orbit)


def audit(folder):
    raw = (folder/'manifest.json').read_bytes()
    manifest = json.loads(raw)
    source = ROOT/manifest['source_summary']
    assert not Path(manifest['source_summary']).is_absolute()
    assert hashlib.sha256(source.read_bytes()).hexdigest() == manifest['source_summary_sha256']
    assert json.loads(source.read_text())['results'][manifest['source_index']]['quotient'] == manifest['quotient']
    assert manifest['complete_exact_exclusion'] is True
    coverage = manifest['coverage']
    normal, representatives, prefixes = enumerate_direct(manifest['quotient'], coverage['gauge_tree'])
    recorded = {tuple(row['phase_masks']) for row in manifest['records']}
    assert len(recorded) == len(manifest['records']) and recorded == representatives
    assert len(normal) == coverage['normal_triangle_free_assignments']
    assert len(representatives) == coverage['representatives']
    assert coverage['raw_triangle_free_assignments'] == len(normal)*5**5
    orbit_sizes = Counter(len(units(normal[representative])) for representative in representatives)
    assert {str(size): count for size, count in orbit_sizes.items()} == coverage['unit_orbit_sizes']
    assert sum(size*count for size, count in orbit_sizes.items()) == len(normal)
    bounds = []
    for row in manifest['records']:
        edges = graph_edges(normal[tuple(row['phase_masks'])])
        assert decode_graph6(row['graph6']) == edges
        path = Path(row['certificate'])
        assert not path.is_absolute() and '..' not in path.parts
        certificate_raw = gzip.decompress((folder/path).read_bytes())
        assert hashlib.sha256(certificate_raw).hexdigest() == row['certificate_sha256']
        certificate = json.loads(certificate_raw)
        assert certificate['n'] == 30 and certificate['edges'] == [list(edge) for edge in edges]
        bound = verify_certificate(certificate)
        assert bound == Fraction(row['exact_theta_upper']) < 3
        bounds.append(bound)
    translations = check_gauge(normal[min(representatives)], coverage['gauge_tree'])
    return {'source_index': manifest['source_index'], 'manifest_sha256': hashlib.sha256(raw).hexdigest(),
            'fixed_order_prefixes': prefixes, 'normal_triangle_free_assignments': len(normal),
            'unit_representatives': len(representatives),
            'unit_orbit_sizes': {str(k): v for k, v in sorted(orbit_sizes.items())},
            'relative_translations_checked': translations, 'exact_certificates_verified': len(bounds),
            'uniform_exact_upper': str(max(bounds))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    records = [audit(HERE/f'index_{index:03}') for index in (1, 14, 75, 107, 124)]
    report = {'method': 'independent fixed-order residue-set DFS; all exact duals replayed',
              'classes': records,
              'total_certificates': sum(row['exact_certificates_verified'] for row in records)}
    rendered = json.dumps(report, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered)

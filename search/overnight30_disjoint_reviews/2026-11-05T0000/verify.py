"""Verify universal disjoint-neighborhood bounds and conditional quotient proofs."""
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import runpy
import sys

HERE = Path(__file__).resolve().parent
SEARCH = HERE.parents[1]
sys.path.insert(0, str(SEARCH))
from event30_disjoint_probability import (bound_is_strict, moment_coefficients,
                                        verify_bound_certificate, verify_quotient_certificate)
from overnight30_screen import exact_psd


def h33_control():
    data = (SEARCH / 'h33.g6').read_bytes().splitlines()[0]
    if data.startswith(b'>>graph6<<'):
        data = data[10:]
    assert data[0]-63 == 33
    adjacency = [set() for _ in range(33)]
    for v in range(1, 33):
        for u in range(v):
            bit = v*(v-1)//2+u
            if (data[1+bit//6]-63) & (1 << (5-bit % 6)):
                adjacency[u].add(v)
                adjacency[v].add(u)
    parts = [list(range(9)), list(range(9, 21)), list(range(21, 33))]
    weights = [F(1, 9) if i < 9 else F(1, 12) for i in range(33)]
    K = [[weights[i] if i == j else
          F(1, 36) if j in adjacency[i] and (i < 9 or j < 9) else
          F(1, 24) if j in adjacency[i] else F(0)
          for j in range(33)] for i in range(33)]
    subset = next(subset for subset in itertools.combinations(parts[1], 3)
                  if all(not (adjacency[i] & adjacency[j] & set(parts[0]))
                         for i, j in itertools.combinations(subset, 2)))
    m = len(subset)
    PP = {(i, j): sum((K[i][a]*K[a][j]/weights[a] for a in parts[0]), F(0))
          for i in subset for j in subset}
    QQ = {(i, j): sum((K[i][a]*K[a][j]/weights[a] for a in parts[2]), F(0))
          for i in subset for j in subset}
    PQ = {(i, j): sum((K[i][a]*K[a][b]*K[b][j]/(weights[a]*weights[b])
                      for a in parts[0] for b in parts[2]), F(0))
          for i in subset for j in subset}
    assert all(PP[i, j] == 0 for i in subset for j in subset if i != j)
    assert all(PQ[i, i] == 0 for i in subset)
    parameters = [sum((weights[i] for i in subset), F(0))/m,
                  sum((PP[i, i] for i in subset), F(0))/m,
                  sum((QQ[i, i] for i in subset), F(0))/m,
                  sum((QQ[i, j] for i in subset for j in subset if i != j), F(0))/(m*(m-1)),
                  sum((PQ[i, j] for i in subset for j in subset if i != j), F(0))/(m*(m-1))]
    coefficients = moment_coefficients(m)
    matrices = []
    for size, basis in zip((3, 4), coefficients):
        matrix = [[sum((coefficient*basis[k][i][j] for k, coefficient in enumerate(parameters)), F(0))
                   for j in range(size)] for i in range(size)]
        if size == 4:
            matrix[0][0] = F(1)
        assert exact_psd(matrix)
        matrices.append(matrix)
    certificate = json.loads((HERE / 'disjoint_03_mean_bound.json').read_text())
    assert parameters[0] < verify_bound_certificate(certificate)
    return subset


def main():
    hashes = json.loads((HERE / 'artifact_sha256.json').read_text())
    for name, digest in hashes.items():
        assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == digest
    summary = json.loads((HERE / 'summary.json').read_text())
    source = SEARCH.parent / summary['source_quotient_summary']
    saturation_source = SEARCH.parent / summary['source_saturation_summary']
    assert hashlib.sha256(source.read_bytes()).hexdigest() == summary['source_quotient_sha256']
    assert hashlib.sha256(saturation_source.read_bytes()).hexdigest() == summary['source_saturation_sha256']
    prior, saturation = json.loads(source.read_text()), json.loads(saturation_source.read_text())
    helpers = runpy.run_path(str(source.parent / 'verify.py'))
    canonical, graph_quotient = helpers['canonical'], helpers['graph_quotient']
    assert len(summary['bounds']) == 9
    for row in summary['bounds']:
        certificate = json.loads((HERE / row['certificate']).read_text())
        upper = verify_bound_certificate(certificate)
        assert str(upper) == row['mean_upper'] and str(certificate['subset_size']*upper) == row['total_upper']
        assert bound_is_strict(certificate) == row['strict'] is True
        assert row['forbidden_at30'] == (certificate['subset_size']*upper <= F(certificate['subset_size']-1, 9))
    covered = {canonical(r['quotient']) for r in prior['results'] if r['exact_lp_excluded']}
    covered.update(canonical(r['quotient']) for r in saturation['results'])
    previous = len(covered)
    for row in summary['results']:
        certificate = json.loads((HERE / row['certificate']).read_text())
        assert certificate['quotient'] == row['quotient'] == graph_quotient(row['graph6'])
        assert str(verify_quotient_certificate(certificate)) == row['margin']
        assert row['quotient'] == prior['results'][row['source_index']]['quotient']
        key = canonical(row['quotient'])
        assert key not in covered
        covered.add(key)
    assert len(covered)-previous == summary['certified'] == 10
    assert len(covered) == summary['combined_probability_certified'] == 141
    assert prior['reviewed']-len(covered) == summary['remaining_probability_unresolved'] == 9
    subset = h33_control()
    print('Verified nine universal moment bounds, ten exact quotient certificates, and actual C5 quotient fixtures.')
    print('Known H33 moment control passed on subset', list(subset))


if __name__ == '__main__':
    main()

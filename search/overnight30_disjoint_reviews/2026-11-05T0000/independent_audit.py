"""Independent determinant-based audit of all moment certificates.

This recomputes every principal minor and all six affine coefficients using
separate formulas rather than the production Schur-complement verifier.
"""
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from event30_disjoint_probability import bound_is_strict, verify_bound_certificate


def determinant(matrix):
    n, value = len(matrix), F(0)
    for permutation in itertools.permutations(range(n)):
        sign = (-1)**sum(permutation[i] > permutation[j] for i in range(n) for j in range(i+1, n))
        term = F(sign)
        for i in range(n):
            term *= matrix[i][permutation[i]]
        value += term
    return value


def principal_minors(matrix):
    return [determinant([[matrix[i][j] for j in indices] for i in indices])
            for size in range(1, len(matrix)+1)
            for indices in itertools.combinations(range(len(matrix)), size)]


def coefficients(m, c, p, q, r, t, constant=0):
    S = [[c, p, q-r], [p, p, -t], [q-r, -t, q-r]]
    Q, R = (q+(m-1)*r)/m, (m-1)*t/m
    T = [[F(constant), c, c, c], [c, c/m, p/m, Q], [c, p/m, p/m, R], [c, Q, R, Q]]
    return S, T


def inner(A, B):
    return sum((x*y for ar, br in zip(A, B) for x, y in zip(ar, br)), F(0))


def audit():
    rows = []
    for m in range(2, 11):
        certificate = json.loads((HERE / f'disjoint_{m:02}_mean_bound.json').read_text())
        A = [[F(v) for v in row] for row in certificate['standard_dual']]
        B = [[F(v) for v in row] for row in certificate['trivial_dual']]
        assert A == list(map(list, zip(*A))) and B == list(map(list, zip(*B)))
        am, bm = principal_minors(A), principal_minors(B)
        assert all(v >= 0 for v in am) and all(v > 0 for v in bm)
        upper = F(certificate['mean_probability_upper'])
        for i in range(6):
            basis = [F(int(j == i)) for j in range(6)]
            S, T = coefficients(m, *basis[:5], constant=basis[5])
            assert inner(A, S)+inner(B, T) == (-1 if i == 0 else upper if i == 5 else 0)
        assert verify_bound_certificate(certificate) == upper and bound_is_strict(certificate)
        rows.append({'m': m, 'mean_upper': str(upper), 'total_upper': str(m*upper),
                     'standard_principal_minors_nonnegative': len(am),
                     'trivial_principal_minors_positive': len(bm),
                     'strict': True, 'forbidden_n30': m*upper <= F(m-1, 9)})
    partitions = {'5+5': 2*F(15, 32), '5+4+1': F(15, 32)+F(2, 5)+F(1, 9),
                  '5+3+2': F(15, 32)+F(5, 16)+F(7, 32)}
    assert partitions == {'5+5': F(15, 16), '5+4+1': F(1411, 1440), '5+3+2': F(1)}
    return {'certificates': rows, 'partition_strict_upper_bounds': {k: str(v) for k, v in partitions.items()},
            'proof': 'All nine trivial duals are positive definite by every principal minor; T00=1 makes every moment bound strict.'}


if __name__ == '__main__':
    report = audit()
    assert report == json.loads((HERE / 'independent_audit.json').read_text())
    print('Independent audit passed:198 principal minors,54 affine coefficients, and all partition totals.')

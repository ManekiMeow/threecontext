"""Fast exact theta upper bounds from rounded factors and integer residuals.

Discovery API: certify_from_matrix(edges, D, upper) -> certificate or None.
The numerical symmetric matrix D is a candidate theta dual: entries on the
diagonal and support edges must equal one. Discovery enforces those entries,
rounds the others, and factors a slightly enlarged numerical upper bound.

Verification uses only integers. If X = t I - D = R R^T + E and E + delta I
is diagonally dominant with nonnegative diagonal, then (t+delta) I - D is
positive semidefinite. No numerical eigenvalue or factorization is trusted.
The bound applies to the saved graph and every subgraph on the same vertices.
"""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def residual_slack(matrix, factor, upper_numerator, denominator):
    """Smallest integer diagonal shift making the exact residual DD."""
    n = len(matrix)
    shift = 0
    for i in range(n):
        diagonal, off_sum = 0, 0
        for j in range(n):
            residual = denominator*((upper_numerator if i == j else 0)-matrix[i][j])
            residual -= sum(factor[i][k]*factor[j][k] for k in range(n))
            if i == j:
                diagonal = residual
            else:
                off_sum += abs(residual)
        shift = max(shift, off_sum-diagonal)
    return shift


def verify_certificate(cert):
    require(cert.get("format") == "overnight30-factor-dual-v1", "unknown certificate format")
    n, q, p = cert["n"], cert["denominator"], cert["base_upper_numerator"]
    require(type(n) is int and n > 0, "invalid order")
    require(type(q) is int and q > 0 and type(p) is int, "invalid denominator or bound")
    matrix, factor = cert["dual_numerators"], cert["factor_numerators"]
    for rows in (matrix, factor):
        require(len(rows) == n and all(len(row) == n for row in rows), "invalid matrix shape")
        require(all(type(x) is int for row in rows for x in row), "noninteger matrix entry")
    require(all(matrix[i][j] == matrix[j][i] for i in range(n) for j in range(i)), "asymmetric dual")
    require(all(matrix[i][i] == q for i in range(n)), "dual diagonal is not one")
    edges = cert["edges"]
    require(edges == sorted({tuple(e) for e in edges}) or
            edges == [list(e) for e in sorted({tuple(e) for e in edges})], "edges not normalized")
    for u, v in edges:
        require(type(u) is int and type(v) is int and 0 <= u < v < n, "invalid edge")
        require(matrix[u][v] == q, "dual entry on a support edge is not one")
    slack = residual_slack(matrix, factor, p, q)
    require(type(cert["residual_shift_numerator"]) is int and
            cert["residual_shift_numerator"] == slack, "incorrect residual shift")
    bound = Fraction(p*q+slack, q*q)
    require(bound < 3, "certificate does not exclude theta three")
    require(str(bound) == cert["bound"], "incorrect stated bound")
    return bound


def certify_from_matrix(edges, D, upper, denominator=1000000):
    """Discover a certificate, returning None when rounding/factoring fails."""
    import numpy as np
    D = np.asarray(D, dtype=float)
    if D.ndim != 2 or D.shape[0] != D.shape[1] or not np.isfinite(D).all():
        return None
    n, q = D.shape[0], denominator
    if not math.isfinite(upper) or upper >= 3:
        return None
    require(type(q) is int and q > 0, "invalid discovery denominator")
    normalized_edges = set()
    for u, v in edges:
        require(isinstance(u, (int, np.integer)) and not isinstance(u, (bool, np.bool_)) and
                isinstance(v, (int, np.integer)) and not isinstance(v, (bool, np.bool_)),
                "noninteger discovery endpoint")
        normalized_edges.add(tuple(sorted((int(u), int(v)))))
    edges = sorted(normalized_edges)
    require(all(0 <= u < v < n for u, v in edges), "invalid discovery edge")
    D = (D+D.T)/2
    matrix = [[int(round(float(D[i, j])*q)) for j in range(n)] for i in range(n)]
    for i in range(n):
        matrix[i][i] = q
    for u, v in edges:
        matrix[u][v] = matrix[v][u] = q
    rounded = np.asarray(matrix, dtype=float)/q
    for allowance in (1e-5, 1e-4, 1e-3, 1e-2):
        p = math.ceil((upper+allowance)*q)
        if p >= 3*q:
            continue
        try:
            L = np.linalg.cholesky((p/q)*np.eye(n)-rounded)
        except np.linalg.LinAlgError:
            continue
        factor = [[int(round(float(L[i, j])*q)) for j in range(n)] for i in range(n)]
        slack = residual_slack(matrix, factor, p, q)
        bound = Fraction(p*q+slack, q*q)
        if bound >= 3:
            continue
        cert = {"format": "overnight30-factor-dual-v1", "n": n,
                "edges": [list(e) for e in edges], "denominator": q,
                "base_upper_numerator": p, "dual_numerators": matrix,
                "factor_numerators": factor, "residual_shift_numerator": slack,
                "bound": str(bound), "scope": "saved support and all its edge-deleted subgraphs"}
        verify_certificate(cert)
        return cert
    return None


def self_test():
    import copy
    import numpy as np
    cert = certify_from_matrix([], np.eye(5), 1.0)
    require(cert is not None and verify_certificate(cert) < Fraction(10001, 10000), "identity control failed")
    damaged = copy.deepcopy(cert)
    damaged["residual_shift_numerator"] += 1
    try:
        verify_certificate(damaged)
    except ValueError:
        pass
    else:
        raise RuntimeError("altered residual certificate accepted")
    damaged = copy.deepcopy(cert)
    damaged["edges"] = [[0, 1]]
    try:
        verify_certificate(damaged)
    except ValueError:
        pass
    else:
        raise RuntimeError("altered support certificate accepted")
    require(certify_from_matrix([(0, 1), (0, 2), (1, 2)], np.ones((3, 3)), 2.9) is None,
            "false triangle exclusion accepted")
    print("integer-factor dual controls passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", nargs="?", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    elif args.certificate:
        print("verified theta <=", verify_certificate(json.loads(args.certificate.read_text())))
    else:
        parser.error("provide a certificate or --self-test")

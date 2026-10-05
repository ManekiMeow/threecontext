"""Stand-alone exact verification of the new three-context GHZ graphs.

For a graph H (the *non-orthogonality* graph; the exclusivity graph is G = complement of H)
the paper's Theorem 1 requires  alpha(G) = 2, theta(G) = 3, chi(Hbar-bar) = chi(H) = 3, i.e.
  * H triangle-free,
  * H properly 3-coloured (here: by the class index a),
  * theta(Hbar) = 3.
theta(Hbar) <= chi(H) = 3 always, so only theta(Hbar) >= 3 needs a certificate:
  rational c_i > 0 and rational Y (zero diagonal, supported on E(H)) with
  sum_{j in N(i) of colour k} Y_ij = c_i for every i and every colour k != colour(i),
  and B = Y + diag(c) PSD (checked by exact Gaussian elimination over Q).
All arithmetic below is exact (fractions.Fraction).
"""
import json, itertools
from fractions import Fraction as Fr


def tricayley(q, S, removed=()):
    V = [(a, x) for a in range(3) for x in range(q) if (a, x) not in removed]
    Vs = set(V)
    E = set()
    for a in range(3):
        for x in range(q):
            for s in S[a]:
                u, v = (a, x), ((a + 1) % 3, (x + s) % q)
                if u in Vs and v in Vs: E.add((u, v))
    return V, E


def triangle_free(V, E):
    adj = {v: set() for v in V}
    for u, v in E: adj[u].add(v); adj[v].add(u)
    return all(not (adj[u] & adj[v]) for u, v in E), adj


def exact_psd(A):
    n = len(A); A = [row[:] for row in A]
    for i in range(n):
        if A[i][i] < 0: return False
        if A[i][i] == 0:
            if any(A[i][j] != 0 for j in range(i + 1, n)): return False
            continue
        for j in range(i + 1, n):
            if A[j][i] != 0:
                r = A[j][i] / A[i][i]
                for k in range(i, n): A[j][k] -= r * A[i][k]
    return True


def check_certificate(V, E, Y, c):
    """Y: dict {(u,v): Fraction} on edges, c: dict {v: Fraction}."""
    ok_tri, adj = triangle_free(V, E)
    proper = all(u[0] != v[0] for u, v in E)
    assert all(c[v] > 0 for v in V)
    for v in V:
        for k in range(3):
            if k == v[0]: continue
            s = sum(Y.get((v, u), Y.get((u, v), 0)) for u in adj[v] if u[0] == k)
            assert s == c[v], (v, k, s, c[v])
    idx = {v: i for i, v in enumerate(V)}
    B = [[Fr(0)] * len(V) for _ in V]
    for (u, v), y in Y.items(): B[idx[u]][idx[v]] = B[idx[v]][idx[u]] = y
    for v in V: B[idx[v]][idx[v]] = c[v]
    psd = exact_psd(B)
    total = sum(c.values()) + 2 * sum(Y.values())
    return ok_tri, proper, psd, total / sum(c.values())


def tricayley_z2z6(S):
    els = [(a, z) for a in range(2) for z in range(6)]
    add = lambda g, s: ((g[0] + s[0]) % 2, (g[1] + s[1]) % 6)
    V = [(k, g) for k in range(3) for g in els]
    E = {((k, g), ((k + 1) % 3, add(g, s))) for k in range(3) for g in els for s in S[k]}
    return V, E


if __name__ == "__main__":
    # 36-vertex graph over Z2 x Z6: unweighted.  W = A/4 (each vertex has 4 neighbours in each
    # other class), Y = A/(4*12), c = 1/12.
    S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)),
           ((0, 1), (1, 0), (1, 1), (1, 3)))
    V, E = tricayley_z2z6(S36)
    Y = {e: Fr(1, 48) for e in E}
    c = {v: Fr(1, 12) for v in V}
    print("36-vertex T(Z2xZ6; S0=S1={00,03,04,14}, S2={01,10,11,13}): |V|=%d |E|=%d" % (len(V), len(E)),
          "triangle-free, proper 3-colouring, B PSD, sum(B)/tr(B) =", check_certificate(V, E, Y, c))
    # 39-vertex graph: uniform weights 1/4, 1/3, 1/3, c = 1/13
    q, S = 13, ((0, 1, 3, 9), (0, 1, 10), (1, 6, 8))
    V, E = tricayley(q, S)
    wts = [Fr(1, 4), Fr(1, 3), Fr(1, 3)]
    Y = {(u, v): wts[u[0]] * Fr(1, 13) for (u, v) in E}
    c = {v: Fr(1, 13) for v in V}
    print("39-vertex T(13;{0,1,3,9},{0,1,10},{1,6,8}): |V|=%d |E|=%d" % (len(V), len(E)),
          "triangle-free, proper 3-colouring, B PSD, sum(B)/tr(B) =", check_certificate(V, E, Y, c))
    # 37-vertex graph from the json certificate
    d = json.load(open("shr_q13_1_0139_014_1257_cert.json"))
    S = tuple(tuple(s) for s in d["S"])
    V, E = tricayley(13, S, removed={(0, 5), (2, 11)})
    nodes = [tuple(v) for v in d["nodes"]]
    assert sorted(nodes) == sorted(V)
    Y = {}
    for a, b, y in d["edges"]:
        a, b = tuple(a), tuple(b)
        e = (a, b) if (a, b) in E else (b, a)
        assert e in E
        Y[e] = Fr(y)
    assert len(Y) == len(E)
    c = {tuple(v): Fr(x) for v, x in zip(d["nodes"], d["c"])}
    print("37-vertex T(13;{0,1,3,9},{0,1,4},{1,2,5,7}) minus (0,5),(2,11): |V|=%d |E|=%d" % (len(V), len(E)),
          "triangle-free, proper 3-colouring, B PSD, sum(B)/tr(B) =", check_certificate(V, E, Y, c))

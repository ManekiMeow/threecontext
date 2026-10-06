# Exhaustive search of triangle-free Cayley graphs Cay(Gamma, S) on groups of order 30 and 33
# for three-context GHZ graphs.  For a vertex-transitive H, theta(Hbar) = n/theta(H), so
# theta(Hbar) = 3 needs theta(H) = n/3, and theta(H) <= Hoffman bound = n(-l_min)/(k - l_min)
# needs l_min(A) <= -k/2.  Survivors are checked with the full SDP.
import sys, numpy as np, networkx as nx
from groupsearch import perm_group
from theta import theta_comp


def cyclic(n):
    return perm_group([[(i + 1) % n for i in range(n)]])


def dihedral(m):
    return perm_group([[(i + 1) % m for i in range(m)], [(-i) % m for i in range(m)]])


GROUPS = {
    "Z30": lambda: cyclic(30),
    "Z33": lambda: cyclic(33),
    "D15": lambda: dihedral(15),
    "Z3xD5": lambda: perm_group([[1, 2, 0, 3, 4, 5, 6, 7], [0, 1, 2, 4, 5, 6, 7, 3], [0, 1, 2, 3, 7, 6, 5, 4]]),
    "Z5xS3": lambda: perm_group([[1, 2, 3, 4, 0, 5, 6, 7], [0, 1, 2, 3, 4, 6, 7, 5], [0, 1, 2, 3, 4, 6, 5, 7]]),
}


def search(name, maxdeg):
    mul, e = GROUPS[name](); n = len(mul)
    inv = [next(y for y in range(n) if mul[x][y] == e) for x in range(n)]
    classes = []
    seen = set()
    for x in range(n):
        if x == e or x in seen: continue
        cl = tuple(sorted({x, inv[x]})); seen.update(cl); classes.append(cl)
    hits = []; tested = 0; best = 0

    def tri_free(S):
        Ss = set(S)
        return not any(mul[mul[a][b]][c] == e for a in Ss for b in Ss for c in Ss)

    def evaluate(S):
        nonlocal tested, best
        A = np.zeros((n, n))
        for g in range(n):
            for s in S: A[g, mul[g][s]] = 1
        ev = np.linalg.eigvalsh(A); k = len(S)
        if ev[0] > -k / 2 + 1e-9: return
        if ev[0] < -k + 1e-9: return               # lambda_min = -k: bipartite component
        if ev[-1] + ev[0] < 1e-9: return          # bipartite (symmetric spectrum): chi = 2
        H = nx.from_numpy_array(A)
        if not nx.is_connected(H) or nx.is_bipartite(H): return
        tested += 1
        t, _ = theta_comp(H)
        best = max(best, t)
        if t > 3 - 1e-6:
            hits.append(S); print("HIT", name, S, t, flush=True)

    def dfs(start, S):
        if S: evaluate(S)
        for i in range(start, len(classes)):
            T = S + list(classes[i])
            if len(T) > maxdeg: continue
            if tri_free(T): dfs(i + 1, T)
    dfs(0, [])
    print(name, "order", n, "SDPs", tested, "best theta(Hbar)", best, "hits", len(hits), flush=True)


if __name__ == "__main__":
    search(sys.argv[1], int(sys.argv[2]))

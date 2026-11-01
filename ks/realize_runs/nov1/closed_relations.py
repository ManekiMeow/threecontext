"""Count 'closed' completeness relations inside the 1-ball {u} + N(u): pairs x, x' in N(u) and a basis b
(not containing x or x') such that every common neighbour of x and x' in b lies in N(u)."""
import os
import sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))  # ks/
from realize_search import load_pattern
X = load_pattern(sys.argv[1]); n = len(X); d = n // 4; part = np.arange(n) // d
u = int(sys.argv[2]) if len(sys.argv) > 2 else 0
N = np.where(X[u])[0]; inN = np.zeros(n, bool); inN[N] = True
stats = {}
for a in range(len(N)):
    for b_ in range(a + 1, len(N)):
        x, xp = N[a], N[b_]
        kind = "adj" if X[x, xp] else ("same" if part[x] == part[xp] else "orth")
        for b in range(4):
            if b in (part[x], part[xp]):
                continue
            com = np.where(X[x] & X[xp] & (part == b))[0]
            if b == part[u]:
                closed = len(com) > 0 and all(c == u for c in com)
            else:
                closed = len(com) > 0 and inN[com].all()
            key = (kind, "u-basis" if b == part[u] else "local", len(com), bool(closed))
            stats[key] = stats.get(key, 0) + 1
for k in sorted(stats):
    print(k, stats[k])

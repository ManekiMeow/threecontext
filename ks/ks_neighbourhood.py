# For the 21-ray KS set (d=6): take a ray w as the state and keep the rays non-orthogonal to w
# (the "GHZ-type paradox derived from a KS set").  Is the result a GHZ graph in the sense of
# Theorem 1 (alpha(G) = n-1, theta(G) = n, n = clique-cover number)?
import sys, itertools, numpy as np, networkx as nx
sys.path.insert(0, "../search")
from ks21 import rays, G, ip, Z, bases
import cvxpy as cp
def theta(Gr):   # Lovasz theta of the exclusivity graph Gr (edges = orthogonal pairs)
    n = Gr.number_of_nodes(); idx = {v: i for i, v in enumerate(Gr)}
    B = cp.Variable((n, n), symmetric=True)
    cons = [B >> 0, cp.trace(B) == 1] + [B[idx[u], idx[v]] == 0 for u, v in Gr.edges()]
    return cp.Problem(cp.Maximize(cp.sum(B)), cons).solve(solver="CLARABEL")
w = 0
N = [v for v in G if v != w and v not in G[w]]                  # rays non-orthogonal to w
Gw = G.subgraph(N).copy()                                       # exclusivity graph of the derived paradox
alpha = max(len(c) for c in nx.find_cliques(nx.complement(Gw)))
cover = next(k for k in range(1, 11) if any(
    all(nx.complement(Gw).subgraph(p).number_of_edges() == 0 for p in parts) for parts in []) ) if False else None
# clique cover number = chromatic number of complement
Hc = nx.complement(Gw)
def colourable(Gr, k):
    order = list(Gr); col = {}
    def bt(i):
        if i == len(order): return True
        v = order[i]
        for c in range(k):
            if all(col.get(u) != c for u in Gr[v]):
                col[v] = c
                if bt(i + 1): return True
                del col[v]
        return False
    return bt(0)
cover = next(k for k in range(1, 11) if colourable(Hc, k))
ctx = [[v for v in b if v in N] for b in [[rays.index(x) for x in bb] for bb in bases]]
ctx = [c for c in ctx if c and not (w in c)]
print("rays non-orthogonal to w:", len(N))
print("contexts restricted to them (each sums to 1 in state w):", [len(c) for c in ctx])
print("alpha(G_w) =", alpha, "  theta(G_w) = %.4f" % theta(Gw), "  clique-cover number =", cover)
print("quantum sum over these events in state w = %.4f" % sum(abs(complex(*[0,0]) ) for _ in []) if False else "")

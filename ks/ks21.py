# Seven-context 21-ray KS set in d=6 (Lisonek, Badziag, Portillo, Cabello, PRA 89, 042101 (2014)).
# Exact arithmetic in Z[w], w = exp(2 pi i/3): elements a + b w with integers a, b.
import itertools, networkx as nx
W = "w"
def m(x, y):  # multiply (a,b)*(c,d) with w^2 = -1 - w
    a, b = x; c, d = y
    return (a * c - b * d, a * d + b * c - b * d)
def conj(x):  # conj(a + b w) = a + b w^2 = (a - b) - b w
    a, b = x; return (a - b, -b)
def ip(u, v):
    s = (0, 0)
    for x, y in zip(u, v):
        p = m(conj(x), y); s = (s[0] + p[0], s[1] + p[1])
    return s
Z, O, w, w2 = (0, 0), (1, 0), (0, 1), (-1, -1)
def V(s):
    return tuple({"0": Z, "1": O, "w": w, "W": w2}[c] for c in s)   # W = w^2
B = [
 ["100000","010000","001000","000100","000010","000001"],
 ["100000","001111","0101wW","0110Ww","01wW01","01Ww10"],
 ["010000","001111","1001Ww","1010wW","10Ww01","10wW10"],
 ["001000","0101wW","1001Ww","110011","wW0101","Ww0110"],
 ["000100","0110Ww","1010wW","110011","Ww1001","wW1010"],
 ["000010","01wW01","10Ww01","wW0101","Ww1001","111100"],
 ["000001","01Ww10","10wW10","Ww0110","wW1010","111100"],
]
bases = [[V(s) for s in b] for b in B]
for b in bases:   # sanity: each basis orthogonal
    assert all(ip(u, v) == Z for u, v in itertools.combinations(b, 2)), b
rays = sorted({v for b in bases for v in b})
print("rays:", len(rays))
G = nx.Graph(); G.add_nodes_from(range(len(rays)))
for i, j in itertools.combinations(range(len(rays)), 2):
    if ip(rays[i], rays[j]) == Z: G.add_edge(i, j)
cl = [c for c in nx.find_cliques(G)]
print("orthogonal pairs:", G.number_of_edges(), " (line graph of K7 has 105)")
print("maximal cliques sizes:", sorted(len(c) for c in cl))
# clique cover number = chromatic number of complement
from itertools import product
Hc = nx.complement(G)
def colourable(Gr, k):
    n = Gr.number_of_nodes(); order = sorted(Gr, key=lambda v: -Gr.degree(v)); col = {}
    def bt(i):
        if i == n: return True
        v = order[i]
        for c in range(k):
            if all(col.get(u) != c for u in Gr[v]):
                col[v] = c
                if bt(i + 1): return True
                del col[v]
        return False
    return bt(0)
for k in range(3, 8):
    if colourable(Hc, k): print("clique cover number of orthogonality graph:", k); break

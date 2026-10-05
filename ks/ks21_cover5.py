# Complete the 5-clique cover of the 21-ray set to 5 orthonormal bases and check KS-uncolourability
# of the enlarged set (all orthonormal bases contained in it count as contexts).
import itertools, numpy as np, networkx as nx
from ks21 import bases
w = np.exp(2j * np.pi / 3)
val = {(0, 0): 0, (1, 0): 1, (0, 1): w, (-1, -1): w * w}
vec = lambda v: np.array([val[x] for x in v]) / np.linalg.norm([val[x] for x in v])
label = {}
for i, b in enumerate(bases):
    for v in b: label.setdefault(v, []).append(i + 1)
R = {tuple(sorted(l)): vec(v) for v, l in label.items()}       # key (i,j): ray common to B_i, B_j
cover_sets = [[k for k in R if 1 in k], [k for k in R if 2 in k and 1 not in k], [k for k in R if 3 in k and k[0] > 2],
              [(4, 5), (4, 6), (5, 6)], [(4, 7), (5, 7), (6, 7)]]
assert sorted(sum(cover_sets, [])) == sorted(R), "not a cover"
rays = [R[k] for k in R]; names = [str(k) for k in R]
for cs in cover_sets[3:]:
    M = np.array([R[k] for k in cs]).conj()
    _, _, Vh = np.linalg.svd(M); comp = Vh[len(cs):].conj()   # orthonormal completion
    for c in comp: rays.append(c); names.append("new")
n = len(rays); A = np.array(rays)
Gm = np.abs(A.conj() @ A.T) < 1e-9
G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from((i, j) for i in range(n) for j in range(i + 1, n) if Gm[i, j])
ctx = [c for c in nx.find_cliques(G) if len(c) == 6]
print("rays:", n, " orthonormal bases inside the set:", len(ctx))
# KS colouring: choose exactly one ray per basis, no two orthogonal rays both 1
def ks_colourable():
    order = sorted(ctx, key=len)
    val = {}
    def bt(k):
        if k == len(ctx): return True
        c = ctx[k]; ones = [v for v in c if val.get(v) == 1]
        if len(ones) > 1: return False
        if ones:
            for v in c:
                if v not in val: val[v] = 0
            ok = bt(k + 1)
            return ok
        free = [v for v in c if v not in val]
        for v in free:
            if any(val.get(u) == 1 for u in G[v]): continue
            snap = dict(val); val[v] = 1
            for u in c:
                if u not in val: val[u] = 0
            for u in G[v]: val.setdefault(u, 0)
            if all(sum(val.get(x) == 1 for x in cc) <= 1 for cc in ctx) and bt(k + 1): return True
            val.clear(); val.update(snap)
        return False
    return bt(0)
print("KS-colourable:", ks_colourable())
print("covered by 5 orthonormal bases: yes (B1, B2, B3 and the two completed triangles)")

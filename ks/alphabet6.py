# Rays in C^6 with entries in {0,1,w,w^2} (first non-zero entry = 1), exact orthogonality in Z[w].
import itertools, pickle, networkx as nx
from ks21 import ip, Z
E = [(0, 0), (1, 0), (0, 1), (-1, -1)]          # 0, 1, w, w^2
rays = []
for v in itertools.product(range(4), repeat=6):
    nz = [i for i in v if i]
    if not nz or nz[0] != 1: continue
    rays.append(tuple(E[i] for i in v))
n = len(rays); print("rays", n)
G = nx.Graph(); G.add_nodes_from(range(n))
for i in range(n):
    for j in range(i + 1, n):
        if ip(rays[i], rays[j]) == Z: G.add_edge(i, j)
print("orthogonal pairs", G.number_of_edges())
bases = [tuple(sorted(c)) for c in nx.find_cliques(G) if len(c) == 6]
# find_cliques gives maximal cliques; a basis is a maximal 6-clique (6 = max possible)
print("bases (6-cliques)", len(bases))
pickle.dump((rays, G, bases), open("alph6.pkl", "wb"))

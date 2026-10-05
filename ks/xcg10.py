# The 10-event GHZ-type proof of Xu-Chen-Guhne (arXiv:2001.07656), derived from the 18-vector CEG set
# with state psi_0.  Compute its exclusivity graph, clique-cover number, alpha and theta.
import itertools, numpy as np, networkx as nx, cvxpy as cp
v = {0:"1000",1:"0001",2:"0110",3:"01-10",4:"1001",5:"111-1",6:"-1111",7:"11-11",8:"1010",9:"010-1",
     10:"10-10",11:"1-11-1",12:"1111",13:"11-1-1",14:"1-100",15:"001-1",16:"0011",17:"0100"}
def parse(s):
    out, i = [], 0
    while i < len(s):
        if s[i] == "-": out.append(-int(s[i + 1])); i += 2
        else: out.append(int(s[i])); i += 1
    return np.array(out, float)
V = {k: parse(s) for k, s in v.items()}
ctx = [[0,1,2,3],[3,4,5,6],[6,7,8,9],[9,10,11,12],[12,13,14,15],[15,16,17,0],[17,1,8,10],[2,4,11,13],[5,16,14,7]]
for c in ctx: assert all(abs(V[a] @ V[b]) < 1e-12 for a, b in itertools.combinations(c, 2)), c
psi = V[0] / np.linalg.norm(V[0])
ev = [k for k in V if k != 0 and abs(V[k] @ psi) > 1e-12]
print("events non-orthogonal to psi_0:", ev, " (paper: 4-8, 10-14)")
G = nx.Graph(); G.add_nodes_from(ev)
G.add_edges_from((a, b) for a, b in itertools.combinations(ev, 2) if abs(V[a] @ V[b]) < 1e-12)
Hc = nx.complement(G)
def colourable(Gr, k):
    order = list(Gr); col = {}
    def bt(i):
        if i == len(order): return True
        x = order[i]
        for c in range(k):
            if all(col.get(u) != c for u in Gr[x]):
                col[x] = c
                if bt(i + 1): return True
                del col[x]
        return False
    return bt(0)
cover = next(k for k in range(1, 11) if colourable(Hc, k))
alpha = max(len(c) for c in nx.find_cliques(Hc))
n = len(ev); idx = {x: i for i, x in enumerate(ev)}
B = cp.Variable((n, n), symmetric=True)
th = cp.Problem(cp.Maximize(cp.sum(B)), [B >> 0, cp.trace(B) == 1] + [B[idx[a], idx[b]] == 0 for a, b in G.edges()]).solve(solver="CLARABEL")
probs = {k: round((V[k] @ psi) ** 2 / (V[k] @ V[k]), 4) for k in ev}
print("clique-cover number:", cover, " alpha:", alpha, " theta: %.4f" % th, " dimension: 4")
print("quantum probabilities in psi_0:", probs, " total:", round(sum(probs.values()), 4))

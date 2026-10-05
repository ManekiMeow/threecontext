# dim(S_k ∩ S_l) for the GHZ realisations (H36, H39, Perkel: unique/uniform; H37: interior-point optimum)
import sys, json, numpy as np, networkx as nx
sys.path.insert(0, "../search")
from ghzgeom import dim as _dim
dim = lambda M: _dim(M, tol=1e-5)
from theta import theta_comp
def report(name, H, cls_of):
    nodes = sorted(H); n = len(nodes); G = nx.Graph(); G.add_nodes_from(range(n))
    idx = {v: i for i, v in enumerate(nodes)}; G.add_edges_from((idx[u], idx[v]) for u, v in H.edges())
    t, B = theta_comp(G)
    w, U = np.linalg.eigh(B); keep = w > 1e-7 * w.max(); X = U[:, keep] * np.sqrt(w[keep])
    cls = [cls_of(v) for v in nodes]
    S = [X[[i for i in range(n) if cls[i] == k]] for k in range(3)]
    out = [f"{name}: theta={t:.6f} rank={X.shape[1]} |I_k|={[len(s) for s in S]}"]
    for k, l in [(0, 1), (0, 2), (1, 2)]:
        out.append(f"dim(S{k}∩S{l})={dim(S[k]) + dim(S[l]) - dim(np.vstack([S[k], S[l]]))}")
    print("  ".join(out), flush=True)
from verify_all import tricayley_z2z6, tricayley
from perkel import perkel
S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
V, E = tricayley_z2z6(S36); H = nx.Graph(); H.add_edges_from(E); report("H36", H, lambda v: v[0])
V, E = tricayley(13, ((0, 1, 3, 9), (0, 1, 10), (1, 6, 8))); H = nx.Graph(); H.add_edges_from(E); report("H39", H, lambda v: v[0])
V, E = tricayley(13, ((0, 1, 3, 9), (0, 1, 4), (1, 2, 5, 7)), removed={(0, 5), (2, 11)}); H = nx.Graph(); H.add_edges_from(E); report("H37", H, lambda v: v[0])
report("Perkel", perkel(), lambda v: v[0])

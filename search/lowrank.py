# Search for a low-dimensional real realisation of a three-context GHZ graph H:
# unit vectors g_i in R^r with g_i . g_j = 0 for every non-edge (i != j) of H, and a unit psi
# with sum_{i in I_k} (psi . g_i)^2 = 1 for each colour class k.
import sys, json, torch, numpy as np, networkx as nx
torch.set_default_dtype(torch.float64)


def load(fn):
    d = json.load(open(fn)); H = nx.Graph()
    for a, b, _ in d["edges"]: H.add_edge(tuple(a), tuple(b))
    return H


def run(H, r, seed, steps=3000):
    nodes = sorted(H); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
    mask = torch.ones(n, n)
    for u, v in H.edges(): mask[idx[u], idx[v]] = mask[idx[v], idx[u]] = 0
    mask.fill_diagonal_(0)
    cls = torch.tensor([v[0] for v in nodes])
    torch.manual_seed(seed)
    X = torch.randn(n, r, requires_grad=True); p = torch.randn(r, requires_grad=True)

    def loss():
        G = X / X.norm(dim=1, keepdim=True); psi = p / p.norm()
        M = G @ G.T
        a2 = (G @ psi) ** 2
        return (mask * M ** 2).sum() + sum((1 - a2[cls == k].sum()) ** 2 for k in range(3))
    opt = torch.optim.Adam([X, p], lr=0.02)
    for _ in range(steps):
        opt.zero_grad(); L = loss(); L.backward(); opt.step()
    opt = torch.optim.LBFGS([X, p], max_iter=5000, tolerance_grad=1e-16, tolerance_change=1e-20,
                            line_search_fn="strong_wolfe")

    def closure():
        opt.zero_grad(); L = loss(); L.backward(); return L
    for _ in range(4): opt.step(closure)
    return loss().item()


if __name__ == "__main__":
    H = load(sys.argv[1])
    for r in map(int, sys.argv[2].split(",")):
        best = min(run(H, r, s) for s in range(int(sys.argv[3])))
        print("r", r, "best loss", best, flush=True)

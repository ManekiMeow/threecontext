import torch, numpy as np, networkx as nx, lowrank
from lowrank import load
H = load("shr_q13_1_0139_014_1257_cert.json")
nodes = sorted(H); n = len(nodes); idx = {v: i for i, v in enumerate(nodes)}
best = None
for seed in range(3):
    torch.manual_seed(seed)
    # reuse lowrank.run internals by re-implementing to keep the vectors
    mask = torch.ones(n, n)
    for u, v in H.edges(): mask[idx[u], idx[v]] = mask[idx[v], idx[u]] = 0
    mask.fill_diagonal_(0); cls = torch.tensor([v[0] for v in nodes])
    X = torch.randn(n, 24, requires_grad=True); p = torch.randn(24, requires_grad=True)
    def loss():
        G = X / X.norm(dim=1, keepdim=True); psi = p / p.norm(); a2 = (G @ psi) ** 2
        return (mask * (G @ G.T) ** 2).sum() + sum((1 - a2[cls == k].sum()) ** 2 for k in range(3))
    opt = torch.optim.Adam([X, p], lr=0.02)
    for _ in range(3000): opt.zero_grad(); L = loss(); L.backward(); opt.step()
    opt = torch.optim.LBFGS([X, p], max_iter=5000, tolerance_grad=1e-16, tolerance_change=1e-20, line_search_fn="strong_wolfe")
    def cl():
        opt.zero_grad(); L = loss(); L.backward(); return L
    for _ in range(6): opt.step(cl)
    if best is None or loss().item() < best[0]: best = (loss().item(), X.detach().clone(), p.detach().clone())
L, X, p = best
G = (X / X.norm(dim=1, keepdim=True)).numpy(); psi = (p / p.norm()).numpy()
a = G @ psi; B = (a[:, None] * G) @ (a[:, None] * G).T / 3   # h_i = a_i g_i, tr B = 1
M = (G @ G.T); nonedge = max(abs(M[i, j]) for i in range(n) for j in range(i + 1, n) if not H.has_edge(nodes[i], nodes[j]))
w = np.linalg.eigvalsh(B)
print("loss %.2e  max |<g_i|g_j>| on non-edges %.1e  context sums" % (L, nonedge),
      [round(float((a[[v[0] == k for v in nodes]] ** 2).sum()), 8) for k in range(3)],
      " rank(B) =", int((w > 1e-8 * w.max()).sum()), " min |a_i|^2 = %.4f" % (a ** 2).min())

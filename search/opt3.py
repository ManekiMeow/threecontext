# Continuous search for three-context GHZ paradoxes.
# Three orthonormal sets U1,U2,U3 (d x m_k) in R^d whose spans all contain psi=e0.
# Loss = sum over transversal triples of |<g_i|g_j>|^2 |<g_j|g_l>|^2 |<g_i|g_l>|^2  (zero <=> triangle-free)
#      + penalty forcing psi into each span.
import sys, torch, numpy as np
torch.set_default_dtype(torch.float64)


def orth(X):
    Q, R = torch.linalg.qr(X)
    return Q * torch.sign(torch.diagonal(R))


def loss_fn(Xs, wpen=10.0):
    U = [orth(X) for X in Xs]
    G12 = (U[0].T @ U[1]) ** 2
    G23 = (U[1].T @ U[2]) ** 2
    G13 = (U[0].T @ U[2]) ** 2
    tri = torch.sum(G23 * (G12.T @ G13))
    pen = sum((1 - torch.sum(u[0] ** 2)) for u in U)
    return tri, pen, U


def run(d, ms, seed, steps=6000, lr=0.02):
    torch.manual_seed(seed)
    Xs = [torch.randn(d, m, requires_grad=True) for m in ms]
    opt = torch.optim.Adam(Xs, lr=lr)
    for it in range(steps):
        opt.zero_grad()
        tri, pen, _ = loss_fn(Xs)
        L = tri + 10 * pen
        L.backward()
        opt.step()
    # polish with LBFGS
    opt = torch.optim.LBFGS(Xs, max_iter=2000, tolerance_grad=1e-14, tolerance_change=1e-16,
                            line_search_fn="strong_wolfe")

    def closure():
        opt.zero_grad()
        tri, pen, _ = loss_fn(Xs)
        L = tri + 100 * pen
        L.backward()
        return L
    for _ in range(5):
        opt.step(closure)
    tri, pen, U = loss_fn(Xs)
    return tri.item(), pen.item(), [u.detach().numpy() for u in U]


if __name__ == "__main__":
    d = int(sys.argv[1]); ms = list(map(int, sys.argv[2:5])); seeds = int(sys.argv[5])
    for s in range(seeds):
        tri, pen, U = run(d, ms, s)
        print(d, ms, "seed", s, "tri", tri, "pen", pen, flush=True)
        np.save(f"sol_{d}_{'_'.join(map(str, ms))}_{s}.npy", np.concatenate(U, axis=1))

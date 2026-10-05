# Numerical search for a four-context Kochen-Specker set:
# four orthonormal bases of R^D (or C^D) such that no 4-tuple (one vector per basis)
# is pairwise non-orthogonal (the 4-partite non-orthogonality graph is K4-free).
# Loss = sum over 4-tuples of the product of the six squared overlaps.
import sys, torch, numpy as np
torch.set_default_dtype(torch.float64)


def unitary(X, cplx):
    if cplx:
        Z = torch.complex(X[0], X[1])
    else:
        Z = X
    Q, R = torch.linalg.qr(Z)
    return Q


def loss(Xs, cplx):
    U = [unitary(X, cplx) for X in Xs]
    G = {}
    for a in range(4):
        for b in range(a + 1, 4):
            G[a, b] = (U[a].conj().T @ U[b]).abs() ** 2
    return torch.einsum("ab,ac,ad,bc,bd,cd->", G[0, 1], G[0, 2], G[0, 3], G[1, 2], G[1, 3], G[2, 3]), U


def run(D, seed, cplx, steps=4000):
    torch.manual_seed(seed)
    shape = (2, D, D) if cplx else (D, D)
    Xs = [torch.randn(*shape, requires_grad=True) for _ in range(4)]
    opt = torch.optim.Adam(Xs, lr=0.01)
    for it in range(steps):
        opt.zero_grad(); L, _ = loss(Xs, cplx); L.backward(); opt.step()
    opt = torch.optim.LBFGS(Xs, max_iter=3000, tolerance_grad=1e-15, tolerance_change=1e-18,
                            line_search_fn="strong_wolfe")

    def closure():
        opt.zero_grad(); L, _ = loss(Xs, cplx); L.backward(); return L
    for _ in range(3):
        opt.step(closure)
    L, U = loss(Xs, cplx)
    return L.item()


if __name__ == "__main__":
    cplx = sys.argv[1] == "c"
    for D in map(int, sys.argv[2].split(",")):
        for s in range(int(sys.argv[3])):
            # normalise by the value for Haar-random bases ~ D^4 * D^-6 = D^-2
            print(D, "complex" if cplx else "real", "seed", s, "loss", run(D, s, cplx), "random-scale", D ** -2.0, flush=True)

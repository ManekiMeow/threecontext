# Route 1 with a (non-abelian) group Gamma acting regularly on all four bases.
# B_4 = {delta_g}, B_k = {lambda_g b_k : g in Gamma}, b_k unitary in C[Gamma] (b~_k * b_k = delta_e).
# F_kl = b~_k * b_l,  F_k4(g) = conj(b_k(g^-1)).  Rays (k,g),(l,h) non-orthogonal iff F_kl(g^-1 h) != 0.
# Loss = sum_{g2,g3,g4} prod over the 6 pairs of |F|^2 with g1 = e (K4 count, smooth) + unitarity penalty.
import sys, numpy as np, torch
sys.path.insert(0, "../search")
from groupsearch import perm_group
torch.set_default_dtype(torch.float64)

def cyc(n): return [(i + 1) % n for i in range(n)]
def semidirect(p, r, q=3):
    # Z_p x| Z_q acting by x -> r x, as permutations of Z_p (affine maps x -> r^j x + i)
    return perm_group([[(x + 1) % p for x in range(p)], [(r * x) % p for x in range(p)]])
def sl23():
    pts = [(a, b) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
    act = lambda M: [pts.index(((M[0][0] * a + M[0][1] * b) % 3, (M[1][0] * a + M[1][1] * b) % 3)) for a, b in pts]
    return perm_group([act([[1, 1], [0, 1]]), act([[0, 2], [1, 0]])])
def heis27():
    # upper unitriangular 3x3 over F_3 acting on F_3^3 (27 points) by x -> M x  (faithful; regular orbit not needed)
    pts = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]
    def act(M):
        return [pts.index(tuple(sum(M[i][j] * p[j] for j in range(3)) % 3 for i in range(3))) for p in pts]
    return perm_group([act([[1, 1, 0], [0, 1, 0], [0, 0, 1]]), act([[1, 0, 0], [0, 1, 1], [0, 0, 1]])])
GROUPS = {
    "SL(2,3)": lambda: sl23(),
    "D13": lambda: perm_group([cyc(13), [(-i) % 13 for i in range(13)]]),
    "Heis27": lambda: heis27(),
    "S4": lambda: perm_group([[1, 2, 3, 0], [1, 0, 2, 3]]),
    "D12": lambda: perm_group([cyc(12), [(-i) % 12 for i in range(12)]]),
    "Z13:Z3 (H39)": lambda: semidirect(13, 3),
    "S3xS3": lambda: perm_group([[1, 2, 0, 3, 4, 5], [1, 0, 2, 3, 4, 5], [0, 1, 2, 4, 5, 3], [0, 1, 2, 4, 3, 5]]),
    "A4xZ3": lambda: perm_group([[1, 2, 0, 3, 4, 5, 6], [1, 0, 3, 2, 4, 5, 6], [0, 1, 2, 3, 5, 6, 4]]),
    "Z7:Z3xZ2": lambda: perm_group([[(x + 1) % 7 for x in range(7)] + [7, 8], [(2 * x) % 7 for x in range(7)] + [7, 8], list(range(7)) + [8, 7]]),
    "GL(2,3)-like S4xZ2": lambda: perm_group([[1, 2, 3, 0, 4, 5], [1, 0, 2, 3, 4, 5], [0, 1, 2, 3, 5, 4]]),
    "Z19:Z3 (Perkel)": lambda: semidirect(19, 7),
    "A5": lambda: perm_group([[1, 2, 3, 4, 0], [1, 2, 0, 3, 4]]),
}

def run(name, seed, steps=1200, lr=0.03):
    mul, e = GROUPS[name](); D = len(mul)
    mul = torch.tensor(mul); inv = torch.tensor([int((mul[x] == e).nonzero()[0]) for x in range(D)])
    L = mul[inv]                                  # L[g,h] = g^-1 h
    Q = mul[inv].T                                # Q[g,x] = x^-1 g   (matrix of f -> f * b is b[Q])
    torch.manual_seed(seed)
    P = [torch.randn(2, D, requires_grad=True) for _ in range(3)]
    def unitary_elements(P):
        out = []
        for p in P:
            b = torch.complex(p[0], p[1]); M = b[Q]
            X = M / torch.linalg.matrix_norm(M)         # Newton-Schulz iteration -> polar factor (stays in the
            for _ in range(40):                         # convolution algebra, differentiable, no SVD)
                X = 1.5 * X - 0.5 * X @ X.conj().T @ X
            U = X
            out.append(U[:, e])                         # u = U delta_e, unitary element of C[Gamma]
        return out
    def F_all(b):
        F = {}
        for k in range(3):
            F[k, 3] = b[k][inv].conj()
            for l in range(k + 1, 3): F[k, l] = b[k].conj() @ b[l][mul]      # (b~_k * b_l)(g) = sum_x conj b_k(x) b_l(xg)
        return F
    def k4(F):
        A = {kl: F[kl].abs() ** 2 for kl in F}
        return torch.einsum("b,c,d,bc,bd,cd->", A[0, 1], A[0, 2], A[0, 3], A[1, 2][L], A[1, 3][L], A[2, 3][L])
    opt = torch.optim.Adam(P, lr=lr)
    first = None
    for it in range(steps):
        opt.zero_grad(); l = k4(F_all(unitary_elements(P)))
        if first is None: first = l.item()
        l.backward(); opt.step()
    b = unitary_elements(P); F = F_all(b)
    unit = max(((bk.conj() @ bk[mul] - torch.nn.functional.one_hot(torch.tensor(e), D)).abs().max().item()) for bk in b)
    supp = {kl: int((F[kl].abs() > 1e-4).sum()) for kl in F}
    return D, first, k4(F).item(), unit, supp

if __name__ == "__main__":
    for name in sys.argv[1].split("|"):
        try:
            mul, e = GROUPS[name](); D = len(mul)
            ab = all(mul[x][y] == mul[y][x] for x in range(D) for y in range(D))
        except Exception as ex:
            print(name, "group build failed", ex); continue
        for s in range(int(sys.argv[2])):
            D, l0, l, u, supp = run(name, s)
            print(f"{name} |G|={D} seed {s}: K4-loss {l0:.3e} -> {l:.3e} (ratio {l/l0:.2e})  unitarity err {u:.1e}  |T_kl| {sorted(supp.values())}", flush=True)

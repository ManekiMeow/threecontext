# Explicit coordinates for the unique realization of H33 (RESULTS.md 1.7) in C^20, adapted to its symmetry.
# Points a in Z3^2, classes i = 0..3 of parallel lines with normals k_i, line (i, c) = {a : k_i.a = c}.
# C^20 = K0 (span u_0..u_3, the four class sums) + 8 two-dimensional blocks C^2_k (k != 0), block basis (chi_k, eps_k).
#   point a      = (1/6) sum_i u_i + (1/3) sum_{k != 0} w^{-k.a} chi_k                 (w = e^{2 pi i/3})
#   B-line (i,c) = (1/sqrt3) u_i + (1/sqrt3) sum_{t=+-1} w^{-tc} e_{t k_i},   e_k = chi_k/2 + (sqrt3/2) eps_k
#   C-line (i,c) = (1/sqrt3) u_i + (1/sqrt3) sum_{t=+-1} w^{-tc} f_{t k_i},   f_k = chi_k/2 - (sqrt3/2) eps_k
#   psi          = (1/2) sum_i u_i.
# In each block the directions chi (0 deg), gamma (30), e (60), eps (90), f (120), beta (150) are the lines of S_1,
# S_3-completion, S_2, S_1-completion, S_3, S_2-completion.
import itertools, numpy as np
w = np.exp(2j * np.pi / 3)
pts = list(itertools.product(range(3), repeat=2))
K = [(1, 0), (0, 1), (1, 1), (1, 2)]                          # class normals
nz = [k for k in pts if k != (0, 0)]
dot = lambda k, a: (k[0] * a[0] + k[1] * a[1]) % 3
neg = lambda k: ((-k[0]) % 3, (-k[1]) % 3)
cls_of = {}
for i, k in enumerate(K):
    cls_of[k] = (i, 1); cls_of[neg(k)] = (i, -1)                # block k = t * k_i
D = 20
U0 = {i: i for i in range(4)}                                 # coordinate of u_i
CHI = {k: 4 + 2 * j for j, k in enumerate(nz)}; EPS = {k: 5 + 2 * j for j, k in enumerate(nz)}
def vec(**parts):
    v = np.zeros(D, complex)
    for idx, c in parts.items(): v[idx] += c
    return v
def unit(i): v = np.zeros(D, complex); v[i] = 1; return v
u = [unit(U0[i]) for i in range(4)]
chi = {k: unit(CHI[k]) for k in nz}; eps = {k: unit(EPS[k]) for k in nz}
s3 = np.sqrt(3)
e = {k: chi[k] / 2 + s3 / 2 * eps[k] for k in nz}; f = {k: chi[k] / 2 - s3 / 2 * eps[k] for k in nz}
beta = {k: s3 / 2 * chi[k] - eps[k] / 2 for k in nz}; gamma = {k: s3 / 2 * chi[k] + eps[k] / 2 for k in nz}
psi = sum(u) / 2
point = {a: sum(u) / 6 + sum(w ** (-dot(k, a)) * chi[k] for k in nz) / 3 for a in pts}
def kt(i, t): return K[i] if t == 1 else neg(K[i])
Bline = {(i, c): (u[i] + sum(w ** (-t * c) * e[kt(i, t)] for t in (1, -1))) / s3 for i in range(4) for c in range(3)}
Cline = {(i, c): (u[i] + sum(w ** (-t * c) * f[kt(i, t)] for t in (1, -1))) / s3 for i in range(4) for c in range(3)}
A = [point[a] for a in pts]; B = [Bline[l] for l in sorted(Bline)]; C = [Cline[l] for l in sorted(Cline)]

if __name__ == "__main__":
    V = A + B + C; G = np.array([[np.vdot(x, y) for y in V] for x in V])
    on = lambda S: np.allclose(np.array([[np.vdot(x, y) for y in S] for x in S]), np.eye(len(S)))
    assert on(A) and on(B) and on(C)
    for a in pts:
        for (i, c) in sorted(Bline):
            inc = dot(K[i], a) == c
            assert np.isclose(np.vdot(point[a], Bline[(i, c)]), 1 / np.sqrt(12) if inc else 0)
            assert np.isclose(np.vdot(point[a], Cline[(i, c)]), 1 / np.sqrt(12) if inc else 0)
    for l, m in itertools.product(sorted(Bline), repeat=2):
        par = l[0] == m[0] and l[1] != m[1]
        assert np.isclose(np.vdot(Bline[l], Cline[m]), 0.5 if par else 0)
    for X in (A, B, C): assert np.isclose(sum(abs(np.vdot(x, psi)) ** 2 for x in X), 1)
    print("explicit model reproduces the H33 Gram matrix; contexts orthonormal; psi has probability 1 in each")

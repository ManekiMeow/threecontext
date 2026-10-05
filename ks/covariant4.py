# Group-covariant four-basis search over G = Z_2^n.
# Basis a = { X^g v_a : g in G },  v_a = (1/sqrt D) sum_chi u_a(chi) |chi>,  u_a = (-1)^{g_a(chi)} * i^{h_a(chi)}.
# <X^g v_a | X^h v_b> = FT(conj(u_a) u_b)(h - g)/D  ->  edge set T_ab = supp FT(conj(u_a) u_b).
# The union is a KS set (4 contexts) if the 4-partite Cayley graph has no K4.
import sys, random, itertools, numpy as np
n = int(sys.argv[1]); seed = int(sys.argv[2]); iters = int(sys.argv[3]); maxdeg = int(sys.argv[4]) if len(sys.argv) > 4 else 3
D = 2 ** n; rng = random.Random(seed)
X = np.array([[(x >> i) & 1 for i in range(n)] for x in range(D)])
monos = [m for k in range(1, maxdeg + 1) for m in itertools.combinations(range(n), k)]
H = np.array([[(-1) ** bin(x & y).count("1") for y in range(D)] for x in range(D)])  # Walsh-Hadamard

def phase(A, Q):
    """A: set of monomials for the +-1 part, Q: dict variable -> {0..3} for an i^(sum q_j x_j) linear Z4 part."""
    g = np.zeros(D, int)
    for m in A: g ^= np.prod(X[:, list(m)], axis=1)
    h = np.zeros(D, int)
    for j, q in Q.items(): h = (h + q * X[:, j]) % 4
    return ((-1) ** g) * (1j ** h)

def supports(state):
    u = [phase(*s) for s in state] + [np.ones(D)]
    T = {}
    for a in range(4):
        for b in range(4):
            if a == b: continue
            f = np.conj(u[a]) * u[b]
            T[a, b] = np.abs(H @ f) > 1e-9          # Walsh transform support (indexed by h - g = h xor g)
    return T

def k4_count(T):
    # count (x2, x3, x4) with x1 = 0 forming a K4 (differences are XOR in Z_2^n)
    c = 0
    S12 = np.nonzero(T[0, 1])[0]; S13 = np.nonzero(T[0, 2])[0]; S14 = T[0, 3]
    for x2 in S12:
        for x3 in S13:
            if not T[1, 2][x2 ^ x3]: continue
            ok = S14 & T[1, 3][np.arange(D) ^ x2] & T[2, 3][np.arange(D) ^ x3]
            c += int(ok.sum())
    return c

def rand_state():
    return [(set(rng.sample(monos, rng.randint(1, 4))), {}) for _ in range(3)]

state = rand_state(); T = supports(state); best = cur = k4_count(T)
print("n", n, "D", D, "initial K4s through vertex", cur, flush=True)
for it in range(iters):
    new = [(set(A), dict(Q)) for A, Q in state]
    a = rng.randrange(3)
    if rng.random() < 0.8:
        m = rng.choice(monos); new[a][0].symmetric_difference_update({m})
    else:
        j = rng.randrange(n); new[a][1][j] = (new[a][1].get(j, 0) + 1) % 4
    T2 = supports(new); c = k4_count(T2)
    if c <= cur or rng.random() < 0.01:
        state, cur = new, c
        if c < best:
            best = c; degs = [int(T2[0, 3].sum()), int(T2[1, 3].sum()), int(T2[2, 3].sum())]
            print("it", it, "K4s", c, "|T_a4|", degs, "state", [(sorted(A), Q) for A, Q in state], flush=True)
            if c == 0: print("K4-FREE: four-context KS set found", flush=True); break
print("best", best)

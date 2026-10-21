# Four Pauli (Lagrangian) eigenbases never form a KS set: numerical check for three qubits.
# Theorem (RESULTS.md 3.7): the Pauli group preserves every Lagrangian eigenbasis and acts transitively on it, so
# by the averaging lemma the union is KS only if it is K4-free; K4-freeness is excluded by the commuting-span
# obstruction.  Here: (1) the Mermin star (5 bases) is KS but none of its 4-subsets is; (2) random 4-tuples of the
# 135 three-qubit Lagrangians (biased towards large intersections) are all KS-colourable, the uniform mixture of
# their K4s has uniform marginals, and no K4 ever misses a 2-piece internal basis.
import sys, itertools, random, numpy as np
from fourbasis import ortho_graph, all_bases, colouring, k4s, pieces

n = 3; D = 2 ** n
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]]); Zm = np.diag([1, -1])
def pauli(x, z):   # Hermitian representative i^{x.z} X^x Z^z
    op = np.eye(1)
    for a, b in zip(x, z):
        op = np.kron(op, np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Zm, b))
    return (1j) ** (sum(a * b for a, b in zip(x, z)) % 4) * op
vecs6 = [tuple(v) for v in itertools.product((0, 1), repeat=2 * n) if any(v)]
omega = lambda u, v: (sum(u[i] * v[n + i] + u[n + i] * v[i] for i in range(n))) % 2
add = lambda u, v: tuple((a + b) % 2 for a, b in zip(u, v))
def span(gens):
    S = {tuple([0] * 2 * n)}
    for g in gens: S |= {add(s, g) for s in S}
    return frozenset(S)
lag = set()
for g1, g2, g3 in itertools.combinations(vecs6, 3):
    if omega(g1, g2) or omega(g1, g3) or omega(g2, g3): continue
    S = span([g1, g2, g3])
    if len(S) == D: lag.add(S)
lag = sorted(lag, key=sorted); print("Lagrangians:", len(lag))
def generators(L):
    gens = []
    for v in sorted(L):
        if any(v) and v not in span(gens): gens.append(v)
    return gens
rng = np.random.default_rng(1)
def eigenbasis(L):
    P = [pauli(g[:n], g[n:]) for g in generators(L)]; out = []
    for s in itertools.product((1, -1), repeat=n):
        Pi = np.eye(D)
        for si, p in zip(s, P): Pi = Pi @ (np.eye(D) + si * p) / 2
        v = Pi @ (rng.normal(size=D) + 1j * rng.normal(size=D)); out.append(v / np.linalg.norm(v))
    return out
EB = {L: eigenbasis(L) for L in lag}

def config(Ls):
    vecs, cover = [], []
    for L in Ls:
        cover.append(list(range(len(vecs), len(vecs) + D))); vecs += EB[L]
    G = ortho_graph(vecs); B = all_bases(G, D)
    return vecs, cover, G, B

def P(s):   # e.g. "XYY" -> its Z2^6 vector
    x = tuple(int(c in "XY") for c in s); z = tuple(int(c in "ZY") for c in s); return x + z
mermin = [span([P("XII"), P("IXI"), P("IIX")]), span([P("XII"), P("IYI"), P("IIY")]),
          span([P("YII"), P("IXI"), P("IIY")]), span([P("YII"), P("IYI"), P("IIX")]),
          span([P("XXX"), P("XYY"), P("YXY")])]
assert all(L in EB for L in mermin)
def check(Ls):
    vecs, cover, G, B = config(Ls)
    assert len(vecs) == 4 * D and all(len(set(c)) == D for c in cover)
    T = k4s(G, cover); assert T, "K4-free four Pauli bases would contradict d >= 17"
    cnt = np.zeros(len(vecs))
    for t in T:
        for r in t: cnt[r] += 1
    uniform = all(np.ptp(cnt[c]) == 0 for c in cover)
    internal = [b for b in B if b not in [frozenset(c) for c in cover]]
    two = [b for b in internal if pieces(b, cover) == 2]; more = [b for b in internal if pieces(b, cover) >= 3]
    assert all(any(r in b for r in t) for b in two for t in T)        # 2-piece lemma
    blocked = sum(1 for t in T if any(not any(r in b for r in t) for b in more))
    return colouring(G, B) is not None, uniform, len(internal), len(more), blocked, len(T)

if __name__ == "__main__":
    vecs, cover, G, B = config(mermin)
    print("Mermin star: rays", len(vecs), "bases", len(B), "KS:", colouring(G, B) is None)
    for sub in itertools.combinations(mermin, 4):
        vecs, cover, G, B = config(sub); f = colouring(G, B)
        print("  4 of the 5 bases: rays", len(vecs), "bases", len(B), "KS:", f is None)
    random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    tally = {"colourable": 0, "uniform K4 marginals": 0, "with >=3-piece internal basis": 0, "some K4 blocked": 0}
    for it in range(N):
        if it % 2 == 0: Ls = random.sample(lag, 4)
        else:      # walk: each next Lagrangian shares a non-identity Pauli with the previous one
            Ls = [random.choice(lag)]
            while len(Ls) < 4:
                Ls.append(random.choice([L for L in lag if L not in Ls and len(L & Ls[-1]) >= 2]))
        col, uni, ni, n3, nb, nT = check(Ls)
        tally["colourable"] += col; tally["uniform K4 marginals"] += uni
        tally["with >=3-piece internal basis"] += n3 > 0; tally["some K4 blocked"] += nb > 0
    print("random 4-tuples:", N, tally)

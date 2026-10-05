# Local search for a KS set covered by four orthonormal bases, over a finite ray alphabet.
# State: 4 bases from the alphabet.  V = union.  Contexts = ALL alphabet bases contained in V.
# Score = number of KS colourings of V (one 1 per context, no two orthogonal 1s), capped.
import sys, pickle, random, itertools
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType

rays, G, bases = pickle.load(open(sys.argv[1], "rb"))
seed = int(sys.argv[2]); iters = int(sys.argv[3]); CAP = 200
rng = random.Random(seed)
bset = [frozenset(b) for b in bases]
ray_bases = {}
for i, b in enumerate(bset):
    for r in b: ray_bases.setdefault(r, []).append(i)
adj = {v: set(G[v]) for v in G}

def contexts(V):
    cand = set()
    for r in V: cand.update(ray_bases.get(r, []))
    return [bset[i] for i in cand if bset[i] <= V]

def count_colourings(V, cap=CAP):
    V = sorted(V); var = {v: i + 1 for i, v in enumerate(V)}; top = len(V)
    s = Cadical153()
    ctx = contexts(set(V))
    for c in ctx:
        lits = [var[v] for v in c]
        s.add_clause(lits)
        for a, b in itertools.combinations(lits, 2): s.add_clause([-a, -b])
    for a, b in itertools.combinations(V, 2):
        if b in adj[a]: s.add_clause([-var[a], -var[b]])
    n = 0
    while n < cap and s.solve():
        m = s.get_model(); n += 1
        s.add_clause([-l for l in m if l > 0 and l <= top] or [0])
    s.delete()
    return n, len(ctx)

# neighbour bases: bases sharing at least one ray with current union (to create internal bases)
def random_basis_touching(V):
    r = rng.choice(sorted(V)); return bset[rng.choice(ray_bases[r])]

def key(V, B):
    if len(set(B)) < 4: return (10**9, 0, 0)
    c, nc = count_colourings(V)
    return (-nc, c, nc)          # lexicographic: more internal contexts first, then fewer colourings
B = [bset[rng.randrange(len(bset))]]
while len(set(B)) < 4: B.append(random_basis_touching(set().union(*B))); B = list(dict.fromkeys(B))
V = set().union(*B); kk = key(V, B); score = kk; best = (kk, B)
T = 2.0
for it in range(iters):
    k = rng.randrange(4); others = [b for j, b in enumerate(B) if j != k]
    U = set().union(*others)
    nb = random_basis_touching(U) if rng.random() < 0.9 else bset[rng.randrange(len(bset))]
    B2 = others[:k] + [nb] + others[k:]
    V2 = set().union(*B2); k2 = key(V2, B2)
    if k2 <= score or rng.random() < 0.02:
        B, score = B2, k2
        if score < best[0]:
            best = (score, B); print("it", it, "internal contexts", score[2], "colourings(cap 200)", score[1], "rays", len(V2), flush=True)
            if score[1] == 0:
                print("KS SET FOUND", [sorted(b) for b in B], flush=True)
                pickle.dump(B, open(f"ks4_found_{seed}.pkl", "wb")); break
    T = max(0.05, T * 0.9995)
print("best", best[0])

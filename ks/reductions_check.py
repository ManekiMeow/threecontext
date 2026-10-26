# Checks for the reduction lemmas of RESULTS.md 3.7.
# (a) Deletion/monotonicity: rays outside every basis are removable; adding rays keeps KS.
# (b) The 2-piece lemma is sharp: a 3-piece internal basis CAN block a K4 (d = 3 example).
# (c) Peres' 33-ray set: which rays lie in no orthogonal triad (they are removable by (a)).
import itertools, numpy as np
from fourbasis import ortho_graph, all_bases, colouring, k4s, pieces
from pauli4 import config, mermin

rng = np.random.default_rng(7)
# (a) Mermin star (KS, d = 8) plus extra rays lying in no basis
vecs, cover, G, B = config(mermin)
extra = []
for _ in range(5):                                   # generic rays: orthogonal to nothing
    v = rng.normal(size=8) + 1j * rng.normal(size=8); extra.append(v / np.linalg.norm(v))
for a, b in [(0, 9), (3, 17)]:                       # rays orthogonal to a few Mermin rays, completing no basis
    M = np.array([vecs[a], vecs[b]]).conj(); _, _, Vh = np.linalg.svd(M)
    v = Vh[2:].conj().T @ (rng.normal(size=6) + 1j * rng.normal(size=6)); extra.append(v / np.linalg.norm(v))
V2 = vecs + extra; G2 = ortho_graph(V2); B2 = all_bases(G2, 8)
outside = [r for r in range(len(V2)) if not any(r in b for b in B2)]
print("(a) Mermin star + 7 extra rays: bases", len(B2), "rays in no basis", len(outside),
      "KS:", colouring(G2, B2) is None)
keep = [r for r in range(len(V2)) if r not in outside]; V3 = [V2[r] for r in keep]
G3 = ortho_graph(V3); B3 = all_bases(G3, 8)
print("    after deleting the rays in no basis: rays", len(V3), "KS:", colouring(G3, B3) is None)

# (b) d = 3: internal basis {e1, e2, e3} with one ray in each of B1, B2, B3 blocks the K4 (a1, a2, a3, a4)
s = 1 / np.sqrt(2); e = np.eye(3)
B1 = [e[0], s * (e[1] + e[2]), s * (e[1] - e[2])]
B2 = [e[1], s * (e[0] + e[2]), s * (e[0] - e[2])]
B3 = [e[2], s * (e[0] + e[1]), s * (e[0] - e[1])]
U, _ = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))); B4 = list(U.T)
vecs = B1 + B2 + B3 + B4; cover = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11]]
G = ortho_graph(vecs); B = all_bases(G, 3)
internal = [b for b in B if b not in [frozenset(c) for c in cover]]
T = k4s(G, cover); t = (1, 4, 7, next(r for r in cover[3] if all(not G.has_edge(r, x) for x in (1, 4, 7))))
print("(b) internal bases and their pieces:", [(sorted(b), pieces(b, cover)) for b in internal])
print("    K4", t, "is a K4:", t in T, " misses the 3-piece basis {0,3,6}:", not set(t) & {0, 3, 6})

# (c) Peres' 33 rays in R^3: permutations/sign changes of (0,0,1), (0,1,1), (0,1,sqrt2), (1,1,sqrt2)
r2 = np.sqrt(2); raw = set()
for base in [(0, 0, 1), (0, 1, 1), (0, 1, r2), (1, 1, r2)]:
    for p in itertools.permutations(base):
        for sg in itertools.product((1, -1), repeat=3):
            v = np.array(p) * np.array(sg); v = v / np.linalg.norm(v)
            k = next(i for i in range(3) if abs(v[i]) > 1e-9); v = v * np.sign(v[k])
            raw.add(tuple(np.round(v, 12)))
P = [np.array(v) for v in sorted(raw)]; GP = ortho_graph(P); BP = all_bases(GP, 3)
noT = [i for i in range(len(P)) if not any(i in b for b in BP)]
print("(c) Peres: rays", len(P), "orthogonal pairs", GP.number_of_edges(), "triads", len(BP),
      "KS:", colouring(GP, BP) is None, " rays in no triad:", len(noT))
if noT:
    Q = [P[i] for i in range(len(P)) if i not in noT]; GQ = ortho_graph(Q)
    print("    after deleting them: rays", len(Q), "KS:", colouring(GQ, all_bases(GQ, 3)) is None)

# (d) Critical cores can have incomplete covering cliques: the 21-ray set (d = 6) is critical, and its clique cover
# number is 5 with two triangles among the cliques (ks21.py, ks21_cover5.py).
from ks21 import bases as B21
w = np.exp(2j * np.pi / 3); val = {(0, 0): 0, (1, 0): 1, (0, 1): w, (-1, -1): w * w}
R21 = sorted({v for b in B21 for v in b})
V21 = [np.array([val[x] for x in v]) / np.linalg.norm([val[x] for x in v]) for v in R21]
G21 = ortho_graph(V21)
crit = all(colouring(ortho_graph(V21[:r] + V21[r + 1:]), all_bases(ortho_graph(V21[:r] + V21[r + 1:]), 6))
           is not None for r in range(len(V21)))
print("(d) 21-ray set: KS:", colouring(G21, all_bases(G21, 6)) is None, " critical:", crit)

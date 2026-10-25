# Checks for RESULTS.md 3.10 (can the orthogonality graph be strongly regular / regular / vertex-transitive?).
# G = orthogonality graph of four disjoint bases of C^d (n = 4d), H = complement (non-orthogonality graph).
# (A) SRG parameter analysis: if G is an SRG then every basis is a Delsarte clique, so every ray outside a basis is
#     orthogonal to the same number e of its rays; k = d-1+3e, r = d-1-e, s = -k/(d-1).  Enumerate the integral
#     feasible cases and confirm: s in {-2,-3}; s = -2 only for d in {4, 7, 16}; s = -3 always has lambda_H = 2.
# (B) Delsarte-Hoffman intersection on concrete SRGs: rook graph L2(4) with a Latin-square partition into 4
#     cocliques (every K4 of H meets every 4-coclique once) and Shrikhande (K4-free, local graphs C6).
# (C) Regular-graph spectral conditions lambda_min(H) <= -D/3 and lambda_2(H) >= (3d-D)/(d-1) on realized examples.
import itertools, math, numpy as np, networkx as nx
from fractions import Fraction as Fr

print("(A) SRG parameters with a spread of four Delsarte d-cliques (d <= 3000)")
fam = {}
for d in range(3, 3001):
    n = 4 * d
    for e in range(1, d - 1):                       # e <= d-2 (disjoint bases)
        k = d - 1 + 3 * e; r = d - 1 - e; s = Fr(-k, d - 1)
        if s.denominator != 1: continue
        s = int(s); mu = k + r * s; lam = mu + r + s
        f = Fr((n - 1) * (-s) - k, r - s); g = Fr(k + (n - 1) * r, r - s)
        if f.denominator != 1 or g.denominator != 1 or f < 0 or g < 0 or mu <= 0: continue
        lamH = n - 2 * k + mu - 2; muH = n - 2 * k + lam
        fam.setdefault(s, []).append((d, e, (n, k, lam, mu), (n, n - 1 - k, lamH, muH)))
for s, rows in sorted(fam.items()):
    print(f"  s = {s}: {len(rows)} parameter sets; d values {[r[0] for r in rows]}")
    print(f"           lambda_H values {sorted({r[3][2] for r in rows})}; examples G, H = {rows[:2]}")
assert set(fam) <= {-2, -3}
assert [r[0] for r in fam[-2]] == [4, 7, 16] and {r[3][2] for r in fam[-3]} == {2}

print("(B) theta of complements of cycles (local graphs when lambda_H = 2):",
      {m: round(1 + 1 / math.cos(math.pi / m), 4) if m % 2 else 2.0 for m in range(3, 10)}, "(m = 3 gives 3)")
def srg_params(Gr):
    A = nx.to_numpy_array(Gr, nodelist=sorted(Gr)); n = len(A); k = A.sum(1)
    if np.ptp(k): return None
    A2 = A @ A; lam = {int(A2[i, j]) for i in range(n) for j in range(n) if i != j and A[i, j]}
    mu = {int(A2[i, j]) for i in range(n) for j in range(n) if i != j and not A[i, j]}
    return (n, int(k[0]), lam, mu)
def report(name, H, bases):
    G = nx.complement(H); K4 = [c for c in nx.find_cliques(H) if len(c) == 4]
    cocl = [set(c) for c in nx.find_cliques(G) if len(c) == len(bases[0])]
    inter = {len(set(T) & S) for T in K4 for S in cocl}
    loc = sorted({tuple(sorted(len(c) for c in nx.cycle_basis(H.subgraph(H[w])))) for w in H})
    print(f"  {name}: H srg {srg_params(H)}, G srg {srg_params(G)}; K4s in H {len(K4)}; d-cliques of G {len(cocl)};"
          f" |K4 ∩ d-clique| values {inter or '-'}; local graph cycle lengths {loc}")
Z4 = range(4)
rook = nx.Graph([((i, j), (a, b)) for i, j, a, b in itertools.product(Z4, Z4, Z4, Z4)
                 if (i, j) < (a, b) and (i == a or j == b)])
report("rook L2(4), bases {(i,i+k)}", rook, [[(i, (i + k) % 4) for i in Z4] for k in Z4])
S = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
shr = nx.Graph([((i, j), (a, b)) for i, j, a, b in itertools.product(Z4, Z4, Z4, Z4)
                if (i, j) < (a, b) and ((a - i) % 4, (b - j) % 4) in S])
cos = [[((x + t) % 4, (y + 2 * t) % 4) for t in Z4] for (x, y) in [(0, 0), (0, 1), (0, 2), (0, 3)]]
assert all(not shr.has_edge(u, v) for b in cos for u, v in itertools.combinations(b, 2))
report("Shrikhande, bases = cosets of <(1,2)>", shr, cos)

print("(C) regular-graph conditions on realized configurations")
from pauli4 import config, mermin
for name, Ls in [("4 Mermin product bases (d=8)", mermin[:4]), ("Mermin product+GHZ bases (d=8)", mermin[1:])]:
    vecs, cover, Gq, B = config(Ls); Hq = nx.complement(Gq); d = 8
    A = nx.to_numpy_array(Hq, nodelist=range(len(vecs))); deg = A.sum(1); ev = np.sort(np.linalg.eigvalsh(A))
    D = deg[0]
    print(f"  {name}: regular {np.ptp(deg) == 0} D = {D:.0f}; lambda_min(H) = {ev[0]:.3f} <= -D/3 = {-D/3:.3f};"
          f" lambda_2(H) = {ev[-2]:.3f} >= (3d-D)/(d-1) = {(3*d-D)/(d-1):.3f}")

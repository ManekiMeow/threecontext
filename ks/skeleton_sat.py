# Combinatorial skeleton of a Route-A (K4-free) four-basis KS set in C^d, as a SAT instance.
# Variables: e[k][l][i][j] = 1 iff ray i of basis k is non-orthogonal to ray j of basis l (k < l).
# Necessary conditions encoded:
#   (K4)  no transversal 4-clique;
#   (C2)  a ~ b  =>  a and b have a common neighbour in each other basis   (<a|b> = sum_c <a|c><c|b>);
#   (C4)  a ⊥ b  =>  a and b have 0 or >= 2 common neighbours in each other basis (same reason);
#         this applies to rays of the same basis and to non-adjacent rays of different bases;
#   (DEG) every ray has >= PMIN (= 3: a GHZ graph has no class of size <= 2) neighbours in each other basis and >= NMIN neighbours in total (Section 2: n >= 28);
#   (TRI) every triple of bases contains a triangle (Re tr(U12 U23 U31) = d > 0).
# UNSAT for some d proves that no K4-free four-basis KS set exists in that dimension.
import sys, itertools, time
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

d = int(sys.argv[1]); NMIN = int(sys.argv[2]) if len(sys.argv) > 2 else 28; PMIN = int(sys.argv[3]) if len(sys.argv) > 3 else 3
pool = IDPool(); cls = []
def E(k, i, l, j):
    if k > l: k, l, i, j = l, k, j, i
    return pool.id(("e", k, i, l, j))
parts = range(4); R = range(d)
# aux t(k,i,l,j,m,c) <-> e(k,i,m,c) & e(l,j,m,c): "c in basis m is a common neighbour of (k,i) and (l,j)"
def T(k, i, l, j, m, c):
    key = ("t", k, i, l, j, m, c); new = key not in pool.obj2id
    v = pool.id(key)
    if new:
        a, b = E(k, i, m, c), E(l, j, m, c)
        cls.extend([[-v, a], [-v, b], [v, -a, -b]])
    return v
# (K4)
for i, j, m, n in itertools.product(R, repeat=4):
    cls.append([-E(0, i, 1, j), -E(0, i, 2, m), -E(0, i, 3, n), -E(1, j, 2, m), -E(1, j, 3, n), -E(2, m, 3, n)])
# (C2) and (C4) for pairs of rays in different bases
for k, l in itertools.combinations(parts, 2):
    for m in parts:
        if m in (k, l): continue
        for i, j in itertools.product(R, R):
            ts = [T(k, i, l, j, m, c) for c in R]
            cls.append([-E(k, i, l, j)] + ts)                                  # adjacent => >= 1 common nb
            for c in R:                                                        # non-adjacent & c common => another one
                cls.append([E(k, i, l, j), -ts[c]] + [t for t in ts if t != ts[c]])
# (C4) for pairs of rays in the same basis (always orthogonal)
for k in parts:
    for m in parts:
        if m == k: continue
        for i, j in itertools.combinations(R, 2):
            ts = [T(k, i, k, j, m, c) for c in R]
            for c in R:
                cls.append([-ts[c]] + [t for t in ts if t != ts[c]])
# (DEG)
for k in parts:
    for i in R:
        tot = []
        for l in parts:
            if l == k: continue
            lits = [E(k, i, l, j) for j in R]; tot += lits
            cls.extend(CardEnc.atleast(lits, bound=PMIN, vpool=pool, encoding=EncType.seqcounter).clauses)
        cls.extend(CardEnc.atleast(tot, bound=NMIN, vpool=pool, encoding=EncType.seqcounter).clauses)
# (SYM) symmetry breaking: rays within a basis are interchangeable.  Sound scheme: among all equivalent
# configurations take the lex-largest (M12, then M23, then M34); it has non-increasing rows and columns of M12
# (double-lex) and non-increasing columns of M23 and of M34.
def lex_geq(X, Y):
    """clauses for X >=_lex Y (lists of literals of equal length)"""
    prev = pool.id(("lex", id(X), id(Y), -1)); cls.append([prev])
    for i, (x, y) in enumerate(zip(X, Y)):
        e = pool.id(("lex", id(X), id(Y), i))
        cls.extend([[-prev, x, -y], [-e, prev], [-e, -x, y], [-e, x, -y], [-prev, -x, -y, e], [-prev, x, y, e]])
        prev = e
rows12 = [[E(0, i, 1, j) for j in R] for i in R]; cols12 = [[E(0, i, 1, j) for i in R] for j in R]
cols23 = [[E(1, j, 2, c) for j in R] for c in R]; cols34 = [[E(2, c, 3, n) for c in R] for n in R]
if "--nosym" not in sys.argv:
    for seq in (rows12, cols12, cols23, cols34):
        for i in range(d - 1): lex_geq(seq[i], seq[i + 1])
# (TRI)
for k, l, m in itertools.combinations(parts, 3):
    tri = []
    for i, j, c in itertools.product(R, R, R):
        v = pool.id(("tri", k, i, l, j, m, c)); tri.append(v)
        cls.extend([[-v, E(k, i, l, j)], [-v, E(k, i, m, c)], [-v, E(l, j, m, c)]])
    cls.append(tri)
print(f"d={d} NMIN={NMIN}: {pool.top} vars, {len(cls)} clauses", flush=True)
t0 = time.time(); s = Cadical153(bootstrap_with=cls)
res = s.solve()
print(f"d={d}: {'SAT' if res else 'UNSAT'} in {time.time()-t0:.1f}s", flush=True)
if res:
    model = set(v for v in s.get_model() if v > 0)
    degs = [[sum(E(k, i, l, j) in model for l in parts if l != k for j in R) for i in R] for k in parts]
    print("degree range per basis:", [(min(x), max(x)) for x in degs])

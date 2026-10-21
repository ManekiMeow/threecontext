# SAT relaxation of the four-context KS problem in small dimension d (Route B; d <= 16).
# Variables: o(r,s) = "rays r and s (in different covering bases) are orthogonal".  Every clause below is a
# NECESSARY condition for four pairwise disjoint orthonormal bases of C^d whose union is KS (RESULTS.md 3.7-3.8):
#   KS     every K4 (pairwise non-orthogonal transversal) misses an internal basis meeting >= 3 covering bases
#   K4     every ray lies in a K4 (d <= 16: otherwise its neighbourhood is a 3-context GHZ graph)
#   deg    a ray is orthogonal to at most d-2 rays of another basis (disjointness)
#   supp   x not-perp y => N_l(x), N_l(y) meet;  x perp y => |N_l(x) & N_l(y)| != 1   (expand <x|y> in basis B_l)
#   compl  two distinct rays are never orthogonal to a common (d-1)-clique
#   span   internal basis B' with piece A_k (|A_k| >= 2): B_k\A_k and B'\A_k span the same space, so a ray is
#          orthogonal to all of one iff it is orthogonal to all of the other
#   split  no proper subspace is spanned by subsets of all four bases (Lemma 10; no example in dim <= 4)
# plus lex-leader symmetry breaking for adjacent transpositions inside each basis.
# UNSAT would prove that no four-context KS set exists in C^d.  A model is only a candidate graph.
import sys, itertools, time
from pysat.solvers import Cadical153

d = int(sys.argv[1]) if len(sys.argv) > 1 else 5
OFF = {a[5:] for a in sys.argv if a.startswith('--no-')}   # e.g. --no-span disables a constraint group
NB = 4; R = NB * d; bas = lambda r: r // d
B = [list(range(k * d, (k + 1) * d)) for k in range(NB)]
TRUE, FALSE = "T", "F"
nv = 0
def new():
    global nv; nv += 1; return nv
ovar = {}
for r, s in itertools.combinations(range(R), 2):
    if bas(r) != bas(s): ovar[(r, s)] = new()
def o(r, s):
    if r > s: r, s = s, r
    return TRUE if bas(r) == bas(s) else ovar[(r, s)]
def neg(l): return FALSE if l == TRUE else TRUE if l == FALSE else -l
clauses = []
def add(lits):
    out = []
    for l in lits:
        if l == TRUE: return
        if l != FALSE: out.append(l)
    clauses.append(out)
def AND(lits):            # aux variable equivalent to the conjunction of lits
    lits = [l for l in lits if l != TRUE]
    if FALSE in lits: return FALSE
    if not lits: return TRUE
    if len(lits) == 1: return lits[0]
    a = new()
    for l in lits: add([-a, l])
    add([a] + [neg(l) for l in lits]); return a

# internal bases: q[S] <-> S is a d-clique
q = {}
for S in itertools.combinations(range(R), d):
    q[S] = AND([o(r, s) for r, s in itertools.combinations(S, 2)])
npieces = lambda S: len({bas(r) for r in S})
# K4s
k4 = {}
for T in itertools.product(*B):
    k4[T] = AND([neg(o(r, s)) for r, s in itertools.combinations(T, 2)])
t0 = time.time()
for T, kt in (k4.items() if 'ks' not in OFF else []):                        # KS
    Ts = set(T)
    add([neg(kt)] + [q[S] for S in q if npieces(S) >= 3 and not Ts & set(S)])
for r in (range(R) if 'k4' not in OFF else []):                              # K4
    add([k4[T] for T in k4 if r in T])
for r in range(R):                                                           # deg
    for l in range(NB):
        if l != bas(r):
            for F in itertools.combinations(B[l], d - 1): add([neg(o(r, s)) for s in F])
for l in (range(NB) if 'supp' not in OFF else []):                           # supp
    out = [r for r in range(R) if bas(r) != l]
    for x, y in itertools.combinations(out, 2):
        c = {z: AND([neg(o(x, z)), neg(o(y, z))]) for z in B[l]}
        add([o(x, y)] + list(c.values()))
        for z in B[l]: add([neg(o(x, y)), neg(c[z])] + [c[w] for w in B[l] if w != z])
for K in (itertools.combinations(range(R), d - 1) if 'compl' not in OFF else []):  # compl
    rest = [r for r in range(R) if r not in K]
    for x, y in itertools.combinations(rest, 2):
        add([neg(q[tuple(sorted(K + (x,)))]), neg(q[tuple(sorted(K + (y,)))])])
for S, qs in (q.items() if 'span' not in OFF else []):                       # span
    if qs == TRUE or npieces(S) < 2: continue
    for k in {bas(r) for r in S}:
        A = [r for r in S if bas(r) == k]
        if len(A) < 2: continue
        X1 = [r for r in B[k] if r not in A]; X2 = [r for r in S if bas(r) != k]
        for z in range(R):
            if z in X1 or z in X2 or z in A: continue
            for Y, Yp in ((X1, X2), (X2, X1)):
                for yp in Yp: add([neg(qs)] + [neg(o(z, y)) for y in Y] + [o(z, yp)])
for m in (range(2, min(d // 2, 4) + 1) if 'split' not in OFF else []):       # split (parts of dim <= 4)
    for X in itertools.product(*[list(itertools.combinations(b, m)) for b in B]):
        add([neg(o(x, y)) for k, l in itertools.combinations(range(NB), 2)
             for x in X[k] for y in B[l] if y not in X[l]])
def lexleq(pairs):   # (v_1..v_m) <=_lex (w_1..w_m)
    e = TRUE
    for v, w in pairs:
        if v == w: continue
        add([neg(e), neg(v), w])
        e2 = new(); add([neg(e), v, w, e2]); add([neg(e), neg(v), neg(w), e2]); e = e2
order = sorted(ovar, key=lambda p: ovar[p])
def transposition(a, b): return lambda r: b if r == a else a if r == b else r
CASE = sys.argv[sys.argv.index("--case") + 1] if "--case" in sys.argv else None
if CASE:
    # WLOG one K4 is the last ray of every basis, and an internal basis blocking it has piece sizes CASE
    # (e.g. "311": 3 rays of B_0, 1 of B_1, 1 of B_2), made of the first rays of those bases.
    sizes = [int(c) for c in CASE] + [0] * (NB - len(CASE)); assert sum(sizes) == d
    K4 = [B[k][-1] for k in range(NB)]; blk = [r for k in range(NB) for r in B[k][:sizes[k]]]
    for r, s in itertools.combinations(K4, 2): add([neg(o(r, s))])
    for r, s in itertools.combinations(blk, 2): add([o(r, s)])
    classes = [B[k][:sizes[k]] for k in range(NB)] + [B[k][sizes[k]:-1] for k in range(NB)]
    gens = [transposition(c[i], c[i + 1]) for c in classes for i in range(len(c) - 1)]
    for k, l in itertools.combinations(range(NB), 2):        # swap two bases playing the same role
        if sizes[k] == sizes[l]:
            gens.append(lambda r, k=k, l=l: r + (l - k) * d if bas(r) == k else r - (l - k) * d if bas(r) == l else r)
else:
    gens = [transposition(B[k][i], B[k][i + 1]) for k in range(NB) for i in range(d - 1)] if 'sym' not in OFF else []
    if "--swap" in sys.argv:                                                 # also break basis permutations
        gens += [lambda r, k=k: r + d if bas(r) == k else r - d if bas(r) == k + 1 else r for k in range(NB - 1)]
for sw in gens:                                                              # symmetry breaking (lex-leader)
    lexleq([(o(r, s), o(sw(r), sw(s))) for r, s in order])
print(f"d={d}: vars {nv}, clauses {len(clauses)}, built in {time.time() - t0:.1f}s", flush=True)
if "--fix" in sys.argv:          # soundness test: fix the o-variables to a given orthogonality pattern
    import pickle; fixed = set(pickle.load(open(sys.argv[sys.argv.index("--fix") + 1], "rb")))
    clauses += [[v] if p in fixed else [-v] for p, v in ovar.items()]
import pysat.solvers as ps
Solver = getattr(ps, sys.argv[sys.argv.index("--solver") + 1]) if "--solver" in sys.argv else Cadical153
s = Solver(bootstrap_with=clauses); t0 = time.time()
sat = s.solve(); print("SAT" if sat else "UNSAT", f"({time.time() - t0:.1f}s)", flush=True)
if sat:
    m = set(l for l in s.get_model() if l > 0)
    orth = sorted(p for p, v in ovar.items() if v in m)
    import pickle; pickle.dump(orth, open(f"smalld_d{d}_model.pkl", "wb"))
    print("orthogonal cross pairs:", len(orth))
    print("internal bases:", [S for S, v in q.items() if v != TRUE and v in m])

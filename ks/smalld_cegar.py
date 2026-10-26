# Second, independent encoding of the small-d four-context KS conditions (cross-check of smalld_sat.py, and
# scalable towards d = 6).  Differences from smalld_sat.py:
#  * "every K4 is blocked" uses selection variables: for each transversal T, if T is a K4 then some d rays outside
#    T, meeting >= 3 covering bases, are pairwise orthogonal (so there are no clique variables for all d-subsets);
#  * unique completion ("compl") and span equalities ("span") are added lazily: solve, find violated instances in
#    the model by brute force, add exactly those clauses, and repeat (CEGAR).
# Static clauses: deg, supp, split, "every ray lies in a K4", and the WLOG case split with residual lex-leader.
# Usage: python3 smalld_cegar.py d CASE [solver]      e.g.  python3 smalld_cegar.py 5 311 Cadical153
import sys, itertools, time
import pysat.solvers as ps
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool

d = int(sys.argv[1]); CASE = sys.argv[2]; SOLVER = sys.argv[3] if len(sys.argv) > 3 else "Cadical153"
NB = 4; R = NB * d; bas = lambda r: r // d; B = [list(range(k * d, (k + 1) * d)) for k in range(NB)]
pool = IDPool(); cl = []
def o(r, s):            # literal "r perp s" for rays in different bases (None = constant true, same basis)
    if bas(r) == bas(s): return None
    return pool.id(("o", min(r, s), max(r, s)))
def perp_clause(pairs, extra=()):   # clause: NOT(all pairs orthogonal) OR extra
    lits = [-o(r, s) for r, s in pairs if o(r, s) is not None]
    return lits + list(extra)
def AND(name, lits):
    a = pool.id(name)
    for l in lits: cl.append([-a, l])
    cl.append([a] + [-l for l in lits]); return a
T_all = list(itertools.product(*B))
k4 = {T: AND(("k4", T), [-o(r, s) for r, s in itertools.combinations(T, 2)]) for T in T_all}
for T in T_all:                                                   # blocked K4s (selection encoding)
    kt = k4[T]; out = [r for r in range(R) if r not in T]; sel = {r: pool.id(("sel", T, r)) for r in out}
    card = CardEnc.equals(lits=list(sel.values()), bound=d, vpool=pool, encoding=EncType.seqcounter)
    for c in card.clauses: cl.append([-kt] + c)
    for r, s in itertools.combinations(out, 2):
        if o(r, s) is not None: cl.append([-kt, -sel[r], -sel[s], o(r, s)])
    p = {k: pool.id(("piece", T, k)) for k in range(NB)}
    for k in range(NB): cl.append([-kt, -p[k]] + [sel[r] for r in B[k] if r in sel])
    for k, l in itertools.combinations(range(NB), 2): cl.append([-kt, p[k], p[l]])   # >= 3 pieces
for r in range(R): cl.append([k4[T] for T in T_all if r in T])     # every ray lies in a K4 (d <= 16)
for r in range(R):                                                 # deg
    for l in range(NB):
        if l != bas(r):
            for F in itertools.combinations(B[l], d - 1): cl.append([-o(r, s) for s in F])
for l in range(NB):                                                # supp
    out = [r for r in range(R) if bas(r) != l]
    for x, y in itertools.combinations(out, 2):
        c = {z: AND(("c", x, y, z), [-o(x, z), -o(y, z)]) for z in B[l]}
        oxy = o(x, y)
        if oxy is not None: cl.append([oxy] + list(c.values()))
        for z in B[l]: cl.append(([-oxy] if oxy is not None else []) + [-c[z]] + [c[w] for w in B[l] if w != z])
for m in range(2, min(d // 2, 4) + 1):                             # split
    for X in itertools.product(*[list(itertools.combinations(b, m)) for b in B]):
        cl.append([-o(x, y) for k, l in itertools.combinations(range(NB), 2) for x in X[k] for y in B[l] if y not in X[l]])
sizes = [int(c) for c in CASE] + [0] * (NB - len(CASE)); assert sum(sizes) == d and len(CASE) >= 3
K4 = [b[-1] for b in B]; blk = [r for k in range(NB) for r in B[k][:sizes[k]]]
for r, s in itertools.combinations(K4, 2): cl.append([-o(r, s)])
for r, s in itertools.combinations(blk, 2):
    if o(r, s) is not None: cl.append([o(r, s)])
order = sorted(((r, s) for r, s in itertools.combinations(range(R), 2) if bas(r) != bas(s)), key=lambda p: o(*p))
LEX = [0]
def lexleq(perm):
    LEX[0] += 1; e = None
    for r, s in order:
        v, w = o(r, s), o(perm(r), perm(s))
        if v == w: continue
        pre = [] if e is None else [-e]
        cl.append(pre + [-v, w]); e2 = pool.id(("lex", LEX[0], r, s))
        cl.append(pre + [v, w, e2]); cl.append(pre + [-v, -w, e2]); e = e2
classes = [B[k][:sizes[k]] for k in range(NB)] + [B[k][sizes[k]:-1] for k in range(NB)]
for c in classes:
    for i in range(len(c) - 1): lexleq(lambda r, a=c[i], b=c[i + 1]: b if r == a else a if r == b else r)
for k, l in itertools.combinations(range(NB), 2):
    if sizes[k] == sizes[l]: lexleq(lambda r, k=k, l=l: r + (l - k) * d if bas(r) == k else r - (l - k) * d if bas(r) == l else r)

def violations(orth):
    perp = lambda r, s: bas(r) == bas(s) or (min(r, s), max(r, s)) in orth
    cliques = [S for S in itertools.combinations(range(R), d) if all(perp(a, b) for a, b in itertools.combinations(S, 2))]
    new = []
    for K in itertools.combinations(range(R), d - 1):              # compl
        if not all(perp(a, b) for a, b in itertools.combinations(K, 2)): continue
        xs = [x for x in range(R) if x not in K and all(perp(x, k) for k in K)]
        for x, y in itertools.combinations(xs, 2):
            new.append(perp_clause(list(itertools.combinations(K, 2)) + [(x, k) for k in K] + [(y, k) for k in K]))
    covering = [tuple(b) for b in B]
    for S in cliques:                                              # span
        if S in covering: continue
        Sp = list(itertools.combinations(S, 2))
        for k in {bas(r) for r in S}:
            A = [r for r in S if bas(r) == k]
            if len(A) < 2: continue
            X1 = [r for r in B[k] if r not in A]; X2 = [r for r in S if bas(r) != k]
            for z in range(R):
                if z in X1 or z in X2 or z in A: continue
                for Y, Yp in ((X1, X2), (X2, X1)):
                    if all(perp(z, y) for y in Y):
                        for yp in Yp:
                            if not perp(z, yp):
                                new.append(perp_clause(Sp + [(z, y) for y in Y], [o(z, yp)] if o(z, yp) else []))
    return new

t0 = time.time(); s = getattr(ps, SOLVER)(bootstrap_with=cl)
print(f"d={d} case {CASE}: vars {pool.top}, static clauses {len(cl)}", flush=True)
rounds = 0
while True:
    if not s.solve():
        print(f"UNSAT after {rounds} refinement rounds ({time.time() - t0:.1f}s)", flush=True); break
    m = set(l for l in s.get_model() if l > 0)
    orth = {(r, q) for r, q in itertools.combinations(range(R), 2) if bas(r) != bas(q) and o(r, q) in m}
    new = violations(orth); rounds += 1
    if not new:
        import pickle; pickle.dump(sorted(orth), open(f"cegar_d{d}_{CASE}_model.pkl", "wb"))
        print(f"SAT: model satisfies all constraints after {rounds} rounds ({time.time() - t0:.1f}s)", flush=True); break
    for c in new: s.add_clause(c)
    if rounds % 50 == 0: print(f"  round {rounds}: +{len(new)} clauses ({time.time() - t0:.0f}s)", flush=True)

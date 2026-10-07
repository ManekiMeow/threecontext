"""Decide whether a ray universe contains a KS set of size <= T: CEGAR with an incremental
SAT master (CaDiCaL + totalizer cardinality) and a SAT colouring oracle.

Same cut logic as minks.py: every 0/1 assignment c of the universe gives a valid clause
"S contains an edge with both ends 1 under c, or a basis with all 0 under c".
"""
import random
import sys
import time

from pysat.card import CardEnc, EncType, ITotalizer
from pysat.formula import IDPool
from pysat.solvers import Cadical153

from minks import Instance, log, sym_actions


def decide(inst, T, sym=None, per_iter=32, fix=(), verbose=True, seed=0, perms=()):
    random.seed(seed)
    pool = IDPool()
    s = [pool.id(("s", i)) for i in range(inst.n)]
    y = [pool.id(("y", k)) for k in range(len(inst.edges))]
    z = [pool.id(("z", k)) for k in range(len(inst.tris))]
    cl = []
    for k, (i, j) in enumerate(inst.edges):
        cl += [[-y[k], s[i]], [-y[k], s[j]]]
    for k, t in enumerate(inst.tris):
        cl += [[-z[k], s[i]] for i in t]
    for i in range(inst.n):
        # criticality: in a basis inside S, and >= 3 orthogonal rays inside S
        cl.append([-s[i]] + [z[inst.tid[t]] for t in inst.vtris[i]])
        nb = [s[j] for j in inst.nbr[i]]
        if len(nb) < 3:
            cl.append([-s[i]])
        else:
            enc = CardEnc.atleast(nb, bound=3, vpool=pool, encoding=EncType.seqcounter)
            cl += [[-s[i]] + c for c in enc.clauses]
    # Lex-leader symmetry breaking: S >=lex g(S) for each ray permutation g in `perms`.
    # Sound because the set of feasible S (all cuts are added together with their images) is
    # invariant under the group, so each orbit keeps its lex-largest member.
    for gi, img in enumerate(perms):
        inv = [0] * inst.n
        for v, gv in enumerate(img):
            inv[gv] = v
        e = pool.id(("lex", gi, 0))
        cl.append([e])
        for i in range(inst.n):
            x, yv = s[i], s[inv[i]]
            if x == yv:
                continue
            cl.append([-e, -yv, x])
            e2 = pool.id(("lex", gi, i + 1))
            cl.append([-e, -x, -yv, e2])
            cl.append([-e, x, yv, e2])
            e = e2
    tot = ITotalizer(lits=s, ubound=T + 1, top_id=pool.top)
    pool.occupy(pool.top + 1, tot.top_id)
    cl += tot.cnf.clauses
    cl.append([-tot.rhs[T]])  # at most T selected
    cl += [[s[i]] for i in fix]
    master = Cadical153(bootstrap_with=cl)
    seen = set()
    ncuts = 0

    def add_cut(c):
        nonlocal ncuts
        we, wt = inst.witnesses(c)
        key = (tuple(we), tuple(wt))
        if key in seen:
            return
        seen.add(key)
        master.add_clause([y[k] for k in we] + [z[k] for k in wt])
        ncuts += 1

    for _ in range(50):
        c = inst.improve([random.randint(0, 1) for _ in range(inst.n)], set())
        for g in (sym or [None]):
            add_cut(c if g is None else g(c))
    t0 = time.time()
    it = 0
    while True:
        it += 1
        if not master.solve():
            if verbose:
                log(f"iter {it}: no KS set of size <= {T} (fix={list(fix)}), "
                    f"{ncuts} cuts, {time.time()-t0:.1f}s")
            master.delete()
            return None
        mod = master.get_model()
        S = [i for i in range(inst.n) if mod[s[i] - 1] > 0]
        cols = []
        for r in range(per_iter):
            c = inst.colour(S, seed=r if r else None)
            if c is None:
                break
            cols.append(c)
        if not cols:
            log(f"iter {it}: FOUND uncolourable set of size {len(S)}")
            master.delete()
            return S
        for c in cols:
            c = inst.improve(c, set(S))
            for g in (sym or [None]):
                add_cut(c if g is None else g(c))
        if verbose and it % 20 == 0:
            log(f"iter {it}: |S|={len(S)}, cuts {ncuts}, {time.time()-t0:.1f}s")


def ray_perms(inst, group=None):
    from rays import signed_perms
    idx = {r: i for i, r in enumerate(inst.rays)}
    out = []
    for g in (group or signed_perms()):
        img = [idx.get(g(r)) for r in inst.rays]
        if None not in img and img != list(range(inst.n)):
            out.append(img)
    return out


if __name__ == "__main__":
    from rays import integer_rays, prune_to_triangles
    k, T = int(sys.argv[1]), int(sys.argv[2])
    rays = prune_to_triangles(integer_rays(k))
    inst = Instance(rays)
    log(f"height {k}: {inst.n} rays, {len(inst.edges)} pairs, {len(inst.tris)} bases; target <= {T}")
    S = decide(inst, T, sym=sym_actions(inst), perms=ray_perms(inst))
    if S:
        log("KS set:", [rays[i] for i in S])

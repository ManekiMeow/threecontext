"""Minimum Kochen-Specker subset of a finite ray universe in R^3, by implicit hitting sets.

KS colouring (010-colouring) of a ray set S: c: S -> {0,1}, no two orthogonal rays both 1,
every orthogonal basis (triangle) inside S has a ray coloured 1. S is a KS set iff it has none.

Key fact: for ANY 0/1 assignment c of the whole universe V, a KS set S must contain a witness
against c (an edge with both ends 1, or a triangle with all 0). So every assignment gives a valid
cut; assignments with few witnesses give strong cuts. Master (CP-SAT): minimum |S| hitting all
cuts. Oracle (CaDiCaL): is S colourable? If yes its colouring, extended to V greedily, is a new cut.
When the master's optimum S is uncolourable, |S| is the minimum KS set size inside V (exact).
"""
import random
import sys
import time

from ortools.sat.python import cp_model
from pysat.solvers import Cadical153

from rays import integer_rays, prune_to_triangles, structure, signed_perms


def log(*a):
    print(*a, flush=True)


class Instance:
    def __init__(self, rays, edges=None, tris=None):
        self.rays = rays
        self.n = len(rays)
        if edges is None:
            edges, tris = structure(rays)
        self.edges, self.tris = edges, tris
        self.nbr = [[] for _ in range(self.n)]
        for i, j in self.edges:
            self.nbr[i].append(j)
            self.nbr[j].append(i)
        self.vtris = [[] for _ in range(self.n)]
        for t in self.tris:
            for i in t:
                self.vtris[i].append(t)
        self.eid = {e: k for k, e in enumerate(self.edges)}
        self.tid = {t: k for k, t in enumerate(self.tris)}

    def colour(self, S, seed=None):
        """A KS colouring of S (dict) or None. Constraints only among rays of S."""
        Sset = set(S)
        cl = []
        for i, j in self.edges:
            if i in Sset and j in Sset:
                cl.append([-(i + 1), -(j + 1)])
        for t in self.tris:
            if all(x in Sset for x in t):
                cl.append([x + 1 for x in t])
        with Cadical153(bootstrap_with=cl) as s:
            if seed is not None:
                rnd = random.Random(seed)
                s.set_phases([(i + 1) * rnd.choice((1, -1)) for i in S])
            if not s.solve():
                return None
            m = set(l for l in s.get_model() if l > 0)
        return [1 if (i + 1) in m else 0 for i in range(self.n)]

    def witnesses(self, c):
        we = [self.eid[(i, j)] for i, j in self.edges if c[i] and c[j]]
        wt = [self.tid[t] for t in self.tris if not (c[t[0]] or c[t[1]] or c[t[2]])]
        return we, wt

    def improve(self, c, frozen):
        """Greedy local search on rays outside `frozen` to reduce the number of witnesses."""
        c = list(c)

        def cost(v):
            if c[v]:
                return sum(c[u] for u in self.nbr[v])
            return sum(1 for t in self.vtris[v] if not (c[t[0]] or c[t[1]] or c[t[2]]))

        free = [v for v in range(self.n) if v not in frozen]
        changed = True
        while changed:
            changed = False
            random.shuffle(free)
            for v in free:
                before = cost(v)
                c[v] ^= 1
                if cost(v) < before:
                    changed = True
                else:
                    c[v] ^= 1
        return c


def solve(inst, lb=0, per_iter=8, sym=None, time_limit=None, verbose=True, target=None):
    """Minimum KS set in inst; with target=T, only decide whether one of size <= T exists."""
    m = cp_model.CpModel()
    s = [m.NewBoolVar(f"s{i}") for i in range(inst.n)]
    y = [m.NewBoolVar(f"y{k}") for k in range(len(inst.edges))]
    z = [m.NewBoolVar(f"z{k}") for k in range(len(inst.tris))]
    for k, (i, j) in enumerate(inst.edges):
        m.AddImplication(y[k], s[i])
        m.AddImplication(y[k], s[j])
    for k, t in enumerate(inst.tris):
        for i in t:
            m.AddImplication(z[k], s[i])
    # Any minimum KS set is vertex-critical, hence (i) every ray lies in a basis inside S and
    # (ii) every ray has >= 3 orthogonal rays inside S (a ray of degree <= 2 can always be
    # coloured after colouring the rest). So these constraints lose no minimum solution.
    for i in range(inst.n):
        m.AddBoolOr([z[inst.tid[t]] for t in inst.vtris[i]]).OnlyEnforceIf(s[i])
        m.Add(sum(s[j] for j in inst.nbr[i]) >= 3).OnlyEnforceIf(s[i])
    m.Add(sum(s) >= lb)
    if target is None:
        m.Minimize(sum(s))
    else:
        m.Add(sum(s) <= target)
    seen = set()
    ncuts = 0

    def add_cut(c):
        nonlocal ncuts
        we, wt = inst.witnesses(c)
        key = (tuple(we), tuple(wt))
        if key in seen:
            return
        seen.add(key)
        m.AddBoolOr([y[k] for k in we] + [z[k] for k in wt])
        ncuts += 1

    # seed cuts: locally optimised random assignments
    for _ in range(50):
        c = inst.improve([random.randint(0, 1) for _ in range(inst.n)], set())
        for g in (sym or [None]):
            add_cut(c if g is None else g(c))

    t0 = time.time()
    it = 0
    while True:
        it += 1
        solver = cp_model.CpSolver()
        solver.parameters.num_workers = 4
        if time_limit:
            solver.parameters.max_time_in_seconds = max(1.0, time_limit - (time.time() - t0))
        st = solver.Solve(m)
        if target is not None and st == cp_model.INFEASIBLE:
            log(f"iter {it}: no KS set of size <= {target} in this universe "
                f"({ncuts} cuts, {time.time()-t0:.1f}s)")
            return None, target + 1
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE) or (target is None and st != cp_model.OPTIMAL):
            log(f"master status {solver.StatusName(st)} after {it} iterations")
            return None, None
        S = [i for i in range(inst.n) if solver.Value(s[i])]
        bound = len(S)
        if target is None:
            m.Add(sum(s) >= bound)
        m.ClearHints()
        for i in range(inst.n):
            m.AddHint(s[i], solver.Value(s[i]))
        cols = []
        for r in range(per_iter):
            c = inst.colour(S, seed=r if r else None)
            if c is None:
                break
            cols.append(c)
        if not cols:
            if verbose:
                log(f"iter {it}: |S|={bound} UNCOLOURABLE{'' if target else ' -> minimum = '+str(bound)}  "
                    f"({ncuts} cuts, {time.time()-t0:.1f}s)")
            return S, bound
        for c in cols:
            c = inst.improve(c, set(S))
            for g in (sym or [None]):
                add_cut(c if g is None else g(c))
        if verbose:
            log(f"iter {it}: {'candidate size' if target else 'lower bound'} {bound}, "
                f"cuts {ncuts}, {time.time()-t0:.1f}s")


def sym_actions(inst, group=None):
    """Assignment-level actions of the signed permutation group (only those preserving V)."""
    idx = {r: i for i, r in enumerate(inst.rays)}
    acts = []
    for g in (group or signed_perms()):
        img = [idx.get(g(r)) for r in inst.rays]
        if None in img:
            continue
        acts.append(lambda c, img=img: _perm(c, img))
    return acts


def _perm(c, img):
    out = [0] * len(c)
    for v, gv in enumerate(img):
        out[gv] = c[v]
    return out


if __name__ == "__main__":
    k = int(sys.argv[1])
    target = int(sys.argv[2]) if len(sys.argv) > 2 else None
    random.seed(0)
    rays = prune_to_triangles(integer_rays(k))
    inst = Instance(rays)
    log(f"height {k}: {inst.n} rays, {len(inst.edges)} orthogonal pairs, {len(inst.tris)} bases")
    S, b = solve(inst, sym=sym_actions(inst), target=target, per_iter=32)
    if S:
        log("KS set:", [rays[i] for i in S])

# Soundness test for the geometric clauses of smalld_sat.py (deg, supp, compl, span): they must be satisfiable when
# the o-variables are fixed to the orthogonality pattern of ANY four pairwise disjoint bases of C^d.
# Configurations: random 4-tuples of disjoint bases drawn from ray alphabets with many orthogonalities, hill-climbed
# towards many internal bases (--three: internal bases meeting >= 3 covering bases), and designed configurations
# containing a blocking basis of each d = 5 shape.  Usage: python3 smalld_soundness.py d N MININT [--three]
import sys, itertools, pickle, random, subprocess, numpy as np, os
from fourbasis import ortho_graph, all_bases

d = int(sys.argv[1]) if len(sys.argv) > 1 else 5; N = int(sys.argv[2]) if len(sys.argv) > 2 else 6
MININT = int(sys.argv[3]) if len(sys.argv) > 3 else 0     # only test configurations with this many internal bases
here = os.path.dirname(os.path.abspath(__file__)); random.seed(3)
w = np.exp(2j * np.pi / 3)
alphabets = {"{0,±1}": [0, 1, -1], "{0,1,w,w²}": [0, 1, w, w * w]}
for name, A in (alphabets.items() if "--designed-only" not in sys.argv else []):
    rays = {}
    for v in itertools.product(A, repeat=d):
        v = np.array(v, dtype=complex)
        if not v.any(): continue
        v = v / v[np.flatnonzero(v)[0]]; v = v / np.linalg.norm(v); rays[tuple(np.round(v, 9))] = v
    rays = list(rays.values()); G = ortho_graph(rays); Bs = all_bases(G, d)
    print(name, "rays", len(rays), "bases", len(Bs), flush=True)
    def internal(pick):   # number of internal bases (with --three: only those meeting >= 3 covering bases)
        Gs = ortho_graph([rays[r] for b in pick for r in sorted(b)]); Bi = all_bases(Gs, d)
        if "--three" in sys.argv:
            return Gs, sum(1 for b in Bi if sum(1 for c in range(4) if b & set(range(c * d, (c + 1) * d))) >= 3)
        return Gs, len(Bi) - 4
    done = 0; tries = 0
    while done < N and tries < 50:
        tries += 1; pick = []
        for b in random.sample(Bs, len(Bs)):
            if all(not (b & c) for c in pick): pick.append(b)
            if len(pick) == 4: break
        if len(pick) < 4: continue
        Gs, nint = internal(pick)
        for _ in range(300):        # hill-climb towards many internal bases (these exercise the span clauses)
            if nint >= MININT: break
            k = random.randrange(4); b = random.choice(Bs)
            if any(b & c for j, c in enumerate(pick) if j != k): continue
            new = pick[:k] + [b] + pick[k + 1:]; G2, n2 = internal(new)
            if n2 >= nint: pick, Gs, nint = new, G2, n2
        if nint < MININT: continue
        orth = sorted((r, s) for r, s in Gs.edges() if r // d != s // d)
        pickle.dump(orth, open("/tmp/_pattern.pkl", "wb"))
        out = subprocess.run([sys.executable, os.path.join(here, "smalld_sat.py"), str(d), "--no-ks", "--no-k4",
                              "--no-split", "--no-sym", "--fix", "/tmp/_pattern.pkl"],
                             capture_output=True, text=True, cwd="/tmp").stdout
        verdict = "SAT" if "\nSAT" in out else "UNSAT (UNSOUND!)"
        print(f"  config {done}: cross orthogonalities {len(orth)}, internal bases {nint}: {verdict}", flush=True)
        done += 1

# Designed configurations: each contains an internal basis of a prescribed shape meeting >= 3 covering bases
# (W_1 + ... + W_m = C^d with W_k spanned by part of B_k), which exercises the span and completion clauses.
rng = np.random.default_rng(5)
def haar(n):
    Q, R = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))); return Q * (np.diag(R) / abs(np.diag(R)))
def designed(shape):
    F = haar(d); cols, start = [], 0
    for n in shape: cols.append(F[:, start:start + n]); start += n
    bases = []
    for k in range(4):
        if k < len(shape):
            Wk = cols[k] @ haar(shape[k]); comp = np.hstack([c for j, c in enumerate(cols) if j != k])
            bases.append(np.hstack([Wk, comp @ haar(d - shape[k])]))
        else:
            bases.append(haar(d))
    return [b[:, i] for b in bases for i in range(d)]
for shape in [(3, 1, 1), (2, 2, 1), (2, 1, 1, 1)][: 3 if d == 5 else 0]:
    for t in range(2):
        vecs = designed(shape); Gs = ortho_graph(vecs); Bs = all_bases(Gs, d)
        orth = sorted((r, s) for r, s in Gs.edges() if r // d != s // d)
        pickle.dump(orth, open("/tmp/_pattern.pkl", "wb"))
        out = subprocess.run([sys.executable, os.path.join(here, "smalld_sat.py"), str(d), "--no-ks", "--no-k4",
                              "--no-split", "--no-sym", "--fix", "/tmp/_pattern.pkl"],
                             capture_output=True, text=True, cwd="/tmp").stdout
        verdict = "SAT" if "\nSAT" in out else "UNSAT (UNSOUND!)"
        pc = sorted(sum(1 for c in range(4) if b & set(range(c * d, (c + 1) * d))) for b in Bs)
        print(f"  designed {shape} #{t}: cross orthogonalities {len(orth)}, pieces of bases {pc}: {verdict}", flush=True)

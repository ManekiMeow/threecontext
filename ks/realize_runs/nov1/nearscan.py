import os
import sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))  # ks/
import warnings; warnings.filterwarnings("ignore")
import multiplier_classify as mc
e = int(sys.argv[1])
for p in map(int, sys.argv[2:]):
    t0 = time.time()
    pats = mc.k4_free_patterns_fast(p, e)
    best = []
    for pat in pats:
        why, th = mc.full_check(p, e, pat)
        if th is None or min(th) > 3 - 1e-6:
            continue
        best.append((min(th), sum(t > 3 - 1e-6 for t in th), pat, [round(float(t), 4) for t in th]))
    best.sort(reverse=True)
    print(f"p = {p}: {len(pats)} K4-free, {len(best)} near misses; top: {best[:2]} ({time.time()-t0:.0f} s)", flush=True)

"""full_check for every K4-free single-class pattern of a classify log (h >= 9 only: (DEG) needs |S_kl| >= 9)."""
import re, ast, sys, time, collections
sys.path.insert(0, "/home/user/threecontext/ks")
import warnings; warnings.filterwarnings("ignore")
import multiplier_classify as m
e, logf = int(sys.argv[1]), sys.argv[2]
for line in open(logf):
    mm = re.match(r"p=(\d+) h=(\d+): (\d+) K4-free single-class patterns \(S_01 = C_0\) (\[.*\]) \(", line)
    if not mm or int(mm.group(3)) == 0 or int(mm.group(2)) < 9:
        continue
    p, pats = int(mm.group(1)), ast.literal_eval(mm.group(4))
    t0 = time.time()
    stats, passing = collections.Counter(), []
    for pat in pats:
        why, th = m.full_check(p, e, pat)
        stats[why] += 1
        if why is None:
            passing.append((pat, th))
    print(f"p={p} e={e}: {len(pats)} K4-free, {dict(stats)}, PASSING {passing} ({time.time() - t0:.0f} s)", flush=True)

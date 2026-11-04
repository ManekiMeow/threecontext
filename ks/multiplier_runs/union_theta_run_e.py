import re, ast, sys, time
sys.path.insert(0, "/home/user/threecontext/ks")
import warnings; warnings.filterwarnings("ignore")
import multiplier_classify as m
E = int(sys.argv[1]); LOG = sys.argv[2]
DONE = set()
if len(sys.argv) > 3:
    for l in open(sys.argv[3]):
        if ': thetas' in l:
            DONE.add(l.split(': thetas')[0])
SC = "/tmp/claude-0/-home-user-threecontext/5e407493-7848-5a52-a504-a0152de6dc76/scratchpad/night"
cur = None
for line in open(LOG):
    mm = re.match(r"p=(\d+):", line)
    if mm:
        cur = int(mm.group(1)); continue
    if line.startswith("    ("):
        T, z = ast.literal_eval(line.strip())
        if all(len(t) == 1 for t in T) and not any(z):
            continue                                   # single-class: done by full_check
        t0 = time.time()
        key = f"p={cur} T={T} z={z}"
        if key in DONE:
            continue
        try:
            th = m.union_thetas(cur, E, T, z)
        except Exception as ex:
            print(f"{key}: ERROR {type(ex).__name__}", flush=True)
            continue
        print(f"p={cur} T={T} z={z}: thetas {th} {'PASS' if min(th) > 3 - 1e-6 and len(th) == 4 else 'fail'} "
              f"({time.time() - t0:.1f} s)", flush=True)

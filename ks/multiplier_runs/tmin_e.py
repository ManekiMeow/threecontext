"""t_min (smallest number of common neighbours, in one other basis, of a non-orthogonal pair) for every single-class
candidate of a single_full log; t_min < sqrt(h) rules out flat realizations."""
import re, ast, sys, itertools, math
sys.path.insert(0, "/home/user/threecontext/ks")
import cyclo4
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
e, logf = int(sys.argv[1]), sys.argv[2]
tally = {"no flat": 0, "flat possible": 0}
for line in open(logf):
    m = re.match(r"p=(\d+) e=\d: .*PASSING (\[.*\]) \(", line)
    if not m or m.group(2) == "[]":
        continue
    p = int(m.group(1)); C = cyclo4.classes(p, e); h = (p - 1) // e
    for pat, th in ast.literal_eval(m.group(2).replace("np.float64", "")):
        S = {}
        for (k, l), c in zip(PAIRS, pat):
            S[(k, l)] = set(int(x) for x in C[c]); S[(l, k)] = {(-x) % p for x in C[c]}
        tmin = min(len({c for c in S[(k, b)] if (c - s) % p in S[(l, b)]})
                   for k in range(4) for l in range(4) if l != k for s in S[(k, l)]
                   for b in range(4) if b not in (k, l))
        key = "no flat" if tmin < math.sqrt(h) else "flat possible"
        tally[key] += 1
        print(f"p={p} pattern {pat}: h={h} t_min={tmin} sqrt(h)={math.sqrt(h):.2f} -> {key}", flush=True)
print("TOTAL (patterns, not classes):", tally)

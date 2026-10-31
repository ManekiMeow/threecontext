"""Isomorphism classes (basis permutations x multiplying by the generator) of the PASS lines of union-theta logs."""
import re, ast, sys, itertools, collections
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
def canon(T, z, p, e):
    neg = (p - 1) // 2 % e
    best = None
    for perm in itertools.permutations(range(4)):
        new = {}
        for (k, l), Tk, zk in zip(PAIRS, T, z):
            a, b = perm[k], perm[l]
            if a < b:
                new[(a, b)] = (tuple(sorted(Tk)), zk)
            else:
                new[(b, a)] = (tuple(sorted((x + neg) % e for x in Tk)), zk)
        for r in range(e):
            key = tuple((tuple(sorted((x + r) % e for x in new[pr][0])), new[pr][1]) for pr in PAIRS)
            best = key if best is None or key < best else best
    return best
e = int(sys.argv[1])
classes = collections.defaultdict(set)
for f in sys.argv[2:]:
    for line in open(f):
        m = re.match(r"p=(\d+) T=(.*) z=(.*): thetas .* PASS", line)
        if m:
            p = int(m.group(1)); T = ast.literal_eval(m.group(2)); z = ast.literal_eval(m.group(3))
            classes[p].add(canon(T, z, p, e))
for p in sorted(classes):
    print(f"p={p}: {len(classes[p])} isomorphism classes of union-type candidates: {sorted(classes[p])}")

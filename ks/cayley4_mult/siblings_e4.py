"""(Copied from the session scratch directory; run from ks/.) Siblings of the p = 137 candidate: every S_kl a single cyclotomic class of index e (no 0), all primes p in a range.
Up to multiplying by a generator. Prints patterns passing (K4), (DEG), (C2), (C4), (TRI) and theta = 3 at all bases."""
import sys, itertools, time
sys.path.insert(0, "/home/user/threecontext/ks")
import numpy as np
import cayley4 as c4, cyclo4
PAIRS = c4.PAIRS
e, plo, phi = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
for p in range(plo, phi + 1):
    if p < 3 or any(p % q == 0 for q in range(2, int(p ** 0.5) + 1)) or (p - 1) % e:
        continue
    t = time.time()
    C = cyclo4.classes(p, e)
    vec = [np.isin(np.arange(p), c) for c in C]
    enc = c4.Cayley4.__new__(c4.Cayley4); enc.G, enc.d = cyclo4.Zp(p), p
    seen, stats, found = set(), {}, []
    for combo in itertools.product(range(e), repeat=6):
        canon = min(tuple((i + r) % e for i in combo) for r in range(e))     # multiply by the generator
        if canon in seen:
            continue
        seen.add(canon)
        S = {pr: vec[i] for pr, i in zip(PAIRS, combo)}
        why = cyclo4.check(p, S)
        stats[why] = stats.get(why, 0) + 1
        if why is None:
            sets = cyclo4.to_sets(p, S)
            th = []
            for k in (3, 0, 1, 2):                                    # stop at the first failing basis
                V, A = enc.neighbourhood(sets, k)
                th.append(2.0 if c4.two_colouring(A) is not None else c4.theta_complement(V, A))
                if th[-1] < 3 - 1e-6:
                    break
            if min(th) > 3 - 1e-6:
                found.append(combo)
            else:
                stats["theta<3"] = stats.get("theta<3", 0) + 1
    print(f"p={p} e={e} h={(p-1)//e}: {len(seen)} patterns, {stats}, passing {found} ({time.time()-t:.0f} s)", flush=True)

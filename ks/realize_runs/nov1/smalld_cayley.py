import os
import sys, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))  # ks/
import warnings; warnings.filterwarnings("ignore")
import cayley4 as c4, groups
limit = float(sys.argv[1])
for spec in sys.argv[2:]:
    n, idx = spec.split(":") if ":" in spec else (spec, None)
    gs = groups.groups_of_order(int(n))
    for i, G in enumerate(gs):
        if idx is not None and str(i) not in idx.split(","):
            continue
        t0 = time.time()
        st = c4.search(G, limit, log=lambda s: print("   ", s, flush=True), progress=300, max_candidates=5, realize_restarts=1)
        cands = [{k: v for k, v in c.items() if k != "sets"} | {"sets": c["sets"]} for c in st["candidates"]]
        print(f"RESULT {G.name} (d={G.n}): {st['status']} after {time.time()-t0:.0f} s; models {st['models']}, "
              f"bipartite {st['bipartite_blocks']}, theta blocks {st['theta_blocks']} (max theta {st.get('max_theta')}), "
              f"candidates {len(cands)}", flush=True)
        for c in cands:
            print("CANDIDATE", json.dumps(c), flush=True)

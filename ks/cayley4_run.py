"""Overnight driver for ks/cayley4.py (RESULTS.md §3.12).

Two families of K4-free four-partite graphs, both searched with the same SAT + theta-CEGAR engine:
  * "mult":  Cayley graphs over Z_p whose connection sets are unions of orbits of the e-th powers
             (multiplier-invariant); one job per (p, e). Small instances, usually decided.
  * "group": all Cayley graphs over a group Γ (every group in ks/groups.py of the given orders);
             one job per group, with a time budget and resumable CEGAR state.

usage: python3 cayley4_run.py --stop-at 2026-10-29T06:45:00Z [--mult-primes 19-199] [--mult-limit 600]
                              [--orders 18-30] [--per-group 2700] [--workers 4]

Per-job state lives in cayley4_runs/; decided jobs (UNSAT, realizable candidate, budget used up)
are skipped on the next call. One JSON line per job-run goes to cayley4_results.jsonl and a
readable line to cayley4.log.
"""
import argparse
import datetime
import json
import os
import re
import time
from multiprocessing import Pool

import cayley4 as c4
import groups

HERE = os.path.dirname(os.path.abspath(__file__))
STATE = os.path.join(HERE, "cayley4_runs")
OUT = HERE


def state_path(key):
    return os.path.join(STATE, re.sub(r"[^A-Za-z0-9]+", "_", key) + ".json")


def load_state(key):
    p = state_path(key)
    if os.path.exists(p):
        return json.load(open(p))
    return dict(key=key, seconds=0.0, runs=0, status=None, models=0, bipartite_blocks=0, theta_blocks=0,
                candidates=[], blocks=[])


def utcnow():
    return datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%S")


def _set_state(path):
    global STATE
    STATE = path


def is_prime(n):
    return n > 1 and all(n % q for q in range(2, int(n ** 0.5) + 1))


def job(args):
    kind, a, b, budget, stop_at = args
    key = f"mult_p{a}_e{b}" if kind == "mult" else f"group_{a}"
    st = load_state(key)
    if st["status"] in ("UNSAT", "realizable", "budget"):
        return None
    limit = min(budget - st["seconds"], stop_at - time.time() - 30)
    if limit < 30:
        return None
    if kind == "mult":
        G, orbits = groups.cyclic(a), c4.multiplier_orbits(a, b)
    else:
        G, orbits = next(g for g in groups.groups_of_order(b) if g.name == a), None
    lines = []
    res = c4.search(G, limit, log=lines.append, blocks=st["blocks"], progress=900, orbits=orbits)
    st["seconds"] += res["seconds"]
    st["runs"] += 1
    for k in ("models", "bipartite_blocks", "theta_blocks"):
        st[k] += res[k]
    st["candidates"] += res["candidates"]
    if res["status"] == "UNSAT":
        st["status"] = "UNSAT"
    elif any(c["realize_residual"] < 1e-9 for c in res["candidates"]):
        st["status"] = "realizable"
    elif st["seconds"] >= budget - 60:
        st["status"] = "budget"
    else:
        st["status"] = "paused"
    tmp = state_path(key) + ".tmp"
    json.dump(st, open(tmp, "w"))
    os.replace(tmp, state_path(key))
    summary = {k: v for k, v in st.items() if k not in ("blocks", "candidates")}
    summary.update(kind=kind, group=G.name, d=G.n, orbits=len(orbits) if orbits else G.n, vars=res["vars"],
                   clauses=res["clauses"], candidates=len(st["candidates"]),
                   best_residual=min((c["realize_residual"] for c in st["candidates"]), default=None),
                   learned=len(st["blocks"]), finished=utcnow())
    return summary, lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stop-at", required=True, help="UTC, e.g. 2026-10-29T06:45:00Z")
    ap.add_argument("--mult-primes", default="19-199")
    ap.add_argument("--mult-limit", type=float, default=600)
    ap.add_argument("--orders", default="18-30")
    ap.add_argument("--per-group", type=float, default=2700)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default=HERE, help="directory for cayley4_runs/, cayley4.log, cayley4_results.jsonl")
    a = ap.parse_args()
    global STATE
    STATE = os.path.join(a.out, "cayley4_runs")
    stop_at = datetime.datetime.fromisoformat(a.stop_at.replace("Z", "+00:00")).timestamp()
    os.makedirs(STATE, exist_ok=True)
    plo, phi = map(int, a.mult_primes.split("-"))
    jobs = [("mult", p, e, a.mult_limit, stop_at) for p in range(plo, phi + 1) if is_prime(p)
            for e in range(2, p) if (p - 1) % e == 0 and (p - 1) // e >= 2]
    lo, hi = map(int, a.orders.split("-"))
    jobs += [("group", G.name, n, a.per_group, stop_at) for n in range(lo, hi + 1) for G in groups.groups_of_order(n)]
    with open(os.path.join(a.out, "cayley4.log"), "a") as log, \
            open(os.path.join(a.out, "cayley4_results.jsonl"), "a") as out, \
            Pool(a.workers, maxtasksperchild=1, initializer=_set_state, initargs=(STATE,)) as pool:
        log.write(f"# run started {utcnow()}Z: {len(jobs)} jobs, "
                  f"mult primes {a.mult_primes} ({a.mult_limit:.0f} s each), orders {a.orders} "
                  f"({a.per_group:.0f} s each), stop at {a.stop_at}\n")
        log.flush()
        for r in pool.imap_unordered(job, jobs):
            if r is None:
                continue
            s, lines = r
            out.write(json.dumps(s) + "\n")
            out.flush()
            log.write(f"{s['key']:>22s} (d={s['d']}, {s['orbits']} orbits): {s['status']:>8s} after "
                      f"{s['seconds']:.0f} s; models {s['models']}, bipartite blocks {s['bipartite_blocks']}, "
                      f"theta blocks {s['theta_blocks']}, candidates {s['candidates']}, "
                      f"best residual {s['best_residual']}\n")
            for line in lines:
                log.write(line + "\n")
            log.flush()


if __name__ == "__main__":
    main()

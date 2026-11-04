"""Recompute the minimal unembeddable graphs of small order (RESULTS.md §7.4).

A graph is *embeddable* if its vertices can be given pairwise non-collinear vectors in R^3 with
adjacent vertices orthogonal (i.e. it is a subgraph of an orthogonality graph). It is *minimal
unembeddable* (MUG) if it is not embeddable but every proper subgraph is. Li–Bright–Ganesh
(arXiv:2306.13319, Table 3) count 0, 2, 5, 10 squarefree MUGs on <= 9, 10, 11, 12 vertices; every
KS candidate their search rejects contains one. Their list is in mugs.json; this script
reproduces it up to order 10 (orders 11-12 are slow here because Z3 times out on some graphs).

Pipeline for each order n (graphs from nauty: connected, squarefree, minimum degree >= 2):
  1. a graph containing a MUG of smaller order (as a subgraph) is unembeddable and not minimal;
  2. a numerical search (least squares, random restarts) certifies "embeddable";
  3. the rest are decided exactly by Z3 (nonlinear real arithmetic), or by hand (see RESULTS.md);
  4. an unembeddable graph is minimal iff no single-edge deletion is unembeddable.

usage: GENG=/path/to/nauty/geng python3 unembeddable.py [max_order, default 10]
       (writes mugs_recomputed.json)
"""
import json
import os
import subprocess
import sys
from multiprocessing import Pool

import networkx as nx
import numpy as np
import z3
from networkx.algorithms import isomorphism
from scipy.optimize import least_squares

GENG = os.environ.get("GENG", "geng")
HERE = os.path.dirname(os.path.abspath(__file__))


def squarefree_g6(n):
    """Connected squarefree graphs with minimum degree >= 2 on n vertices (nauty geng -cfd2)."""
    out = subprocess.run([GENG, "-cfd2q", str(n)], capture_output=True, text=True, check=True).stdout
    return out.split()


def core2(G):
    """Delete vertices of degree <= 1 until none is left; such a vertex can always be embedded last."""
    G = G.copy()
    while True:
        low = [v for v in G if G.degree(v) <= 1]
        if not low:
            return G
        G.remove_nodes_from(low)


def numeric_embed(G, tries=10, seed=0, sep=1e-3):
    nodes = list(G)
    n = len(nodes)
    if n == 0:
        return True
    idx = {v: i for i, v in enumerate(nodes)}
    E = np.array([(idx[a], idx[b]) for a, b in G.edges()]).reshape(-1, 2)
    P = np.array([(i, j) for i in range(n) for j in range(i + 1, n)]).reshape(-1, 2)
    rng = np.random.default_rng(seed)

    def res(x):
        U = x.reshape(n, 3)
        U = U / np.linalg.norm(U, axis=1, keepdims=True)
        c = np.einsum("ij,ij->i", U[P[:, 0]], U[P[:, 1]])
        r1 = np.einsum("ij,ij->i", U[E[:, 0]], U[E[:, 1]])
        return np.concatenate([r1, 10 * np.maximum(0, np.abs(c) - (1 - sep)), [0.0]])

    for _ in range(tries):
        r = least_squares(res, rng.normal(size=3 * n), xtol=1e-15, ftol=1e-15, gtol=1e-15,
                          max_nfev=3000)
        if np.max(np.abs(r.fun)) < 1e-11:
            return True
    return False


def z3_embed(G, timeout_s=600):
    """Exact decision (Z3, QF_NRA). Returns True, False, or None (timeout)."""
    nodes = list(G)
    n = len(nodes)
    idx = {v: i for i, v in enumerate(nodes)}
    X = [[z3.Real(f"x{i}_{k}") for k in range(3)] for i in range(n)]
    dot = lambda a, b: sum(a[k] * b[k] for k in range(3))
    s = z3.SolverFor("QF_NRA")
    s.set("timeout", timeout_s * 1000)
    b0 = idx[next(iter(G[nodes[0]]))]
    s.add(X[0][0] == 0, X[0][1] == 0, X[0][2] == 1)                  # rotate: v_0 = e3 and
    s.add(X[b0][0] == 1, X[b0][1] == 0, X[b0][2] == 0)               # one neighbour = e1
    for i in range(n):
        s.add(dot(X[i], X[i]) == 1)
    for a, b in G.edges():
        s.add(dot(X[idx[a]], X[idx[b]]) == 0)
    for i in range(n):
        for j in range(i + 1, n):                                   # distinct rays: v_i != +-v_j
            s.add(z3.Or([X[i][k] != X[j][k] for k in range(3)]))
            s.add(z3.Or([X[i][k] != -X[j][k] for k in range(3)]))
    r = s.check()
    return True if r == z3.sat else False if r == z3.unsat else None


def contains(G, H):
    """H is (isomorphic to) a subgraph of G, not necessarily induced."""
    return isomorphism.GraphMatcher(G, H).subgraph_is_monomorphic()


def stage1(args):
    s, small = args
    G = nx.from_graph6_bytes(s.encode())
    if any(contains(G, nx.from_graph6_bytes(h.encode())) for h in small):
        return s, "contains"
    if numeric_embed(G) or numeric_embed(G, tries=100, seed=1):
        return s, "embeddable"
    return s, "open"


def stage2(s):
    r = z3_embed(nx.from_graph6_bytes(s.encode()))
    return s, {True: "embeddable", False: "unembeddable", None: "unknown"}[r]


def mugs_of_order(n, small, procs=4):
    """Returns (number of graphs, status per graph, [(MUG, how it was decided)])."""
    gs = squarefree_g6(n)
    with Pool(procs) as p:
        st = dict(p.map(stage1, [(s, small) for s in gs], chunksize=4))
        open_ = [s for s, v in st.items() if v == "open"]
        st.update(dict(p.map(stage2, open_, chunksize=1)))
    unemb = {s for s, v in st.items() if v in ("unembeddable", "unknown", "contains")}
    mugs = []
    for s in sorted(unemb):
        if st[s] == "contains":
            continue
        G = nx.from_graph6_bytes(s.encode())
        minimal = True
        for a, b in list(G.edges()):
            H = core2(nx.Graph([e for e in G.edges() if e != (a, b)]))
            # a smaller 2-core contains no smaller MUG (else G would), so it is embeddable; an
            # order-n one is again in the list, so look its status up
            if H.number_of_nodes() == n and nx.is_connected(H) and any(
                    nx.is_isomorphic(H, nx.from_graph6_bytes(u.encode())) for u in unemb):
                minimal = False
                break
        if minimal:
            mugs.append((s, st[s]))
    return len(gs), st, mugs


if __name__ == "__main__":
    top = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    found = {}
    small = []
    for n in range(4, top + 1):
        total, st, mugs = mugs_of_order(n, small)
        counts = {k: sum(v == k for v in st.values()) for k in
                  ("embeddable", "contains", "unembeddable", "unknown")}
        print(f"order {n}: {total} graphs; {counts}; minimal unembeddable: {len(mugs)}", flush=True)
        for s, kind in mugs:
            how = "Z3 unsat" if kind == "unembeddable" else "Z3 timeout, see the hand proof in RESULTS.md"
            print(f"   MUG {s}  ({how})", flush=True)
        found[n] = [s for s, _ in mugs]
        small += found[n]
    json.dump(found, open(os.path.join(HERE, "mugs_recomputed.json"), "w"), indent=1)

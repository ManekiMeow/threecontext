"""Bounded SAT/SDP search for uniform-probability (10,10,10) constructions.

The SAT stage enforces triangle-freeness and the six cross-degree pairs allowed
by Lemma14. Independent sets of size11 and disconnected context pairs are cut
lazily. Connectivity is necessary in dimension19, but imposed here as an ansatz
when searching for a uniform example of arbitrary rank. The numerical SDP fixes
all event probabilities to1/10 while allowing distinct real overlap weights.
No numerical failure or time-limited SAT run is a nonexistence certificate.
"""
import argparse
import itertools
import json
from pathlib import Path
import time

import cvxpy as cp
import networkx as nx
import numpy as np
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Cadical153
from scipy.linalg import null_space


ALLOWED = {(2,4), (3,3), (3,4), (4,2), (4,3), (4,4)}


class UniformMargin:
    def __init__(self, parts):
        self.n = sum(map(len, parts))
        C = np.zeros((self.n, 3))
        for k, part in enumerate(parts):
            C[part, k] = 1
        Q = null_space(C.T)
        self.G = cp.Variable((self.n, self.n), symmetric=True)
        self.mask = cp.Parameter((self.n, self.n))
        self.tau = cp.Variable()
        self.problem = cp.Problem(cp.Maximize(self.tau), [
            cp.diag(self.G) == 1, self.G @ C == np.ones((self.n, 3)),
            cp.multiply(self.mask, self.G) == 0,
            Q.T @ self.G @ Q - self.tau*np.eye(self.n-3) >> 0])

    def solve(self, graph):
        mask = np.ones((self.n, self.n))-nx.to_numpy_array(graph, nodelist=range(self.n))
        np.fill_diagonal(mask, 0)
        self.mask.value = mask
        try:
            self.problem.solve(solver="CLARABEL", max_iter=200)
        except cp.error.SolverError:
            # Sparse, degenerate supports can make the dense masked model
            # numerically awkward. This fallback affects only heuristic scores.
            self.problem.solve(solver="SCS", eps=1e-6, max_iters=20000)
        if self.problem.value is None or not np.isfinite(self.problem.value):
            raise RuntimeError(f"uniform margin solve failed: {self.problem.status}")
        return float(self.problem.value)


def build(seed):
    pool, clauses = IDPool(), []
    parts = [list(range(10*k,10*k+10)) for k in range(3)]
    def edge(u, v):
        assert u//10 != v//10
        return pool.id(("e", min(u,v), max(u,v)))
    for a,b,c in itertools.product(*parts):
        clauses.append([-edge(a,b),-edge(b,c),-edge(a,c)])
    for u in range(30):
        selections = []
        for k in range(3):
            if k == u//10:
                continue
            ds = {d: pool.id(("degree",u,k,d)) for d in (2,3,4)}
            selections.append(ds)
            clauses.append(list(ds.values()))
            for d,e in itertools.combinations(ds,2):
                clauses.append([-ds[d],-ds[e]])
            for d in ds:
                enc=CardEnc.equals([edge(u,v) for v in parts[k]],bound=d,
                                  vpool=pool,encoding=EncType.seqcounter)
                clauses.extend([[-ds[d]]+cl for cl in enc.clauses])
        for d,e in itertools.product((2,3,4),repeat=2):
            if (d,e) not in ALLOWED:
                clauses.append([-selections[0][d],-selections[1][e]])
    # Sound: label B and C independently so the first A-event's neighbours
    # come first in each. This does not constrain any other row.
    for part in parts[1:]:
        for u,v in zip(part,part[1:]):
            clauses.append([edge(0,u),-edge(0,v)])
    pairs=[(u,v) for u,v in itertools.combinations(range(30),2) if u//10!=v//10]
    solver=Cadical153(bootstrap_with=clauses)
    rng=np.random.default_rng(seed)
    solver.set_phases([edge(u,v)*(1 if rng.random()<0.34 else -1) for u,v in pairs])
    return solver, edge, pairs, parts


def control():
    from verify_all import tricayley_z2z6
    S=(((0,0),(0,3),(0,4),(1,4)),((0,0),(0,3),(0,4),(1,4)),
       ((0,1),(1,0),(1,1),(1,3)))
    V,E=tricayley_z2z6(S)
    index={v:i for i,v in enumerate(V)}
    H=nx.Graph(); H.add_nodes_from(range(36)); H.add_edges_from((index[u],index[v]) for u,v in E)
    value=UniformMargin([list(range(12*k,12*k+12)) for k in range(3)]).solve(H)
    if abs(value)>1e-6:
        raise RuntimeError(f"H36 positive control failed: {value}")
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds",type=float,default=165)
    parser.add_argument("--seed",type=int,default=1)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--verify",action="store_true",help="exactly exclude the best saved graph")
    args=parser.parse_args()
    if args.verify:
        from dualcert import exact_upper
        best=json.loads(args.output.read_text())["best"]
        if best is None:
            raise RuntimeError("no saved graph to verify")
        H=nx.Graph(); H.add_nodes_from(range(30)); H.add_edges_from(best["edges"])
        assert not any(nx.triangles(H).values())
        assert max(map(len,nx.find_cliques(nx.complement(H))))==10
        for u in H:
            ds=tuple(sum(v//10==k for v in H[u]) for k in range(3) if k!=u//10)
            assert ds in ALLOWED
        for a,b in itertools.combinations(best["parts"],2):
            assert nx.is_connected(H.subgraph(a+b))
        bound,numeric=exact_upper(H)
        if bound is None or bound>=3:
            raise RuntimeError("no strict exact theta bound certified")
        print("best saved graph: exact theta upper",bound,"numeric dual",numeric)
        return
    result={"scope":"bounded heuristic SAT/SDP search; not a nonexistence proof",
            "ansatz":"uniform probabilities; connected pairs; triangle-free; alpha<=10",
            "seed":args.seed,"seconds_budget":args.seconds,"h36_control_margin":control(),
            "models":0,"scored":0,"independence_cuts":0,"connectivity_cuts":0,"best":None}
    solver,edge,pairs,parts=build(args.seed)
    margin=UniformMargin(parts)
    start=time.monotonic()
    def save():
        result["elapsed_seconds"]=time.monotonic()-start
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    save()
    while time.monotonic()-start<args.seconds:
        solver.conf_budget(5000)
        sat=solver.solve_limited()
        if sat is None:
            save(); continue
        if not sat:
            result["termination"]="SAT encoding exhausted; no independent proof log"
            break
        result["models"]+=1
        mdl={x for x in solver.get_model() if x>0}
        H=nx.Graph(); H.add_nodes_from(range(30)); H.add_edges_from((u,v) for u,v in pairs if edge(u,v) in mdl)
        bad=False
        for a,b in itertools.combinations(parts,2):
            K=H.subgraph(a+b)
            components=list(nx.connected_components(K))
            if len(components)>1:
                X=components[0]
                solver.add_clause([edge(u,v) for u in X for v in set(K)-X if u//10!=v//10])
                result["connectivity_cuts"]+=1; bad=True
        count=0
        for independent in nx.find_cliques(nx.complement(H)):
            if len(independent)>10:
                X=independent[:11]
                solver.add_clause([edge(u,v) for u,v in itertools.combinations(X,2) if u//10!=v//10])
                count+=1; bad=True
                if count>=50: break
        result["independence_cuts"]+=count
        if bad:
            save(); continue
        value=margin.solve(H)
        result["scored"]+=1
        if result["best"] is None or value>result["best"]["uniform_margin_numeric"]:
            result["best"]={"uniform_margin_numeric":value,"edges":sorted(map(list,H.edges())),"parts":parts}
            print("new best",value,"scored",result["scored"],"models",result["models"],flush=True)
        # Block only this complete graph assignment; negative numerical margins
        # are not used to learn unsound stronger clauses.
        solver.add_clause([-edge(u,v) if edge(u,v) in mdl else edge(u,v) for u,v in pairs])
        save()
        if value>=-1e-7:
            result["termination"]="numerical candidate; needs exact certification"
            break
    else:
        result["termination"]="time budget"
    save(); solver.delete()
    print(json.dumps({k:v for k,v in result.items() if k!="best"}),flush=True)


if __name__=="__main__":
    main()

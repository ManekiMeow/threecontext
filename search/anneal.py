# Local search over edge-maximal tripartite triangle-free graphs, maximising theta(Hbar).
import sys, random, pickle, math, networkx as nx
from theta import theta_comp
def can_add(H,col,u,v):
    return col[u]!=col[v] and not H.has_edge(u,v) and not (set(H[u])&set(H[v]))
def maximize(H,col,rng):
    nodes=list(H); pairs=[(u,v) for i,u in enumerate(nodes) for v in nodes[i+1:] if col[u]!=col[v]]
    rng.shuffle(pairs)
    for u,v in pairs:
        if can_add(H,col,u,v): H.add_edge(u,v)
    return H
def run(H,col,seed,iters,T0=0.01,kick=3,tag=""):
    rng=random.Random(seed)
    H=maximize(H.copy(),col,rng); t=theta_comp(H)[0]; best=(t,H.copy())
    for it in range(iters):
        H2=H.copy(); E=list(H2.edges())
        for e in rng.sample(E,min(kick,len(E))): H2.remove_edge(*e)
        maximize(H2,col,rng); t2=theta_comp(H2)[0]
        T=T0*(1-it/iters)+1e-6
        if t2>=t or rng.random()<math.exp((t2-t)/T):
            H,t=H2,t2
            if t>best[0]:
                best=(t,H.copy()); print(tag,"it",it,"best",t,flush=True)
                pickle.dump((col,list(H.edges())),open(f"ann_best_{tag}.pkl","wb"))
                if t>3-1e-6: print(tag,"HIT n=",len(H),flush=True); return best
    return best
if __name__=="__main__":
    mode=sys.argv[1]; seed=int(sys.argv[2]); iters=int(sys.argv[3])
    if mode.startswith("from:"):
        import json
        d=json.load(open(mode[5:])); H=nx.Graph()
        for a,b,_ in d["edges"]: H.add_edge(tuple(a),tuple(b))
        rng=random.Random(seed); v=rng.choice(sorted(H)); H.remove_node(v)
        col={u:u[0] for u in H}
        run(H,col,seed,iters,tag=f"from{len(H)}_{seed}")
    else:
        m=list(map(int,mode.split(",")))
        col={}; H=nx.Graph()
        for k,mk in enumerate(m):
            for i in range(mk): col[(k,i)]=k; H.add_node((k,i))
        run(H,col,seed,iters,tag=f"{mode}_{seed}")

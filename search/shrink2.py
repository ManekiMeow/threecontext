import sys, random, json, pickle, networkx as nx
from tsearch_graph import graph
from theta import theta_comp
q=int(sys.argv[1]); S=json.loads(sys.argv[2]); seed=int(sys.argv[3]); random.seed(seed)
H=graph(q,S); TOL=3-2e-6
while True:
    vs=list(H); random.shuffle(vs); done=False
    for v in vs:
        H2=H.copy(); H2.remove_node(v)
        if theta_comp(H2)[0]>TOL: H=H2; done=True; print("n=",len(H),flush=True); break
    if not done: break
tag=f"q{q}_{seed}_"+"_".join("".join(map(str,s)) for s in S)
print("FINAL",len(H),flush=True); pickle.dump((q,S,sorted(H)),open(f"shr_{tag}.pkl","wb"))

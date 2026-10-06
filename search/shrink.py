import sys, random, networkx as nx
from tsearch_graph import graph
from theta import theta_comp
seed=int(sys.argv[1]); random.seed(seed)
q=13; S=((0,1,3,9),(0,1,10),(1,6,8))
H=graph(q,S)
TOL=3-2e-6
while True:
    vs=list(H); random.shuffle(vs); done=False
    for v in vs:
        H2=H.copy(); H2.remove_node(v)
        if min(dict(H2.degree()).values(), default=0)<2: pass
        t,_=theta_comp(H2)
        if t>TOL:
            H=H2; done=True; print(seed,"n=",len(H),flush=True); break
    if not done: break
print("FINAL",seed,len(H),sorted(H),flush=True)
nx.write_graph6(nx.convert_node_labels_to_integers(H,ordering="sorted"),f"shrunk_{seed}.g6")
import pickle; pickle.dump(sorted(H),open(f"shrunk_{seed}.pkl","wb"))

import sys, itertools, networkx as nx
from tcount import gen
from theta import theta_comp
def graph(q,S):
    H=nx.Graph()
    for a in range(3):
        for x in range(q):
            for s in S[a]: H.add_edge((a,x),((a+1)%3,(x+s)%q))
    return H
q=int(sys.argv[1]); sizes=tuple(map(int,sys.argv[2:5]))
best=(0,None); cnt=0
for S in gen(q,sizes):
    H=graph(q,S)
    if len(H)<3*q: continue
    t,_=theta_comp(H); cnt+=1
    if t>best[0]: best=(t,S)
    if t>3-1e-5: print("HIT",q,S,t,flush=True)
print(q,sizes,"tested",cnt,"best",best,flush=True)

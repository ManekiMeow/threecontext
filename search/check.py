import sys, networkx as nx
from theta import theta_comp
best=0
for line in sys.stdin:
    H=nx.from_graph6_bytes(line.strip().encode())
    t,_=theta_comp(H)
    best=max(best,t)
    if t>3-1e-5: print("HIT",line.strip(),t,flush=True)
print("max theta(Hbar)=",best)

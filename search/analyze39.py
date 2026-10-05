import numpy as np, networkx as nx
from tsearch_graph import graph
from theta import theta_comp
q=13; S=((0,1,3,9),(0,1,10),(1,6,8))
H=graph(q,S); nodes=sorted(H)
idx={v:i for i,v in enumerate(nodes)}
W=np.zeros((39,39))
wts=[1/4,1/3,1/3]
for a in range(3):
    for x in range(q):
        for s in S[a]:
            i,j=idx[(a,x)],idx[((a+1)%3,(x+s)%q)]; W[i,j]=W[j,i]=wts[a]
ev=np.linalg.eigvalsh(W); print("spectrum of W:",np.round(ev,4))
print("rank of W+I (=dimension):",(ev>-1+1e-9).sum())
for v in [(0,0),(1,0),(2,0)]:
    H2=H.copy(); H2.remove_node(v); print("delete",v,"theta=",theta_comp(H2)[0])
# edge orbit deletion
for a in range(3):
    for s in S[a]:
        H2=H.copy()
        for x in range(q): H2.remove_edge((a,x),((a+1)%3,(x+s)%q))
        print("delete edge orbit",a,s,"theta=",round(theta_comp(H2)[0],6))

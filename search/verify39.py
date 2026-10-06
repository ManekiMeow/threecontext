import networkx as nx, numpy as np
from tsearch_graph import graph
from theta import theta_comp
for S in [((0,1,3,9),(0,1,10),(1,6,8)),((0,1,3,9),(0,1,4),(1,2,5,7)),((0,1,2,4),(0,1,3,9),(1,5,7)),((0,1,3,9),(0,1,3,9),(2,5,6))]:
    H=graph(13,S)
    tri=sum(nx.triangles(H).values())//3
    t,B=theta_comp(H)
    ev=np.linalg.eigvalsh(B); r=(ev>1e-6*ev.max()).sum()
    print(S,"n",len(H),"edges",H.number_of_edges(),"triangles",tri,"degrees",sorted(set(dict(H.degree()).values())),"theta",t,"rank",r, "iso-class", nx.weisfeiler_lehman_graph_hash(H))

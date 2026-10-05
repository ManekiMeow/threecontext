import sys, random, pickle, json, networkx as nx
from groupsearch import groups, tricayley
from theta import theta_comp
from gcert import rational_certificate
mul,e=groups()["Z2xZ6"]
S=((0,1,3,7),(0,1,3,7),(4,6,9,10))
H=tricayley(mul,S)
print("n",len(H),"edges",H.number_of_edges(),"triangles",sum(nx.triangles(H).values())//3,"degrees",sorted(set(dict(H.degree()).values())),"theta",theta_comp(H)[0],flush=True)
for v in [(0,0),(1,0),(2,0)]:
    H2=H.copy(); H2.remove_node(v); print("delete",v,theta_comp(H2)[0],flush=True)
cert,tau=rational_certificate(H,{v:v[0] for v in H})
print("margin",tau,"cert",cert is not None)
if cert:
    nodes,E,w,c=cert
    json.dump({"group":"Z2xZ6","mul":mul,"identity":e,"S":S,"nodes":[list(v) for v in nodes],"edges":[[list(nodes[a]),list(nodes[b]),x] for (a,b),x in zip(E,w)],"c":c},open("z2z6_36_cert.json","w"))
    print(sorted(set(w)), sorted(set(c)))

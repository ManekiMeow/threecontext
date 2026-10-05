import pickle, json, sys
from tsearch_graph import graph
from gcert import rational_certificate
q,S,keep=pickle.load(open(sys.argv[1],"rb"))
H=graph(q,[tuple(s) for s in S]).subgraph([tuple(v) for v in keep]).copy()
colour={v:v[0] for v in H}
cert,tau=rational_certificate(H,colour)
print("n",len(H),"edges",H.number_of_edges(),"margin",tau,"certificate:",cert is not None, "removed:",sorted(set(graph(q,[tuple(s) for s in S]))-set(H)))
if cert:
    nodes,E,w,cw=cert
    json.dump({"q":q,"S":S,"nodes":[list(v) for v in nodes],"edges":[[list(nodes[a]),list(nodes[b]),x] for (a,b),x in zip(E,w)],"c":cw},open(sys.argv[1].replace(".pkl","_cert.json"),"w"))

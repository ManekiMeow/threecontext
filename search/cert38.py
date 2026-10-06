import pickle, json, sys
from tsearch_graph import graph
from gcert import rational_certificate
q=13; S=((0,1,3,9),(0,1,10),(1,6,8))
H=graph(q,S)
keep=pickle.load(open(sys.argv[1],"rb"))
H=H.subgraph(keep).copy()
colour={v:v[0] for v in H}
cert,tau=rational_certificate(H,colour)
print("n",len(H),"margin",tau,"certificate found:",cert is not None)
if cert:
    nodes,E,w,cw=cert
    json.dump({"nodes":[list(v) for v in nodes],"edges":[[list(nodes[a]),list(nodes[b]),x] for (a,b),x in zip(E,w)],"c":cw},open(sys.argv[1].replace(".pkl","_cert.json"),"w"))
    print(sorted(set(w)))

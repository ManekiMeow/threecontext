import json, networkx as nx
from dualcert import exact_upper
d=json.load(open("shr_q13_1_0139_014_1257_cert.json")); H=nx.Graph()
for a,b,_ in d["edges"]: H.add_edge(tuple(a),tuple(b))
print("n",len(H),"edges",H.number_of_edges(), "triangles", sum(nx.triangles(H).values()))
worst=0
for v in sorted(H):
    H2=H.copy(); H2.remove_node(v); tq,t=exact_upper(H2)
    worst=max(worst,float(tq) if tq else 9)
    print(v,"theta(Hbar - v) <=",tq,"(numeric",round(t,6),")",flush=True)
print("max certified upper bound after deleting one vertex:",worst)

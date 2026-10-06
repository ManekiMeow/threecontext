import networkx as nx
def graph(q,S):
    H=nx.Graph()
    for a in range(3):
        for x in range(q):
            for s in S[a]: H.add_edge((a,x),((a+1)%3,(x+s)%q))
    return H

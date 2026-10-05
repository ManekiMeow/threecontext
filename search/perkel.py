import itertools, networkx as nx, numpy as np
def build(S, m):
    G = nx.Graph()
    for a in range(3):
        for x in range(19):
            for s in S:
                d = (s * pow(m, a, 19)) % 19
                G.add_edge((a, x), ((a+1) % 3, (x+d) % 19))
    return G
def perkel():
    for S in itertools.combinations(range(1,19),3):
        for m in range(1,19):
            G = build(S, m)
            if G.number_of_edges()!=171 or nx.girth(G)!=5: continue
            if nx.girth(G)!=5 or nx.diameter(G)!=3: continue
            # check intersection array
            D = dict(nx.all_pairs_shortest_path_length(G))
            v=(0,0); ok=True
            arr=set()
            for u in G:
                d=D[v][u]
                b=sum(1 for w in G[u] if D[v][w]==d+1); c=sum(1 for w in G[u] if D[v][w]==d-1)
                arr.add((d,b,c))
            if arr=={(0,6,0),(1,5,1),(2,2,1),(3,0,3)}:
                return G
if __name__=="__main__":
    G=perkel(); print(sorted(set(np.round(np.linalg.eigvalsh(nx.to_numpy_array(G)),4))))
    nx.write_graph6(G,"perkel.g6")

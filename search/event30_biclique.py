"""A dimension-independent exact obstruction for thirty-event GHZ supports.

No bipartite graph between two contexts may have a complete bipartite connected
component K_(r,s), with r,s>0. In a theta-three realization write f_i=sqrt(c_i)g_i
so the vectors in each context sum to psi. For a pair-support component X,Y,
    v=sum_(x in X)f_x=sum_(y in Y)f_y.
Indeed the squared norms of the two sums and their inner product all coincide,
by completing psi in the other context. Their common norm is positive when all
events have positive probability, as is forced at30 events by the event bound.
If X times Y is complete, triangle-freeness means every event z in the third
context misses X or misses Y. The corresponding expansion of v gives <g_z,v>=0.
Thus v is orthogonal to psi, contrary to <psi,v>=sum_(x in X)c_x>0.

The condition is stronger than balanced-component recolouring: even K_(1,1)
is excluded. It remains sound when allowed edges acquire zero overlaps.

Standard-Python command: python search/event30_biclique.py graph.json
Accepts an edges/parts record or an edges record with the default contiguous
three contexts of size10. Prints a directly checkable witness or null.
"""
import argparse
import json
from itertools import combinations
from pathlib import Path


def graph_data(edges, parts):
    assert len(parts)==3 and all(len(part)==10 for part in parts)
    colour={v:k for k,part in enumerate(parts) for v in part}
    assert len(colour)==30
    adj={v:set() for v in colour}
    for u,v in edges:
        assert u in adj and v in adj and colour[u]!=colour[v]
        adj[u].add(v);adj[v].add(u)
    assert all(not (adj[u] & adj[v]) for u in adj for v in adj[u])
    return adj,colour


def verify(edges, parts, witness):
    adj,colour=graph_data(edges,parts)
    a,b=witness["context_pair"]
    assert 0<=a<b<3
    X,Y=set(witness["left"]),set(witness["right"])
    assert X and Y and X<=set(parts[a]) and Y<=set(parts[b])
    assert all(adj[x] & set(parts[b])==Y for x in X)
    assert all(adj[y] & set(parts[a])==X for y in Y)
    return True


def find(edges, parts):
    adj,colour=graph_data(edges,parts)
    for a,b in combinations(range(3),2):
        pair=set(parts[a]) | set(parts[b])
        unseen=set(pair)
        while unseen:
            start=min(unseen)
            component={start};pending=[start]
            while pending:
                u=pending.pop()
                for v in (adj[u] & pair)-component:
                    component.add(v);pending.append(v)
            unseen-=component
            X=sorted(component & set(parts[a]));Y=sorted(component & set(parts[b]))
            if X and Y and all(y in adj[x] for x in X for y in Y):
                witness={"kind":"complete_bipartite_pair_component",
                         "context_pair":[a,b],"left":X,"right":Y}
                assert verify(edges,parts,witness)
                return witness
    return None


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("graph",type=Path)
    args=parser.parse_args()
    data=json.loads(args.graph.read_text())
    parts=data.get("parts",[list(range(10*k,10*k+10)) for k in range(3)])
    print(json.dumps(find(data["edges"],parts)))

"""Exact exclusions of large-orbit symmetries in thirty-event GHZ graphs.

Standard Python: python search/event30_large_cycle.py

Assume three contexts of10 positive-probability events and a context-preserving
automorphism GROUP with one transitive7-orbit and three individually fixed
vertices in each context. The group need not be cyclic.
Probabilities are positive and strictly below1/9 by Corollary32. No Gram
averaging is needed: a fixed vertex sees either the whole moving orbit or
none of it, regardless of the probabilities within that orbit.

For a fixed event touching a moving orbit, its neighbourhood weight alpha is
greater than2/3. Lemma14 gives
    c <= h(alpha) = alpha*(1-alpha)/(1+3*alpha) < h(2/3)=2/27.
Here the other two sector weights were optimized by making them equal, and h
decreases for alpha>1/3. Conversely, fixed-event probabilities in each context
sum to more than2/9 (seven moving events each have probability<1/9). Thus every
context has a HEAVY fixed event of probability>2/27, which has no moving
neighbours. Every heavy event needs at least two neighbours in each other fixed
triple: a single neighbour has weight<1/9 and would give c<h(1/9)=2/27.

Two heavy events in different contexts cannot be adjacent, since their two
neighbour sets in the third fixed triple intersect. There is exactly one heavy
event per context: two would leave a heavy event elsewhere at most one allowed
neighbour. Each heavy event therefore sees at most the two nonheavy events of
each other context, giving alpha,beta<=4/27. The harmonic bound in Lemma14 is
increasing in each variable on this rectangle (rho>=19/27>alpha,beta). Hence
    c <= [27/4+27/4+27/19]**(-1) =38/567 <2/27,
a contradiction. Faithfulness on each context implies that any AUTOMORPHISM
of order7 generates a group with the assumed action. Cauchy's theorem excludes7 as a divisor
of the order of the context-preserving automorphism group.

The same method excludes any group action (8,1,1) in ALL THREE contexts. Fixed
probabilities sum>1/9, so every context has a fixed heavy event c>1/18. Touching
the moving orbit would give c<h(7/9)=7/135<1/18. A heavy event must touch another
heavy event in each other fixed pair: one nonheavy neighbour has weight<=1/18,
giving c<=h(1/18)=17/378<1/18. A triangle-free tripartite heavy graph with at most
two vertices per part and no missing colour at any vertex needs two vertices
in every part. Thus all six fixed events are heavy and have no moving
neighbours. The30-event support is disconnected, contradicting the event lower
bound. Mixed order-eight actions involving a2-cycle are not excluded here.

Finally, suppose a group has one transitive9-orbit and one individually fixed
vertex in each context. A fixed vertex's neighbourhood in another context
is nonempty, proper and invariant, so consists of either the sole fixed vertex
or the entire moving orbit. Each moving-to-moving block is nonempty, since
otherwise its nine rows coincide. Thus a fixed vertex cannot touch both other
moving orbits, so it touches at least one other fixed vertex. The graph on the
three fixed vertices is triangle-free with minimum degree1, hence a path P3.
Let its centre be a and endpoints b,c. The neighbourhoods of b have weights
c_a and1-c_c, so its residual weight is c_c-c_a>0 by Lemma14 and c_b>0.
But c has only a as neighbour in A, and Lemma14 gives c_c<c_a. Contradiction.
In particular this excludes an ELEMENT of order9 by faithfulness, and excludes
a transitive C3xC3 action on9 vertices plus one fixed vertex in each context.
It does not exclude divisibility of the automorphism group order by9 in general.
"""
from fractions import Fraction as Q
from itertools import combinations, product


def h(alpha):
    return alpha*(1-alpha)/(1+3*alpha)


def arithmetic():
    threshold7=Q(2,27)
    assert (1-7*Q(1,9))/3 == threshold7
    assert h(Q(2,3)) == threshold7
    assert h(Q(1,9)) == threshold7
    cap=Q(4,27)
    rho=1-2*cap
    assert rho==Q(19,27) and rho>cap
    double_bound=1/(1/cap+1/cap+1/rho)
    assert double_bound==Q(38,567)<threshold7
    print("C7: heavy threshold2/27; final bound38/567<2/27: impossible")

    threshold8=Q(1,18)
    assert (1-8*Q(1,9))/2 == threshold8
    assert h(Q(7,9))==Q(7,135)<threshold8
    assert h(threshold8)==Q(17,378)<threshold8
    print("C8(8,1,1)^3: heavy threshold1/18; all fixed vertices must be disconnected from movement")


def heavy_combinatorics():
    subsets=[set(S) for size in (2,3) for S in combinations(range(3),size)]
    assert all(S & T for S,T in product(subsets,repeat=2))

    counts={}
    for sizes in product((1,2),repeat=3):
        parts=[];offset=0
        for size in sizes:
            parts.append(set(range(offset,offset+size)));offset+=size
        colour={u:k for k,part in enumerate(parts) for u in part}
        edges=[(u,v) for u,v in combinations(range(offset),2) if colour[u]!=colour[v]]
        valid=0
        for mask in range(1 << len(edges)):
            adj=[set() for _ in range(offset)]
            for k,(u,v) in enumerate(edges):
                if mask >> k & 1:
                    adj[u].add(v);adj[v].add(u)
            if any(adj[u] & adj[v] for u in range(offset) for v in adj[u]):
                continue
            if any(not (adj[u] & part) for u in range(offset)
                   for k,part in enumerate(parts) if k!=colour[u]):
                continue
            valid+=1
        if valid:
            counts[sizes]=valid
    assert counts=={(2,2,2):4}
    print("exact heavy-graph enumeration: only(2,2,2), four labelled C6 patterns")


def order_nine():
    possible=[]
    edges=list(combinations(range(3),2))
    for mask in range(8):
        adj=[set() for _ in range(3)]
        for k,(u,v) in enumerate(edges):
            if mask >> k & 1:
                adj[u].add(v);adj[v].add(u)
        if all(adj) and not all(len(N)==2 for N in adj):
            assert sorted(map(len,adj))==[1,1,2]
            centre=next(u for u in range(3) if len(adj[u])==2)
            b,c=sorted(adj[centre])
            assert adj[b]=={centre} and adj[c]=={centre}
            possible.append(centre)
    assert sorted(possible)==[0,1,2]
    # In coordinates(c_a,c_b,c_c), both of these linear forms would have
    # to be positive: rho_b=c_c-c_a and alpha_c-c_c=c_a-c_c.
    rho_b=(-1,0,1)
    single_neighbour_gap=(1,0,-1)
    assert tuple(a+b for a,b in zip(rho_b,single_neighbour_gap))==(0,0,0)
    print("C9: three possible fixed P3 patterns; rho_b=c_c-c_a>0 and c_a-c_c>0 contradict")


if __name__=="__main__":
    arithmetic()
    heavy_combinatorics()
    order_nine()

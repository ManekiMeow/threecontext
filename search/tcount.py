import itertools, sys, math
def canon(q,S):
    # canonical form under: shifts (t0,t1,t2), sum t=0 ; unit multipliers; cyclic rotation of classes; reversal
    best=None
    units=[u for u in range(1,q) if math.gcd(u,q)==1]
    variants=[]
    S=[tuple(s) for s in S]
    rots=[S, [S[1],S[2],S[0]], [S[2],S[0],S[1]]]
    rev=lambda T: [tuple(sorted((-x)%q for x in T[2])),tuple(sorted((-x)%q for x in T[1])),tuple(sorted((-x)%q for x in T[0]))]
    allr=rots+[rev(r) for r in rots]
    for T in allr:
        for u in units:
            U=[sorted((u*x)%q for x in Ta) for Ta in T]
            for t0 in U[0]:
                for t1 in U[1]:
                    A=tuple(sorted((x-t0)%q for x in U[0])); B=tuple(sorted((x-t1)%q for x in U[1]))
                    C=tuple(sorted((x+t0+t1)%q for x in U[2]))
                    key=(A,B,C)
                    if best is None or key<best: best=key
    return best
def gen(q,sizes):
    seen=set()
    a,b,c=sizes
    for S0r in itertools.combinations(range(1,q),a-1):
        S0=(0,)+S0r
        for S1r in itertools.combinations(range(1,q),b-1):
            S1=(0,)+S1r
            forb={(-(x+y))%q for x in S0 for y in S1}
            allowed=[z for z in range(q) if z not in forb]
            for S2 in itertools.combinations(allowed,c):
                k=canon(q,(S0,S1,S2))
                if k not in seen:
                    seen.add(k); yield k
if __name__=="__main__":
    q=int(sys.argv[1]); sizes=tuple(map(int,sys.argv[2:5]))
    n=0
    for k in gen(q,sizes): n+=1
    print(q,sizes,n)

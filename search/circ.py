import numpy as np, itertools, sys
from scipy.optimize import linprog
def theta_circ(n,S):
    # variables b(x) for x in orbits {x,-x}, x=1..n//2 not in S ; b(0)=1
    xs=[x for x in range(1,n//2+1) if x not in S]
    mult=[1 if 2*x==n else 2 for x in xs]
    c=-np.array(mult,float)
    A=[];bb=[]
    for k in range(n):
        row=[-m*np.cos(2*np.pi*k*x/n) for x,m in zip(xs,mult)]
        A.append(row); bb.append(1.0)
    r=linprog(c,A_ub=A,b_ub=bb,bounds=[(None,None)]*len(xs),method="highs")
    return 1-r.fun
def trifree(n,S):
    Sf=set(S)|{(-s)%n for s in S}
    return not any(((a+b)%n) in Sf for a in Sf for b in Sf)
def threecol(n,S):
    import networkx as nx
    Sf=set(S)|{(-s)%n for s in S}
    # backtracking 3-coloring
    adj=[[ (i+s)%n for s in Sf] for i in range(n)]
    col=[-1]*n
    sys.setrecursionlimit(10000)
    def bt(i):
        if i==n: return True
        for c in range(3):
            if all(col[j]!=c for j in adj[i]):
                col[i]=c
                if bt(i+1): return True
                col[i]=-1
        return False
    return bt(0)
maxdeg=int(sys.argv[2]) if len(sys.argv)>2 else 8
for n in range(6,int(sys.argv[1])+1,3):
    for k in range(1,maxdeg//2+1):
        for S in itertools.combinations(range(1,(n+1)//2),k):
            if not trifree(n,S): continue
            t=theta_circ(n,S)
            if abs(t-n/3)<1e-6:
                if threecol(n,S): print(n,S,"theta(H)=",t,flush=True)

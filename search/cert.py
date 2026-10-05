# Find rational Hoffman-type certificate W for theta(Hbar)>=3 on a tri-Cayley graph:
# W supported on edges, W 1 = 2*1, W + I PSD  =>  theta(Hbar) >= 1 - lmax/lmin = 3 (Lovasz 1979, Thm 6).
import numpy as np, cvxpy as cp, sys, json
from fractions import Fraction as Fr
def margin_sdp(q,S):
    b=[cp.Variable(len(Sa)) for Sa in S]; tau=cp.Variable()
    cons=[cp.sum(x)==1 for x in b]
    for k in range(1,q//2+1):
        f=[np.exp(2j*np.pi*k*np.array(Sa)/q)@b[a] for a,Sa in enumerate(S)]
        M=cp.Variable((3,3),hermitian=True)
        cons+=[M>>0, M[0,0]==1-tau,M[1,1]==1-tau,M[2,2]==1-tau, M[0,1]==f[0],M[1,2]==f[1],M[2,0]==f[2]]
    p=cp.Problem(cp.Maximize(tau),cons); p.solve(solver="CLARABEL")
    return tau.value,[x.value for x in b]
def exact_psd(A):
    # exact LDL^T / Gaussian elimination on Fractions; returns True iff PSD
    n=len(A); A=[row[:] for row in A]
    for i in range(n):
        if A[i][i]<0: return False
        if A[i][i]==0:
            if any(A[i][j]!=0 for j in range(i+1,n)): return False
            continue
        for j in range(i+1,n):
            if A[j][i]!=0:
                r=A[j][i]/A[i][i]
                for k in range(i,n): A[j][k]-=r*A[i][k]
    return True
def build_W(q,S,w):
    n=3*q; idx=lambda a,x:a*q+x
    W=[[Fr(0)]*n for _ in range(n)]
    for a in range(3):
        for x in range(q):
            for s,ws in zip(S[a],w[a]):
                i,j=idx(a,x),idx((a+1)%3,(x+s)%q); W[i][j]=W[j][i]=ws
    return W
if __name__=="__main__":
    q=int(sys.argv[1]); S=json.loads(sys.argv[2])
    tau,b=margin_sdp(q,S); print("margin tau=",tau,"weights",b)
    for den in [6,12,24,30,60,120,360,1000,10000]:
        w=[]
        for x in b:
            r=[Fr(round(v*den),den) for v in x[:-1]]; r.append(1-sum(r)); w.append(r)
        W=build_W(q,S,w); n=len(W)
        A=[[W[i][j]+(1 if i==j else 0) for j in range(n)] for i in range(n)]
        ok=exact_psd(A)
        rows=all(sum(W[i])==2 for i in range(n))
        print("den",den,"weights",[[str(v) for v in r] for r in w],"W+I PSD:",ok,"W1=2*1:",rows)
        if ok and rows: break

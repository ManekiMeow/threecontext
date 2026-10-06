# Tri-Cayley graphs over Z_q: vertices (a,x), a in Z3, x in Z_q; (a,x)~(a+1,x+s), s in S_a.
# |S_a|=3 forces uniform weights 1/3; theta(Hbar)=3 iff for all k!=0 the 3x3 Fourier block has lambda_min>=-1.
import numpy as np, itertools, sys
def ok(q,S):
    w=np.exp(2j*np.pi*np.arange(q)/q)
    for k in range(1,q):
        f=[sum(w[(k*s)%q] for s in Sa)/3 for Sa in S]
        M=np.array([[0,f[0],np.conj(f[2])],[np.conj(f[0]),0,f[1]],[f[2],np.conj(f[1]),0]])
        if np.linalg.eigvalsh(M)[0] < -1-1e-9: return False
    return True
for q in range(int(sys.argv[1]),int(sys.argv[2])+1):
    hits=set()
    subs=list(itertools.combinations(range(q),3))
    S0s=[s for s in subs if s[0]==0]
    for S0 in S0s:
        for S1 in S0s:
            sums01={(a+b)%q for a in S0 for b in S1}
            for S2 in subs:
                if any((-c)%q in sums01 for c in S2): continue
                if ok(q,(S0,S1,S2)): hits.add((S0,S1,S2)); print(q,S0,S1,S2,flush=True)
    print("q",q,"hits",len(hits),flush=True)

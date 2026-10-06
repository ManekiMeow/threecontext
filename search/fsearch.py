import sys, itertools
from multiprocessing import Pool
from tcount import gen
from ftheta import ftheta
q=int(sys.argv[1]); sizes=tuple(map(int,sys.argv[2:5]))
def work(S):
    try: return ftheta(q,S),S
    except Exception as e: return -1,S
if __name__=="__main__":
    cands=list(gen(q,sizes))
    with Pool(4) as P: res=P.map(work,cands,chunksize=4)
    res.sort(reverse=True)
    for t,S in res:
        if t>3-1e-6: print("HIT",q,S,t,flush=True)
    print(q,sizes,"tested",len(res),"best",res[0] if res else None,flush=True)

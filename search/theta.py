import numpy as np, cvxpy as cp, networkx as nx
def theta_comp(H, solver="CLARABEL"):
    """Lovasz theta of the complement of H: max sum(B) s.t. B psd, tr B=1, B_ij=0 for non-edges of H."""
    nodes=list(H); idx={v:i for i,v in enumerate(nodes)}; n=len(nodes)
    B=cp.Variable((n,n),symmetric=True)
    mask=np.ones((n,n))
    for u,v in H.edges(): mask[idx[u],idx[v]]=mask[idx[v],idx[u]]=0
    np.fill_diagonal(mask,0)
    cons=[B>>0, cp.trace(B)==1, cp.multiply(mask,B)==0]
    p=cp.Problem(cp.Maximize(cp.sum(B)),cons); p.solve(solver=solver)
    return p.value, B.value
if __name__=="__main__":
    H=nx.read_graph6("perkel.g6"); print("perkel theta(Hbar)=",theta_comp(H)[0])

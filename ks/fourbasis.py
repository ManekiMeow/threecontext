# Shared helpers for unions of four orthonormal bases (four-context KS problem).
# A configuration is a list of unit vectors `vecs` plus `cover`, four lists of indices (the covering bases).
# Contexts of the union = ALL orthonormal bases contained in it (covering bases and internal bases).
import itertools, numpy as np, networkx as nx
from pysat.solvers import Cadical153

TOL = 1e-9

def ortho_graph(vecs, tol=TOL):
    A = np.array(vecs); M = np.abs(A.conj() @ A.T) < tol
    G = nx.Graph(); G.add_nodes_from(range(len(vecs)))
    G.add_edges_from((i, j) for i in range(len(vecs)) for j in range(i + 1, len(vecs)) if M[i, j])
    return G

def all_bases(G, d):
    # an orthonormal basis is a d-clique, which is automatically a maximal clique
    return [frozenset(c) for c in nx.find_cliques(G) if len(c) == d]

def colouring(G, bases):
    """A KS colouring (set of rays coloured 1) or None: exactly one 1 per basis, no orthogonal pair of 1s."""
    s = Cadical153()
    for b in bases:
        lits = [v + 1 for v in b]; s.add_clause(lits)
        for a, c in itertools.combinations(lits, 2): s.add_clause([-a, -c])
    for u, v in G.edges(): s.add_clause([-(u + 1), -(v + 1)])
    ok = s.solve(); m = s.get_model() if ok else None; s.delete()
    return None if m is None else {l - 1 for l in m if l > 0}

def k4s(G, cover):
    """Pairwise non-orthogonal transversals (one ray from each covering basis)."""
    out = []
    for t in itertools.product(*cover):
        if all(not G.has_edge(a, b) for a, b in itertools.combinations(t, 2)): out.append(t)
    return out

def pieces(basis, cover):
    return sum(1 for c in cover if basis & set(c))

def born(vecs, rho):
    A = np.array(vecs); return np.real(np.einsum("ri,ij,rj->r", A.conj(), rho, A))

def mixture_with_marginals(T, n, marg):
    """LP: is there a probability distribution on the transversals T whose one-ray marginals equal marg?
    By the averaging lemma, feasibility with marg = Born probabilities of a state certifies 'not KS'."""
    from scipy.optimize import linprog
    if not T: return False
    A = np.zeros((n, len(T)))
    for j, t in enumerate(T):
        for r in t: A[r, j] = 1
    res = linprog(np.zeros(len(T)), A_eq=A, b_eq=marg, bounds=(0, None), method="highs")
    return res.status == 0

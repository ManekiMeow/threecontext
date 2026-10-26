# Numerical realizability of a prescribed orthogonality pattern for four bases of C^d.
# Basis 3 is the standard basis; bases 0-2 are polar factors of free complex matrices.  Only the prescribed
# zeros are fitted: extra orthogonalities can only add constraints, so a realization of a KS pattern is KS
# (RESULTS.md Lemma 8).  Usage: python3 realize.py d pattern.pkl [restarts]
import sys, pickle, numpy as np
from scipy.optimize import least_squares

def polar(A):
    W, _, Vh = np.linalg.svd(A); return W @ Vh

def unpack(x, d):
    z = x[: 3 * d * d] + 1j * x[3 * d * d:]
    return [polar(z[k * d * d:(k + 1) * d * d].reshape(d, d)) for k in range(3)] + [np.eye(d)]

def vectors(x, d):
    U = unpack(x, d); return [U[r // d][:, r % d] for r in range(4 * d)]

def residuals(x, d, orth):
    v = vectors(x, d); out = []
    for r, s in orth:
        z = np.vdot(v[r], v[s]); out += [z.real, z.imag]
    return np.array(out)

def realize(d, orth, restarts=20, seed=0):
    rng = np.random.default_rng(seed); best = (np.inf, None)
    for _ in range(restarts):
        x0 = rng.normal(size=6 * d * d)
        sol = least_squares(residuals, x0, args=(d, orth), method="trf", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                            max_nfev=4000)
        c = np.sqrt(2 * sol.cost)
        if c < best[0]: best = (c, sol.x)
        if c < 1e-11: break
    return best

if __name__ == "__main__":
    d = int(sys.argv[1]); orth = pickle.load(open(sys.argv[2], "rb"))
    res, x = realize(d, orth, int(sys.argv[3]) if len(sys.argv) > 3 else 20)
    print("prescribed orthogonal pairs:", len(orth), " best residual: %.3e" % res)
    if res < 1e-9:
        from fourbasis import ortho_graph, all_bases, colouring
        v = vectors(x, d); G = ortho_graph(v, 1e-7); Bs = all_bases(G, d)
        print("realized: orthogonal pairs", G.number_of_edges(), "bases", len(Bs), "KS:", colouring(G, Bs) is None)

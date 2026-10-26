# The group Gamma = < AGL(2,3), W > acting on C^20 (RESULTS.md 3.11).  AGL(2,3) acts on the H33 realization by
# permuting points and lines; W is the duality psi <-> a0.  Both preserve B2 = B-lines + betas and B3 = C-lines +
# gammas as sets of rays.  We close the group projectively (matrices up to a phase) and list the orbit of psi.
import itertools, numpy as np
import h33model as M
from h33w import W

def line_unitary(pt_map, line_map):
    """Unitary from a permutation of the AG(2,3) geometry: pt_map on points, line_map (i,c) -> (j,c')."""
    src = [M.psi] + [M.point[x] for x in M.pts] + [M.Bline[l] for l in sorted(M.Bline)] + [M.Cline[l] for l in sorted(M.Cline)]
    tgt = [M.psi] + [M.point[pt_map(x)] for x in M.pts] + [M.Bline[line_map(l)] for l in sorted(M.Bline)] + \
          [M.Cline[line_map(l)] for l in sorted(M.Cline)]
    S = np.array(src).T; T = np.array(tgt).T; U = T @ np.linalg.pinv(S)
    assert np.allclose(U.conj().T @ U, np.eye(20), atol=1e-10) and np.allclose(U @ S, T, atol=1e-10)
    return U

def affine(g, t):
    g = np.array(g); gi = np.round(np.linalg.inv(g) * np.linalg.det(g)).astype(int) * int(round(np.linalg.det(g))) % 3
    def pm(x): return tuple((g @ np.array(x) + np.array(t)) % 3)
    def lm(l):
        i, c = l; k = np.array(M.K[i]); kn = tuple((gi.T @ k) % 3)       # normal of the image line (up to sign)
        for j in range(4):
            for s in (1, -1):
                if tuple((s * np.array(M.K[j])) % 3) == kn:
                    x = next(x for x in M.pts if M.dot(M.K[i], x) == c)     # a point of the line
                    return (j, M.dot(M.K[j], pm(x)))
    return line_unitary(pm, lm)

def key(U):
    """Canonical form of U up to a global phase (normalised at the first entry of maximal modulus)."""
    m = np.round(np.abs(U.ravel()), 6); z = U.ravel()[int(np.argmax(m))]
    V = U * abs(z) / z; return tuple(np.round(V.ravel(), 6).view(float) + 0.0)

gens = [affine([[1, 0], [0, 1]], (1, 0)), affine([[1, 0], [0, 1]], (0, 1)), affine([[0, 2], [1, 0]], (0, 0)),
        affine([[1, 1], [0, 1]], (0, 0)), affine([[2, 0], [0, 1]], (0, 0)), W]

def close(gens, limit=200000):
    I = np.eye(20, dtype=complex); elems = {key(I): I}; frontier = [np.eye(20, dtype=complex)]
    while frontier:
        nxt = []
        for U in frontier:
            for g in gens:
                V = g @ U; kV = key(V)
                if kV not in elems:
                    elems[kV] = V; nxt.append(V)
                    if len(elems) > limit: return elems, False
        frontier = nxt
    return elems, True

if __name__ == "__main__":
    agl, ok = close(gens[:5]); print("AGL(2,3) acting on C^20: projective order", len(agl), "closed:", ok)
    # Gamma itself is large (the orbit of psi has 10000 rays, see h33orbit.py); its order is not needed here.

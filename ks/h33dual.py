# The dual H33 copies (RESULTS.md 3.11).  With the block-pure completions B' = {beta_k}, C' = {gamma_k}, every point a
# of H33 is itself a GHZ state for (B2, B3, B4): its neighbourhood is an exact copy of H33 whose "lines" are the 4
# B-lines through a plus the 8 beta's (and the 4 C-lines through a plus the 8 gamma's), and whose "points" are psi and
# 8 vectors y_x(a) = 2 u_x + 2 w_x - a  (u_x, w_x = projections of a onto the lines through x), which must lie in B4.
import itertools, numpy as np
import h33model as M

def dual_points(a0):
    a = M.point[a0]; L = []
    for i in range(4):
        kp, km = M.K[i], M.neg(M.K[i]); c0 = M.dot(M.K[i], a0)
        # offsets relative to a0: line through a0 = 0, beta_{k_i} = 1, beta_{-k_i} = 2 (labelling found by search)
        L.append(({0: M.Bline[(i, c0)], 1: M.beta[kp], 2: M.beta[km]}, {0: M.Cline[(i, c0)], 1: M.gamma[km], 2: M.gamma[kp]}))
    Y = {}
    for x in M.pts:
        u = sum(np.vdot(L[i][0][M.dot(M.K[i], x)], a) * L[i][0][M.dot(M.K[i], x)] for i in range(4))
        w = sum(np.vdot(L[i][1][M.dot(M.K[i], x)], a) * L[i][1][M.dot(M.K[i], x)] for i in range(4))
        Y[x] = 2 * u + 2 * w - a
    return Y

if __name__ == "__main__":
    rays = []
    def add(v):
        for j, r in enumerate(rays):
            if abs(abs(np.vdot(r, v)) - 1) < 1e-9: return j
        rays.append(v); return len(rays) - 1
    owner = {}
    for a0 in M.pts:
        Y = dual_points(a0); assert abs(abs(np.vdot(Y[(0, 0)], M.psi)) - 1) < 1e-12
        for x, v in Y.items():
            if x != (0, 0): owner.setdefault(add(v), []).append(a0)
    Gm = np.abs(np.array([[np.vdot(r, s) for s in rays] for r in rays]))
    vals = sorted({round(g, 6) for g in Gm[np.triu_indices(len(rays), 1)]})
    print(f"distinct new 'points' over all 9 dual copies: {len(rays)}; each shared by {sorted({len(o) for o in owner.values()})} copies")
    print("pairwise overlaps among them:", vals)

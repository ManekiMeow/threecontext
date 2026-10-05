# Z_q-invariant Lovasz theta of the complement of a tri-Cayley graph.
# Vertices (a,x), a in Z3, x in Z_q; (a,x)~(a+1,x+s) for s in S[a].
# The theta SDP is invariant under translations, so averaging an optimal solution
# gives an invariant optimum: the 3n x 3n PSD constraint splits into q Hermitian 3x3 blocks.
import numpy as np, cvxpy as cp


def ftheta(q, S):
    d = cp.Variable(3)
    b = [cp.Variable(len(Sa)) for Sa in S]
    cons = [q * cp.sum(d) == 1]
    for k in range(q // 2 + 1):
        f = [np.exp(2j * np.pi * k * np.array(Sa) / q) @ b[a] for a, Sa in enumerate(S)]
        M = cp.Variable((3, 3), hermitian=True)
        cons += [M >> 0, M[0, 0] == d[0], M[1, 1] == d[1], M[2, 2] == d[2],
                 M[0, 1] == f[0], M[1, 2] == f[1], M[2, 0] == f[2]]
    obj = q * (cp.sum(d) + 2 * sum(cp.sum(x) for x in b))
    p = cp.Problem(cp.Maximize(obj), cons)
    p.solve(solver="CLARABEL")
    return p.value


if __name__ == "__main__":
    print(ftheta(19, ((0, 1, 8), (0, 2, 5), (2, 8, 12))),
          ftheta(11, ((0, 1, 2, 5), (0, 1, 4), (1, 3, 4))))

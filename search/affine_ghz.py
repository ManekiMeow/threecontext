# Three-context GHZ graphs from affine planes and nets (RESULTS.md 1.8), and the qutrit picture of H33 (1.9).
# Net of order q with r parallel classes: q^2 points, r*q lines.  Graph H(q,r): contexts A = points, B and C = two
# copies of the lines; a point and a line are adjacent iff incident, a B-line and a C-line iff distinct and parallel.
# Triangle-free because a point never lies on two parallel lines.  Gram matrix: s = 1/sqrt(rq) on incidences,
# t = 1/(q-1) on parallel pairs (forced by expanding psi in each context).  PSD iff r(q-2) >= 2(q-1).
import itertools, numpy as np, networkx as nx

def field(q):
    if q == 4:                                        # GF(4) = {0,1,a,a+1}, a^2 = a+1, addition = xor
        mul = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]]
        return list(range(4)), (lambda x, y: x ^ y), (lambda x, y: mul[x][y])
    return list(range(q)), (lambda x, y: (x + y) % q), (lambda x, y: (x * y) % q)

def net(q, r):
    F, add, mul = field(q); pts = [(x, y) for x in F for y in F]
    dirs = ["v"] + F                                   # vertical, then slopes 0, 1, ...
    lines = []
    for m in dirs[:r]:
        for b in F:
            L = frozenset((b, y) for y in F) if m == "v" else frozenset((x, add(mul(m, x), b)) for x in F)
            lines.append((m, L))
    return pts, lines

def gram(q, r):
    pts, lines = net(q, r); s, t = 1 / np.sqrt(r * q), 1 / (q - 1)
    nP, nL = len(pts), len(lines); n = nP + 2 * nL
    G = np.eye(n); H = nx.Graph(); H.add_nodes_from(range(n))
    for i, p in enumerate(pts):
        for j, (m, L) in enumerate(lines):
            if p in L:
                for off in (nP, nP + nL):
                    G[i, off + j] = G[off + j, i] = s; H.add_edge(i, off + j)
    for j, (m, L) in enumerate(lines):
        for k, (m2, L2) in enumerate(lines):
            if m == m2 and j != k:
                G[nP + j, nP + nL + k] = G[nP + nL + k, nP + j] = t; H.add_edge(nP + j, nP + nL + k)
    return G, H, nP, nL

def check(q, r):
    G, H, nP, nL = gram(q, r); n = len(G)
    tri = sum(1 for c in nx.enumerate_all_cliques(H) if len(c) == 3) if n < 120 else None
    ev = np.linalg.eigvalsh(G); psd = ev[0] > -1e-9; d = int(np.sum(ev > 1e-9))
    yA = np.r_[np.full(nP, 1 / q), np.zeros(2 * nL)]
    yB = np.r_[np.zeros(nP), np.full(nL, 1 / np.sqrt(r * q)), np.zeros(nL)]
    yC = np.r_[np.zeros(nP + nL), np.full(nL, 1 / np.sqrt(r * q))]
    ok = psd and np.allclose(G @ (yA - yB), 0, atol=1e-9) and np.allclose(G @ (yA - yC), 0, atol=1e-9) \
        and np.isclose(yA @ G @ yA, 1)
    pred = r * (q - 2) >= 2 * (q - 1)
    form = r + (q - 1) * (2 * r + q + 1) - (r * (q - 1) if r * (q - 2) == 2 * (q - 1) else 0)
    return n, (nP, nL, nL), tri, round(ev[0], 6), d, ok, pred, form

if __name__ == "__main__":
    print(" q  r |   n  contexts      | triangles | min eig | dim | GHZ realization | predicted PSD | predicted dim")
    for q, r in [(2, 3), (3, 3), (3, 4), (4, 2), (4, 3), (4, 4), (4, 5), (5, 2), (5, 3), (5, 6), (7, 3), (7, 8)]:
        n, ctx, tri, mn, d, ok, pred, form = check(q, r)
        print(f"{q:2d} {r:2d} | {n:4d}  {str(ctx):14s}| {str(tri):9s} | {mn:7.3f} | {d:3d} | {str(ok):15s} | {str(pred):13s} | {form if pred else '-'}")

def qutrit_picture(q=3):
    """H33's points and lines inside the 9-dimensional operator space of a qutrit (RESULTS.md 1.9).
    Phase-point operators A_a = D_a P D_a^dag (P = parity), Hesse SIC |f_a> with A_a = I - 2|f_a><f_a|, and the
    12 MUB projectors P_l = (1/3) sum_{a in l} A_a.  Checks Wootters' line sums, the Hesse orthogonalities used by
    Bengtsson-Blanchfield-Cabello, and that e_a -> A_a/sqrt3, (S1-part of a line) -> P_l/2, psi -> I/sqrt3
    reproduces the inner products of H33 inside S1."""
    w = np.exp(2j * np.pi / q); tau = -np.exp(1j * np.pi / q)          # tau^2 = w (odd prime q)
    X = np.roll(np.eye(q), 1, axis=0); Z = np.diag([w ** k for k in range(q)])
    D = lambda a: tau ** (a[0] * a[1]) * np.linalg.matrix_power(X, a[0]) @ np.linalg.matrix_power(Z, a[1])
    Par = np.zeros((q, q)); Par[[(-k) % q for k in range(q)], list(range(q))] = 1   # |x> -> |-x>
    pts = [(x, y) for x in range(q) for y in range(q)]
    A = {a: D(a) @ Par @ D(a).conj().T for a in pts}
    _, lines = net(q, q + 1)
    hs = lambda M, N: np.trace(M.conj().T @ N)
    out = {"tr(A_a A_b) = q delta": all(np.isclose(hs(A[a], A[b]), q * (a == b)) for a in pts for b in pts),
           "sum_a A_a = q I": np.allclose(sum(A.values()), q * np.eye(q))}
    Pl = [sum(A[a] for a in L) / q for _, L in lines]
    out["P_l rank-1 projectors"] = all(np.allclose(P @ P, P) and np.isclose(np.trace(P), 1) for P in Pl)
    out["tr(A_a P_l) = [a in l]"] = all(np.isclose(hs(A[a], P), 1.0 * (a in L)) for a in pts for P, (_, L) in zip(Pl, lines))
    if q == 3:
        f = {}
        for a in pts:
            val, vec = np.linalg.eigh((np.eye(3) - A[a]) / 2); f[a] = vec[:, -1]
            assert np.isclose(val[-1], 1) and np.allclose(val[:-1], 0)    # (I - A_a)/2 is a rank-1 projector
        out["Hesse SIC |<f_a|f_b>|^2 = 1/4"] = all(np.isclose(abs(np.vdot(f[a], f[b])) ** 2, 0.25) for a in pts for b in pts if a != b)
        m = [np.linalg.eigh(P)[1][:, -1] for P in Pl]
        out["<f_a|m_l> = 0 iff a in l (Hesse)"] = all(np.isclose(abs(np.vdot(f[a], ml)), 0) == (a in L) for a in pts for ml, (_, L) in zip(m, lines))
    s_ = 1 / np.sqrt(q * (q + 1))                                         # point-line overlap in H(q, q+1)
    # S1-parts: points e_a themselves; lines Pi_1 g_l = s sum_{a in l} e_a  <->  e_a -> A_a/sqrt(q), line -> s sqrt(q) P_l
    img = [A[a] / np.sqrt(q) for a in pts] + [s_ * np.sqrt(q) * P for P in Pl]
    S1 = np.array([[hs(x, y) for y in img] for x in img])
    line_proj = np.array([[s_ if a in L else 0 for a in pts] for _, L in lines])
    ref = np.block([[np.eye(len(pts)), line_proj.T], [line_proj, line_proj @ line_proj.T]])
    out["S1-part of H(q,q+1) = (A_a/sqrt q, s sqrt q P_l)"] = np.allclose(S1, ref)
    out["psi -> I/sqrt q: overlaps 1/q and s"] = all(np.isclose(hs(np.eye(q) / np.sqrt(q), x), v) for x, v in
                                                       zip(img, [1 / q] * len(pts) + [s_] * len(Pl)))
    return out

if __name__ == "__main__":
    for q in (3, 5, 7):
        print(f"qudit picture, q = {q}" + (" (H33):" if q == 3 else ":"))
        for k, v in qutrit_picture(q).items(): print(f"   {k}: {v}")

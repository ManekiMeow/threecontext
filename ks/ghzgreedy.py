# Combinatorial search for GHZ states of three completed contexts B_1, B_2, B_3 (see ghzstates.py):
# start from a random state in the complement of the current family, then repeatedly
#   - find triangles (pairwise non-orthogonal transversal triples) inside the current support,
#   - zero the ray of smallest overlap in a randomly chosen triangle (project psi onto its complement),
# until the support is triangle-free (GHZ state found) or psi vanishes.  Many randomised restarts.
import sys, numpy as np
from ghzstates import perp, triangle_tensors, ghz_check

def greedy_ghz(B, fam, rng, E, max_iter=400):
    d = B[0].shape[1]; E12, E23, E13 = E
    Z = [np.array(fam).conj()] if fam else []
    psi = rng.standard_normal(d) + 1j * rng.standard_normal(d)
    zeroed = set()
    for it in range(max_iter):
        V = perp(np.vstack(Z) .conj() if Z else np.zeros((0, d)), d) if Z else np.eye(d, dtype=complex)
        if len(V) == 0: return None
        psi = V.T @ (V.conj() @ psi); n = np.linalg.norm(psi)
        if n < 1e-9: return None
        psi /= n
        o = [np.abs(b.conj() @ psi) for b in B]
        s = [x > 1e-7 for x in o]
        # triangles inside the support
        T = (E12 * s[0][:, None] * s[1][None, :])
        cand = np.argwhere((T[:, :, None] * (E23 * s[2][None, :])[None, :, :] * E13[:, None, :]) > 0)
        if len(cand) == 0: return psi
        a, b, c = cand[rng.integers(len(cand))]
        trio = [(o[0][a], 0, a), (o[1][b], 1, b), (o[2][c], 2, c)]
        trio.sort(key=lambda t: t[0] * (0.5 + rng.random()))       # mostly the smallest overlap, with noise
        _, k, i = trio[0]
        Z.append(B[k][i][None, :].conj())
    return None

def family(B, psi0, rng, restarts=300):
    E = triangle_tensors(B); d = B[0].shape[1]
    fam = [psi0]
    while len(fam) < d:
        found = None
        for r in range(restarts):
            psi = greedy_ghz(B, fam, rng, E)
            if psi is not None and ghz_check(psi, B, *E)[0]:
                found = psi; break
        if found is None: break
        fam.append(found)
        print(f"    GHZ state #{len(fam)}: support sizes {ghz_check(found, B, *E)[1]}", flush=True)
    return len(fam)

if __name__ == "__main__":
    sys.path.insert(0, "../search")
    from ghzstates import realisation, completions
    from verify_all import tricayley_z2z6
    S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
    V, E = tricayley_z2z6(S36); idx = {v: i for i, v in enumerate(V)}
    P = np.eye(36)
    for u, v in E: P[idx[u], idx[v]] = P[idx[v], idx[u]] = 0.25
    cls = [v[0] for v in V]
    X, psi0 = realisation(P, cls); d = X.shape[1]
    B = completions(X, cls, "random", np.random.default_rng(1))
    Z = np.zeros((d, d), complex)
    for copies in (2, 3):
        Bc = [np.block([[b if i == j else Z for j in range(copies)] for i in range(copies)]) for b in B]
        psi = np.concatenate([psi0] + [np.zeros(d)] * (copies - 1))
        print(f"control {copies} x H36 (d={copies*d}): family found =", family(Bc, psi, np.random.default_rng(0)), f"(expected >= {copies})", flush=True)

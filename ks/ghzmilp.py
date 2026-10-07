# Exact test: does a GHZ state orthogonal to a given family exist for three completed contexts B_1, B_2, B_3?
# MILP over psi in C^d (real and imaginary parts) and binaries z_r (r in B_1 ∪ B_2 ∪ B_3):
#   z_r = 1  =>  <r|psi> = 0                    (|Re|, |Im| <= M (1 - z_r))
#   z_a + z_b + z_c >= 1 for every triangle      (support is triangle-free)
#   <f|psi> = 0 for every f in the family
#   <r*|psi> = 1 for a fixed r* in B_1, z_{r*} = 0   (normalisation; enumerate r*)
# If psi is a unit GHZ state, its largest overlap with B_1 is >= 1/sqrt(d), so scaling to <r*|psi> = 1 for that r*
# gives |psi| <= sqrt(d) and every |<r|psi>| <= sqrt(d): the big-M bound M = sqrt(d) + 1 is safe.
import sys, numpy as np, scipy.sparse as sp
from scipy.optimize import milp, LinearConstraint, Bounds
from ghzstates import triangle_tensors, ghz_check

def solve(B, fam, rstar, time_limit=300):
    d = B[0].shape[1]; R = np.vstack(B); n = len(R); M = np.sqrt(d) + 1
    E12, E23, E13 = triangle_tensors(B)
    tris = np.argwhere((E12[:, :, None] * E23[None, :, :] * E13[:, None, :]) > 0)
    # variables: x = [Re psi (d), Im psi (d), z (n)]
    # <r|psi> = sum conj(r_j) psi_j  ->  Re = a.x_re - b.x_im... with r = a + ib: conj(r) psi = (a - ib)(x + iy)
    A_re = np.hstack([R.real, R.imag]); A_im = np.hstack([-R.imag, R.real])   # Re<r|psi>, Im<r|psi>
    rows, lo, hi = [], [], []
    def add(row, l, h): rows.append(row); lo.append(l); hi.append(h)
    Zn = np.zeros(n)
    for i in range(n):
        e = np.zeros(n); e[i] = M
        add(np.concatenate([A_re[i], e]), -np.inf, M); add(np.concatenate([-A_re[i], e]), -np.inf, M)
        add(np.concatenate([A_im[i], e]), -np.inf, M); add(np.concatenate([-A_im[i], e]), -np.inf, M)
    off = [0, d, 2 * d]
    for a, b, c in tris:
        e = np.zeros(n); e[off[0] + a] = e[off[1] + b] = e[off[2] + c] = 1
        add(np.concatenate([np.zeros(2 * d), e]), 1, np.inf)
    for f in fam:
        fr = np.hstack([f.real, f.imag]); fi = np.hstack([-f.imag, f.real])
        add(np.concatenate([fr, Zn]), 0, 0); add(np.concatenate([fi, Zn]), 0, 0)
    add(np.concatenate([A_re[rstar], Zn]), 1, 1); add(np.concatenate([A_im[rstar], Zn]), 0, 0)
    A = sp.csr_matrix(np.array(rows))
    lb = np.concatenate([-M * np.ones(2 * d), np.zeros(n)]); ub = np.concatenate([M * np.ones(2 * d), np.ones(n)])
    ub[2 * d + rstar] = 0
    integ = np.concatenate([np.zeros(2 * d), np.ones(n)])
    res = milp(c=np.zeros(2 * d + n), constraints=LinearConstraint(A, lo, hi), bounds=Bounds(lb, ub),
               integrality=integ, options={"time_limit": time_limit, "presolve": True})
    if res.x is None: return res.status, None
    psi = res.x[:d] + 1j * res.x[d:2 * d]
    return res.status, psi / np.linalg.norm(psi)

def extra_ghz_state(B, fam, time_limit=300, verbose=True):
    """Return (found psi | None, list of statuses).  Status 0 = optimal/feasible, 2 = infeasible, 1 = time limit."""
    d = B[0].shape[1]; E = triangle_tensors(B); stats = []
    for rstar in range(d):
        st, psi = solve(B, fam, rstar, time_limit)
        stats.append(st)
        if psi is not None:
            ok, sizes = ghz_check(psi, B, *E)
            if verbose: print(f"      r*={rstar}: feasible, GHZ check {ok}, support sizes {sizes}", flush=True)
            if ok: return psi, stats
    return None, stats

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
    Bc = [np.block([[b, Z], [Z, b]]) for b in B]
    psi = np.concatenate([psi0, np.zeros(d)])
    # put the B_1 rays of the second block first so the enumeration over r* reaches them quickly
    order = list(range(d, 2 * d)) + list(range(d))
    Bc = [b[order] if k == 0 else b for k, b in enumerate(Bc)]
    found, stats = extra_ghz_state(Bc, [psi])
    print("control 2 x H36: extra GHZ state found:", found is not None, "statuses:", stats[:5])

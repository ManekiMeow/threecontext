# Variant: B_4 = {psi} ∪ ONB(S_1 ∩ psi^perp) ∪ R·C_1 with R a random rotation inside S_1^perp,
# so that S_1 / S_1^perp are linked between B_1 and B_4 with *different* spanning sets
# (creates the internal bases I_1 ∪ R C_1 and {psi} ∪ A_1 ∪ C_1).  Repeat over random seeds,
# also linking S_2^perp ∩ psi... where dimensions allow.
import numpy as np, sys
from ghz_gadget import build, ks_analysis, perp, onb
from ghzgeom import realisation_from_gram
sys.path.insert(0, "../search")
from verify_all import tricayley_z2z6
S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
Vv, E = tricayley_z2z6(S36); idx = {v: i for i, v in enumerate(Vv)}
P = np.eye(36)
for u, v in E: P[idx[u], idx[v]] = P[idx[v], idx[u]] = 0.25
cls = [v[0] for v in Vv]
X, psi = realisation_from_gram(P, cls)
rng = np.random.default_rng(0)
for trial in range(5):
    bases = build(X, psi, cls)
    d = X.shape[1]
    C1 = bases[0][12:]                                  # completion of B_1 (ONB of S_1^perp)
    Q, _ = np.linalg.qr(rng.standard_normal((len(C1), len(C1))))
    if trial % 2: Q = np.eye(len(C1)); Q[:2, :2] = [[np.cos(.7), -np.sin(.7)], [np.sin(.7), np.cos(.7)]]  # minimal rotation
    A1 = perp(np.vstack([psi, C1]), d)
    bases[3] = np.vstack([psi, A1, Q @ C1])
    n, nctx, sat, model, V = ks_analysis(bases)
    print("trial", trial, "rays", n, "internal bases", nctx, "KS-colourable:", sat, flush=True)

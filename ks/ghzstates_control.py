# Positive control for ghzstates.find_family: two block-diagonal copies of the completed H36 contexts (d = 52)
# have at least two orthogonal GHZ states, psi0 (+) 0 and 0 (+) psi0.
import sys, numpy as np
sys.path.insert(0, "../search")
from ghzstates import realisation, completions, find_family
from verify_all import tricayley_z2z6
S36 = (((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 0), (0, 3), (0, 4), (1, 4)), ((0, 1), (1, 0), (1, 1), (1, 3)))
V, E = tricayley_z2z6(S36); idx = {v: i for i, v in enumerate(V)}
P = np.eye(36)
for u, v in E: P[idx[u], idx[v]] = P[idx[v], idx[u]] = 0.25
cls = [v[0] for v in V]
X, psi0 = realisation(P, cls); d = X.shape[1]
B = completions(X, cls, "random", np.random.default_rng(1))
Z = np.zeros((d, d), complex)
B2 = [np.block([[b, Z], [Z, b]]) for b in B]                       # rows: 2d vectors in C^{2d}
psi = np.concatenate([psi0, np.zeros(d)])
print("control (2 x H36, d=52): largest orthogonal GHZ family found =", find_family(B2, psi, restarts=40), "(expected >= 2)")

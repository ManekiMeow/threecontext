# Necessary condition for embedding H36 covariantly (Z_D) as the neighbourhood of a ray in a four-basis KS set:
# a unimodular function on Z_D whose Fourier transform has constant modulus on a 12-element support T and
# vanishes elsewhere  <=>  phases c_t (t in T) with zero periodic autocorrelation at every nonzero shift.
import sys, itertools, numpy as np, torch
torch.set_default_dtype(torch.float64)
rng = np.random.default_rng(0)
def best_residual(D, T, restarts=6):
    T = np.array(sorted(T)); best = 9
    for r in range(restarts):
        th = torch.tensor(rng.uniform(0, 2 * np.pi, len(T)), requires_grad=True)
        opt = torch.optim.LBFGS([th], max_iter=2000, line_search_fn="strong_wolfe", tolerance_grad=1e-14, tolerance_change=1e-16)
        F = torch.exp(-2j * np.pi * torch.outer(torch.arange(D, dtype=torch.float64), torch.tensor(T, dtype=torch.float64)) / D)
        def loss():
            f = F @ torch.exp(1j * th) / np.sqrt(len(T))     # f(x), should be unimodular
            return ((f.abs() ** 2 - 1) ** 2).sum()
        def cl():
            opt.zero_grad(); L = loss(); L.backward(); return L
        for _ in range(3): opt.step(cl)
        best = min(best, loss().item())
    return best
for D in [26, 27, 28, 30, 32, 36, 39, 48]:
    res = []
    for trial in range(15):
        T = rng.choice(D, 12, replace=False)
        res.append(best_residual(D, T, 3))
    cos = None
    if D % 12 == 0:
        T = [i * (D // 12) for i in range(12)]           # subgroup of order 12 (coset case)
        cos = best_residual(D, T, 3)
    print("D=%d  random 12-sets: min residual %.2e  (median %.2e)   subgroup support: %s"
          % (D, min(res), np.median(res), "%.1e" % cos if cos is not None else "n/a"), flush=True)

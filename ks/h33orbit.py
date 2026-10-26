# Orbit of a ray under Gamma = <AGL(2,3), W> by breadth-first search on rays (vectors up to phase).
import sys, numpy as np
import h33model as M
from h33group import gens

def rkey(v):
    j = int(np.argmax(np.abs(v) > 1e-7)); z = v[j]; w = v * abs(z) / z
    return tuple(np.round(np.concatenate([w.real, w.imag]), 6) + 0.0)

def orbit(v0, gens, limit=10 ** 6):
    seen = {rkey(v0): v0}; frontier = [v0]
    while frontier:
        nxt = []
        for v in frontier:
            for g in gens:
                w = g @ v; k = rkey(w)
                if k not in seen:
                    seen[k] = w; nxt.append(w)
                    if len(seen) > limit: return list(seen.values()), False
        frontier = nxt
    return list(seen.values()), True

if __name__ == "__main__":
    orb, ok = orbit(M.psi, gens)
    print("orbit of psi:", len(orb), "rays, closed:", ok, flush=True)
    np.save("h33_orbit_psi.npy", np.array(orb))
    for name, v in [("B-line", M.B[0]), ("beta", M.beta[M.nz[0]])]:
        o, ok = orbit(v, gens); print(f"orbit of a {name}:", len(o), "closed:", ok, flush=True)

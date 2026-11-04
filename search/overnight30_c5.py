"""Exact C5 symmetry reduction of the full theta SDP on (10,10,10) graphs.

There are six consecutive vertex orbits of size five. With real invariant
entries B[(a,i),(b,j)] = f_ab(j-i), the five Fourier blocks are
F_k[a,b] = sum_s f_ab(s) exp(2*pi*i*k*s/5). Only F0,F1,F2 need PSD constraints;
the other two are their complex conjugates. Trace(B)=5*sum_a f_aa(0) and
sum(B)=5*sum(F0). No event-probability or overlap-weight ansatz is imposed.

All numerical answers reconstruct and audit the entire real 30x30 primal
matrix. Its numerical rank is not a claim about minimum physical dimension.
"""
import math
import time

ORBIT_PAIRS = [(a, b) for a in range(6) for b in range(a+1, 6) if a//2 != b//2]


def support_shifts(edges):
    """Return the 60 allowed shifts, or None for a non-C5-double support."""
    clean = {tuple(sorted(edge)) for edge in edges}
    if any(u < 0 or v >= 30 or u//10 == v//10 for u, v in clean):
        return None
    selected = [(5*a, 5*b+s) in clean for a, b in ORBIT_PAIRS for s in range(5)]
    recovered = {(5*a+i, 5*b+(i+s) % 5)
                 for index, (a, b) in enumerate(ORBIT_PAIRS) for s in range(5)
                 if selected[5*index+s] for i in range(5)}
    return selected if recovered == clean else None


def reconstruct(values):
    import numpy as np
    matrix = np.zeros((30, 30))
    for a in range(6):
        for i in range(5):
            matrix[5*a+i, 5*a+i] = values[a]
    for index, (a, b) in enumerate(ORBIT_PAIRS):
        for s in range(5):
            entry = values[6+5*index+s]
            for i in range(5):
                u, v = 5*a+i, 5*b+(i+s) % 5
                matrix[u, v] = matrix[v, u] = entry
    return matrix


def fourier_blocks(values):
    import numpy as np
    result = []
    for k in range(5):
        block = np.diag(values[:6]).astype(complex)
        roots = np.exp(2j*np.pi*k*np.arange(5)/5)
        for index, (a, b) in enumerate(ORBIT_PAIRS):
            entry = np.dot(roots, values[6+5*index:11+5*index])
            block[a, b], block[b, a] = entry, entry.conjugate()
        result.append(block)
    return result


def audit_dft(values):
    """Independent unitary-change-of-basis check, also valid for indefinite inputs."""
    import numpy as np
    matrix = reconstruct(values)
    unitary = np.zeros((30, 30), dtype=complex)
    for a in range(6):
        for i in range(5):
            for k in range(5):
                unitary[5*a+i, 6*k+a] = np.exp(2j*np.pi*k*i/5)/math.sqrt(5)
    diagonalized = unitary.conj().T@matrix@unitary
    expected = np.zeros((30, 30), dtype=complex)
    blocks = fourier_blocks(values)
    for k, block in enumerate(blocks):
        expected[6*k:6*k+6, 6*k:6*k+6] = block
    return {"unitary_residual": float(np.max(abs(unitary.conj().T@unitary-np.eye(30)))),
            "block_diagonalization_residual": float(np.max(abs(diagonalized-expected))),
            "trace_identity_residual": float(abs(np.trace(matrix)-5*sum(values[:6]))),
            "objective_identity_residual": float(abs(np.sum(matrix)-5*np.sum(blocks[0])))}


class C5ThetaSolver:
    """Reusable parameterized Fourier SDP; timeouts remain the caller's job."""
    def __init__(self):
        import cvxpy as cp
        import numpy as np
        self.values = cp.Variable(66)
        self.forbidden = cp.Parameter(60, nonneg=True)
        blocks = []
        for k in range(3):
            entries = [[0 for _ in range(6)] for _ in range(6)]
            for a in range(6):
                entries[a][a] = self.values[a]
            roots = np.ones(5) if k == 0 else np.exp(2j*np.pi*k*np.arange(5)/5)
            for index, (a, b) in enumerate(ORBIT_PAIRS):
                entry = cp.sum(cp.multiply(roots, self.values[6+5*index:11+5*index]))
                entries[a][b] = entry
                entries[b][a] = entry if k == 0 else cp.conj(entry)
            blocks.append(cp.bmat(entries))
        self.problem = cp.Problem(cp.Maximize(5*cp.sum(blocks[0])),
                                  [5*cp.sum(self.values[:6]) == 1,
                                   cp.multiply(self.forbidden, self.values[6:]) == 0,
                                   *(block >> 0 for block in blocks)])
        self.unitary = np.zeros((30, 30), dtype=complex)
        for a in range(6):
            for i in range(5):
                for k in range(5):
                    self.unitary[5*a+i, 6*k+a] = np.exp(2j*np.pi*k*i/5)/math.sqrt(5)

    def dual_candidate(self):
        """Reconstruct an uncertified full dual; exact rounding checks are separate."""
        import numpy as np
        upper = float(self.problem.constraints[0].dual_value)
        raw = [np.asarray(constraint.dual_value) for constraint in self.problem.constraints[-3:]]
        # CVXPY's returned complex duals need not be perfectly Hermitian.
        blocks = [(block+block.conj().T)/2 for block in raw]
        spectral = np.zeros((30, 30), dtype=complex)
        ordered = (blocks[0], blocks[1]/2, blocks[2]/2, blocks[2].conj()/2, blocks[1].conj()/2)
        for k, block in enumerate(ordered):
            spectral[6*k:6*k+6, 6*k:6*k+6] = block
        slack = (self.unitary@spectral@self.unitary.conj().T).real
        dual = upper*np.eye(30)-(slack+slack.T)/2
        return dual, upper

    def solve(self, edges):
        import numpy as np
        selected = support_shifts(edges)
        if selected is None:
            return {"error": "support does not have the required C5-double action"}
        self.forbidden.value = 1-np.asarray(selected, dtype=float)
        mask = np.ones((30, 30))
        np.fill_diagonal(mask, 0)
        for u, v in edges:
            mask[u, v] = mask[v, u] = 0
        attempts, start = [], time.monotonic()
        for solver in ("CLARABEL", "SCS"):
            try:
                if solver == "CLARABEL":
                    self.problem.solve(solver=solver, max_iter=250)
                else:
                    self.problem.solve(solver=solver, eps=1e-6, max_iters=20000)
                if (self.problem.value is None or not np.isfinite(self.problem.value)
                        or self.values.value is None):
                    attempts.append(f"{solver}: missing finite primal result")
                    continue
                matrix = reconstruct(np.asarray(self.values.value))
                if not np.isfinite(matrix).all():
                    attempts.append(f"{solver}: nonfinite primal entries")
                    continue
                spectrum = np.linalg.eigvalsh(matrix)
                trace_residual = float(abs(np.trace(matrix)-1))
                zero_residual = float(np.max(abs(mask*matrix)))
                objective_residual = float(abs(np.sum(matrix)-self.problem.value))
                if (self.problem.status not in ("optimal", "optimal_inaccurate")
                        or spectrum[0] < -2e-5 or trace_residual > 2e-5
                        or zero_residual > 2e-5 or objective_residual > 2e-5
                        or self.problem.value > 3.0001):
                    attempts.append(f"{solver}: unreliable full primal residuals/status")
                    continue
                result = {"theta_numeric": float(self.problem.value), "solver": solver,
                        "solver_status": self.problem.status, "reduction": "C5 Fourier: real6 + complex6 + complex6",
                        "trace_residual": trace_residual, "zero_residual": zero_residual,
                        "objective_residual": objective_residual, "min_eigenvalue": float(spectrum[0]),
                        "primal_rank_numeric_1e-7": int(sum(spectrum > 1e-7)),
                        "primal_diagonal": matrix.diagonal().tolist(), "primal_matrix": matrix.tolist(),
                        "numerical_seconds": time.monotonic()-start}
                try:
                    dual, upper = self.dual_candidate()
                    if np.isfinite(dual).all() and np.isfinite(upper):
                        result["theta_dual_candidate_matrix"] = dual.tolist()
                        result["theta_dual_upper_numeric"] = upper
                        result["theta_dual_fixed_residual"] = float(max(
                            [abs(dual[v, v]-1) for v in range(30)]+[abs(dual[u, v]-1) for u, v in edges]))
                except Exception as error:
                    result["dual_reconstruction_error"] = f"{type(error).__name__}: {error}"
                return result
            except Exception as error:
                attempts.append(f"{solver}: {type(error).__name__}: {error}")
        return {"error": "; ".join(attempts)}


def self_test():
    """Compare five fixed supports with an independent, unreduced theta SDP."""
    import cvxpy as cp
    import networkx as nx
    import numpy as np
    fixtures = [
        "]???????@mB[a[pKKq?G??_?@??@???_@tOINg@Su_cXk@GwkD?aINOWBboBGK]?K_to?X@z??",
        "]???????CFgmHkDs@{BC?KO?W_CO_EG??@A??GO??`??`???`??[AG_F@_GCoK_?q?q?B_@c??",
        "]???????EGKOKOEGPABS?kO@W_DO_Eg??@kR?LeOGew?oLo?qHoRka?SyG_An@C?I}?OATsO??",
        "]???????Amd{c{R[SuE]_xw@ro@z_C|_AK?Co_?bA?AECA?ECA?GAS?A?d??OSG?_@O_@?A`??",
        "]???????FQ]C]Cn?ZO?K_@o?B_?B_?H_?BF??Ww?@r??bc?@b_?___oGGIG@@@o?CEB??KCE??",
    ]
    rng = np.random.default_rng(31)
    for _ in range(20):
        if max(audit_dft(rng.normal(size=66)).values()) >= 1e-11:
            raise AssertionError("Fourier unitary/trace/objective identity failed")
    matrix = cp.Variable((30, 30), symmetric=True)
    forbidden = cp.Parameter((30, 30))
    full = cp.Problem(cp.Maximize(cp.sum(matrix)),
                      [matrix >> 0, cp.trace(matrix) == 1, cp.multiply(forbidden, matrix) == 0])
    reduced, differences = C5ThetaSolver(), []
    for encoded in fixtures:
        graph = nx.from_graph6_bytes(encoded.encode())
        edges = list(graph.edges())
        if support_shifts(edges[:-1]) is not None:
            raise AssertionError("noninvariant support was accepted")
        zeros = np.ones((30, 30))
        np.fill_diagonal(zeros, 0)
        for u, v in edges:
            zeros[u, v] = zeros[v, u] = 0
        forbidden.value = zeros
        full.solve(solver="CLARABEL", max_iter=250)
        answer = reduced.solve(edges)
        if "error" in answer or full.value is None:
            raise AssertionError(answer)
        difference = abs(answer["theta_numeric"]-full.value)
        if difference >= 3e-6:
            raise AssertionError(f"reduced/full theta discrepancy: {difference}")
        differences.append(difference)
    permutation = {5*a+i: 5*(a ^ 1)+(2*i+3*a) % 5 for a in range(6) for i in range(5)}
    relabelled = [(permutation[u], permutation[v]) for u, v in edges]
    changed = reduced.solve(relabelled)
    if "error" in changed or abs(changed["theta_numeric"]-answer["theta_numeric"]) >= 3e-6:
        raise AssertionError("cyclic-origin/unit/orbit relabelling invariance failed")
    print("C5 Fourier controls passed: 20 DFT audits, five independent full-SDP comparisons,",
          "noninvariant-support rejection, relabelling; max theta difference", max(differences))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    arguments = parser.parse_args()
    if arguments.self_test:
        self_test()
    else:
        parser.error("provide --self-test, or import C5ThetaSolver")

"""Global Hall condition for two contexts of a thirty-event GHZ support.

For S in context A, let N=N_B(S) be its allowed support neighbors in B.
The vectors within each context are orthonormal, and psi lies in span(B).
All vectors in S are orthogonal to B minus N, even when some allowed edges
have zero overlap. Therefore P_S psi=P_S P_N psi, and contractivity gives

    c(S)=||P_S psi||^2 <= ||P_N psi||^2=c(N).

Equivalently, S together with B minus N is orthonormal, so Bessel gives
c(S)+c(B minus N)<=1. This is the weighted two-context Hall inequality used
by the transportation argument in RESULTS section3.4.

Now both contexts have10 events, every probability is positive and strictly
less than1/9. If1<=s=|S|<=9 then c(S)=1-c(A minus S)>(s-1)/9. A deficient
N has|N|<=s-1 and c(N)<=(s-1)/9, a contradiction (including N empty).
If s10, then c(S)=1 while any deficient N has at most9 terms and c(N)<1.
The empty S is vacuous. Thus all Hall inequalities hold, and every pair
support has a perfect matching, in every dimension. Pair connectivity is
still a separate dimension19 condition and is not claimed globally.

The stored control is a saved triangle-free30-event support whose contexts
pass the earlier contraction/component/bridge screens, but whose B--C pair
has a Hall-deficient subset of size8. This script checks the new global
certificate, historical pair_hall_d19 compatibility, malformed witnesses,
the fast-screen integration, and structural dominance cut replay. It uses
only standard Python and does not invoke SAT, SDP, or NetworkX.
"""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path


def graph6(graph):
    """Minimal independent headerless encoder for the30-vertex control."""
    assert graph['n'] == 30
    edges = {tuple(edge) for edge in graph['edges']}
    bits = [int((u, v) in edges) for v in range(1, 30) for u in range(v)]
    bits += [0]*(-len(bits) % 6)
    return chr(93)+''.join(chr(63+sum(bits[i+j] << (5-j) for j in range(6)))
                          for i in range(0, len(bits), 6))


def main():
    from overnight30_screen import adjacency, fast_screen, hall_deficiency, verify_certificate
    from overnight30_symmetry import dominance_clause, replay_structural_dominance
    inequalities = 0
    for s in range(1, 11):
        for n in range(s):
            lower = F(s-1, 9) if s < 10 else F(1)
            upper = F(n, 9)
            assert lower >= upper and (s < 10 or n > 0 or lower > upper)
            inequalities += 1
    control = json.loads(Path(__file__).with_name('event30_hall_control.json').read_text())
    assert control['method'] == 'pair_hall'
    assert verify_certificate(control)
    legacy = deepcopy(control)
    legacy['method'] = 'pair_hall_d19'
    assert verify_certificate(legacy)
    graph = control['graph']
    adj = adjacency(graph)
    k, l = control['witness']['pair']
    S = set(control['witness']['deficient_subset'])
    N = set().union(*(adj[v] for v in S)) & set(graph['parts'][l])
    assert len(S) == 8 and len(N) == 7
    discovered = hall_deficiency(adj, graph['parts'][k], graph['parts'][l])
    assert discovered == sorted(S)
    screened = fast_screen(graph)
    assert screened['status'] == 'excluded' and screened['certificate'] == control
    bad = deepcopy(control)
    bad['witness']['deficient_subset'] = [min(S)]
    try:
        verify_certificate(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('A nondeficient singleton was accepted')
    orbits = [[(u, v)] for u in range(30) for v in range(u+1, 30) if u//10 != v//10]
    cut = {'clause': dominance_clause(graph, orbits),
           'reason': 'verified_structural_dominance',
           'witness': {'graph6': graph6(graph), 'method': 'pair_hall',
                       'witness': control['witness']}}
    assert replay_structural_dominance(cut, orbits) == cut['clause']
    bad_cut = deepcopy(cut)
    bad_cut['clause'] = bad_cut['clause'][1:]
    try:
        replay_structural_dominance(bad_cut, orbits)
    except ValueError:
        pass
    else:
        raise AssertionError('A modified dominance clause was accepted')
    print(json.dumps({'scope': 'all dimensions for30positive events in10+10+10 contexts',
                      'cardinality_boundary_cases': inequalities,
                      'hall_subset_size': len(S), 'neighbor_size': len(N),
                      'legacy_certificate_verified': True, 'fast_screen_global_exclusion': True,
                      'structural_dominance_verified': True, 'tamper_controls_rejected': 2}, indent=2))


if __name__ == '__main__':
    main()

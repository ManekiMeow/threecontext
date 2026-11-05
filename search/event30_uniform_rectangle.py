"""Exact zero-rectangle obstruction from two uniform ten-event contexts.

MAIN RESULT (all dimensions, complex vectors allowed)
-----------------------------------------------------
In a thirty-event three-context GHZ realization, if every event in A and B
has probability 1/10, their allowed support must contain a 3-by-3 zero
rectangle. Equivalently, some three A events and three B events have no
allowed cross-edge. This statement also holds for the actual nonzero support.
It is conditional on the two uniform contexts, not an unrestricted graph
exclusion.

Suppose no such rectangle exists. Triangle-freeness then forces every C event
to have at most two allowed neighbors in A or at most two in B. The associated
neighbor mass alpha is at most 1/5. The one-sector bound

    c <= h(alpha),   h(a)=a(1-a)/(1+3a),

and strict increase of h on [0,1/5] give c<=1/10. Since C has ten events of
total probability one, all are uniform too. Equality forces alpha=1/5 and
the two other sector masses beta=rho=2/5. Each C event therefore has degree
type (2,4) or (4,2) towards A,B, and saturates the three-sector inequality.

For a type (2,4) event, write S=N_A(c), T=N_B(c). Choose event phases so all
state overlaps are positive. Let psi_S,psi_T be the state projections onto
these neighbor subsets and psi_R=psi-psi_S-psi_T. The three vectors are
orthogonal. Cauchy equality gives

    g_c = sqrt(c) (psi_S/alpha + psi_T/beta - psi_R/rho)
        = (3/4) sum_{a in S} g_a + (1/2) sum_{b in T} g_b
          - (1/4) sum_{a in A} g_a.

Every incident A overlap is 1/2 and every incident B overlap is 1/4. The
opposite type has coefficients (1/2,3/4,-1/4) and overlaps (1/4,1/2).
These are exact complex-vector equality statements, not assumptions that
all quantum overlaps are real.

For two distinct events with the same type (2,4), put
s=|S intersect S'| and t=|T intersect T'|. Their orthogonality reads

    3*s+t=2,

so their two-element A neighborhoods are disjoint. For opposite types,
expanding the first event gives 3*s+4*t=4 and expanding the second gives
4*s+3*t=4. Thus s=t=4/7, impossible for integers. All ten C events must
therefore have the same type, but ten pairwise disjoint two-element subsets
cannot lie in a ten-event context. This proves the main result without any
moment-bound, dimension, graph automorphism or numerical assumption.

UNIFORM SATURATED-TYPE COROLLARIES
--------------------------------
In a fully uniform thirty-event realization, the saturated degree types are
(2,4),(4,2),(4,4). Opposite types (2,4),(4,2) cannot coexist in one context.
Two (4,4) events would require 3*(s+t)=8; a (2,4) event with a (4,4) event
would require 6*s+3*t=8, and the reflected case is analogous. Hence the
saturated events in a context are all of one degree-two orientation, or
there is exactly one (4,4) event and no degree-two event.

The degree-two events have pairwise disjoint neighborhoods on their small
side. Four such uniform events have total probability 2/5, contrary to the
strict m=4 moment bound of Theorem 53. Thus there are at most THREE of them
in a context. This sharpening uses the verified moment certificate; the main
zero-rectangle theorem above does not depend on it.

ISOCLINIC CLOSED-PETERSEN COROLLARY (all dimensions)
---------------------------------------------------
Suppose the A--B allowed pair support is I+P, where P has Petersen parameters
(10,3,0,1), and their reduced rank-nine spaces are isoclinic with cosine u.
Assume first 0<u<1 and write x=u^2, p_i=sqrt(c_i), q_j=sqrt(d_j). The overlap
matrix satisfies

    XX*=xI+(1-x)pp*,   X*X=xI+(1-x)qq*.

For each column j, its three off-diagonal allowed rows are the neighbors
of j in P. Any two are nonadjacent in P and their allowed row supports
intersect only at j. The off-diagonal row Gram identity forces these entries
to be nonzero and of the form

    X_ij=sqrt(1-x)*p_i*eta_j,   |eta_j|=1.

For a fixed row i, any two of its three leaf-neighbor columns j,k likewise
intersect only at i. The column Gram identity gives c_i=sqrt(d_j*d_k).
Using all three pairs of leaf-neighbor columns forces d_j=d_k=d_l=c_i.
The graph linking these equalities is connected, so every A and B probability
is 1/10. Equivalently, every two nonadjacent P vertices have equal B weights,
and the six-regular complement of P is connected. This reasoning permits
allowed zero overlaps: the off-diagonal ones used in the proof are forced
nonzero by the Gram identities.

The endpoint u=0 gives X=pq* with full nonzero support, incompatible with
I+P. If u=1, the two context spans coincide, which is impossible in a
three-context GHZ realization: for any C event, orthogonality of its A and B
neighbor projections gives <g_C,E_A E_B g_C>=0, while E_A=E_B and its state
overlap is positive. Thus the interior argument covers every possible u.

Finally I+P has no 3-by-3 zero rectangle. For any three P vertices, let e be
their induced edge count and t their triple closed-neighborhood intersection.
Their closed-neighborhood union has size 9-e+t. There are no triangles, so
e<=2; if e=2 the path center contributes t>=1. The union therefore has size
at least eight. The main theorem now gives a contradiction. In particular,
an isoclinic closed-Petersen pair cannot occur in a thirty-event realization
in ANY ambient dimension, even with initially nonuniform probabilities.

The helper verifies the rational identities, all saturated-type integer
possibilities, the parameter-graph relations, and the optional strict
four-event moment sharpening. It is not a numerical search.
"""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use python -O.')

from fractions import Fraction as F
import itertools
import json
from pathlib import Path


def h(a):
    a = F(a)
    return a*(1-a)/(1+3*a)


def saturated_coefficients(degrees):
    a, b = (F(q, 10) for q in degrees)
    r, c = 1-a-b, F(1, 10)
    assert min(a, b, r) > 0 and c*(1/a+1/b+1/r) == 1
    return c*(1/a+1/r), c*(1/b+1/r), c/r


def overlap(first, second, s, t):
    a, b, constant = saturated_coefficients(first)
    return a*s/second[0]+b*t/second[1]-constant


def verify():
    assert h(F(1, 5)) == F(1, 10)
    assert all(h(F(k, 10)) < F(1, 10) for k in (0, 1))
    # h'(a)=(1+a)(1-3a)/(1+3a)^2 is positive on [0,1/5].
    assert (1+F(1, 5))*(1-3*F(1, 5)) > 0
    types = [(2, 4), (4, 2), (4, 4)]
    assert [saturated_coefficients(t) for t in types] == [
        (F(3, 4), F(1, 2), F(1, 4)),
        (F(1, 2), F(3, 4), F(1, 4)),
        (F(3, 4), F(3, 4), F(1, 2))]
    compatibility = {}
    for first, second in itertools.combinations_with_replacement(types, 2):
        assert overlap(first, first, *first) == 1
        permitted = [(s, t)
                     for s in range(min(first[0], second[0])+1)
                     for t in range(min(first[1], second[1])+1)
                     if overlap(first, second, s, t) == 0
                     and overlap(second, first, s, t) == 0]
        expected = [(0, 2)] if first == second == (2, 4) else (
            [(2, 0)] if first == second == (4, 2) else [])
        assert permitted == expected
        compatibility[str((first, second))] = permitted
    assert 10*2 > 10  # Ten same-orientation disjoint pairs cannot fit.

    vertices = list(itertools.combinations(range(5), 2))
    P = [{j for j, b in enumerate(vertices) if set(a).isdisjoint(b)}
         for a in vertices]
    assert all(len(row) == 3 for row in P)
    closed = [row | {i} for i, row in enumerate(P)]
    for i, j in itertools.combinations(range(10), 2):
        assert len(P[i] & P[j]) == (0 if j in P[i] else 1)
        assert len(closed[i] & closed[j]) == (2 if j in P[i] else 1)
    leaf_pairs = 0
    relations = [set() for _ in range(20)]
    for j in range(10):
        for i, k in itertools.combinations(P[j], 2):
            assert closed[i] & closed[k] == {j}
            leaf_pairs += 1
        for i in P[j]:
            # Consequence of the three positive equations
            # c_i^2=d_j*d_k=d_j*d_l=d_k*d_l: c_i=d_j.
            relations[i].add(10+j)
            relations[10+j].add(i)
    reached, pending = {0}, [0]
    while pending:
        vertex = pending.pop()
        for neighbor in relations[vertex]-reached:
            reached.add(neighbor)
            pending.append(neighbor)
    assert len(reached) == 20
    unions = [len(set().union(*(closed[i] for i in subset)))
              for subset in itertools.combinations(range(10), 3)]
    assert len(unions) == 120 and min(unions) == 8

    from event30_disjoint_probability import verify_bound_certificate, bound_is_strict
    path = Path(__file__).with_name('overnight30_disjoint_reviews')/'2026-11-05T0000'/'disjoint_04_mean_bound.json'
    moment = json.loads(path.read_text())
    assert moment['subset_size'] == 4
    assert 4*verify_bound_certificate(moment) == F(2, 5)
    assert bound_is_strict(moment)
    return {
        'scope': 'two uniform ten-event contexts; every dimension and complex vectors',
        'saturated_type_compatibility': compatibility,
        'rectangle_proof_requires_moment_bound': False,
        'uniform_degree_two_events_per_context_at_most': 3,
        'petersen_leaf_pair_checks': leaf_pairs,
        'petersen_probability_equality_component_size': len(reached),
        'petersen_triples_checked': len(unions),
        'minimum_closed_neighborhood_union': min(unions),
        'isoclinic_petersen_corollary_scope': 'arbitrary dimension; initially nonuniform probabilities allowed',
        'unrestricted_graph_filter': False,
    }


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))

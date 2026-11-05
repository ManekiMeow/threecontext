"""Exact uniform-probability exclusion of an isoclinic reduced pair in dimension 19.

This file supplies proof notes and elementary rational/combinatorial checks.
It does not search realizations or assume real overlaps. The strengthened
conclusion below excludes even ONE reduced isoclinic pair under uniform
probabilities; it makes no such claim with nonuniform event probabilities.

GENERAL ANGLE RESTRICTION (no uniform-probability assumption)
==========================================================
Let A,B,C be the reduced rank-nine projections on psi-perp (dimension18),
and suppose A,B are isoclinic with principal cosine u. The generic case of
RESULTS Lemma33 has 0<u<1; the nongeneric case is excluded by Theorem35.
The Hermitian operator T=(AB+BA)/2 has eigenvalues

    lambda_minus=(u*u-u)/2, lambda_plus=(u*u+u)/2,

each ninefold. Since Re tr(CAB)=-1 and rank C=9, the variational principle
gives -1>=9*lambda_minus. Hence 1/3<=u<=2/3. Equality at either endpoint
forces C to be precisely the negative eigenspace of T. On each two-dimensional
pair block, if a,b are unit A/B vectors with <a,b>=u, that eigenspace is the
line through a-b; its squared overlap with both a and b is (1-u)/2. Thus all
three reduced pairs would be isoclinic, excluded by RESULTS Theorem36.
Consequently every isoclinic reduced pair has 1/3<u<2/3.

UNIFORM-PROBABILITY RESTRICTION
==============================
Now assume all thirty probabilities equal1/10. Work with the ACTUAL nonzero
overlap support; deleting allowed zero-overlap edges preserves a GHZ realization
and triangle-freeness. Let A--B be isoclinic with x=u*u. After positive-state
rephasing its10x10 overlap matrix X has all row and column sums1, and

    XX*=x I+(1-x)J/10,  X*X=x I+(1-x)J/10.

Every row and column therefore has squared norm (1+9*x)/10<1/2. Cauchy--Schwarz
with sum1 forces at least3 nonzero entries. The exact uniform degree list in
RESULTS Proposition37 restricts every cross-degree to2,3,4, so this pair's
degrees are3 or4. Also, every pair of distinct rows (and columns) has a common
neighbor because its Gram entry (1-x)/10 is positive.

The value x=1/6 is impossible: all rows then have squared norm1/4 and at most4
nonzeros, so equality in Cauchy forces exactly4 entries, each1/4 (also over C).
Writing X=M/4 gives off-diagonal entries (MM*)_ij=4/3, an impossible integer.
In fact the same row bound proves x>1/6.

No A or B event can have degree type(4,2) or(4,4), where the first degree is
toward the isoclinic partner. Both types saturate the uniform three-sector
bound, forcing every overlap along that four-neighbor side to1/4. Its row
norm would be1/4, again x=1/6. The other types with cross-degree2 are already
excluded by the isoclinic row bound or the exact uniform degree list.
Thus at A and B the two cross-degrees are3 or4, and they cannot both be4.

INCIDENCE OBSTRUCTION
====================
Suppose, contrary to the following necessary condition, every pair of C events
has a common A neighbor. Let r be the number of A vertices with B-degree4 and
s the number with C-degree4. Counting pairs of B vertices through A vertices
gives30+3*r>=45 (the isoclinic Gram requires every B pair to meet); counting C
pairs under the supposition gives30+3*s>=45. But r+s<=10 by the preceding
paragraph. Hence r=s=5. There are35 A--B edges, so exactly5 B vertices also
have A-degree4. The pair capacities in BOTH directions are exactly45.
Therefore every pair of A points is on exactlyone B line and every pair of
B lines meets at exactlyone A point. All degrees and line sizes are3 or4,
with both sizes present.

Choose a3-line L and a4-line M. They meet once, so their union contains6 points;
choose a point p outside it. For either line K=L or M, the lines through p are
in bijection with the points of K: unique joins give existence and uniqueness,
and every line through p meets K exactlyonce. Thus the degree of p equals both3
and4, a contradiction. Interchanging A,B proves the analogous statement.

Consequently, with ONE isoclinic pair A--B, C must contain two events whose
ACTUAL A-neighborhoods are disjoint, and two (possibly different) events whose
ACTUAL B-neighborhoods are disjoint. This is not a direct filter on an allowed
support supergraph, whose extra zero-overlap edges may destroy disjointness.

If TWO reduced pairs are isoclinic, relabel them A--B and A--C. The second
pair's positive off-diagonal Gram gives every C pair a common A neighbor,
contradicting the preceding result. Therefore two reduced isoclinic pairs
are impossible in dimension19 under uniform probabilities1/10.

ONE-PAIR PROJECTION IDENTITY
===========================
The reduced pair satisfies AB+BA=A+B-(1-x)I on psi-perp. Restoring the full
context projections E_A,E_B therefore gives

    E_A E_B+E_B E_A=E_A+E_B-(1-x)I+(1-x)|psi><psi|.

For each C event, triangle-freeness makes the anticommutator expectation zero.
Writing p_c=||E_A g_c||^2 and q_c=||E_B g_c||^2, it follows that

    p_c+q_c=(1-x)(1-c_c).

Under uniform probabilities, a saturated type(2,4) or(4,2) event has p_c+q_c=3/4,
forcing x=1/6; a type(4,4) event has p_c+q_c=1/2, forcing x=4/9. Both values are
already excluded. Thus every event in all three contexts has degree type
(3,3),(3,4), or(4,3). If x<7/27, every A--B row has squared norm<1/3 and so must
have four nonzero entries. All A/B events then have three neighbors in C, and
the average degree of C into either other context is three. Since those C
degrees are at least three, they all equal three. The pair supports therefore
have regular degrees(4,3,3). The following stronger argument reaches this
regularity without a restriction on x and then excludes it.

COMPLETE UNIFORM ONE-PAIR EXCLUSION
=================================
Write k=(1-x)/10>0. Every row/column of X has sum1, squared norm1-9*k, and
three or four nonzeros. Every pair of rows intersects, since its Gram entry
is k>0. For each column j, define L_j on its incident rows: connect i,l exactly
when those two rows' supports meet ONLY in column j. Each such edge imposes

    X_ij * conjugate(X_lj)=k.

In any connected component of L_j, all entries therefore have a common phase;
their positive magnitudes alternate R,S with RS=k along a bipartition. An odd
cycle forces every magnitude in that component to sqrt(k).

1. LOCAL EDGE BOUNDS. A three-entry column cannot have a triangle L_j: the
   column sum then forces three entries1/3 and k=1/9, contradicting its squared
   norm1-9*k=0. Thus L_j has at most two edges for a three-entry column.

   For a four-entry column, a connected odd L_j would give four entries1/4,
   hence k=1/16 and squared norm1/4, contradicting1-9*k=7/16. A connected
   balanced bipartite L_j would give entries R,R,S,S, with positive R,S and
   2*(R+S)=1. Its squared norm is1/2-4*k; equality to1-9*k forces k=1/10,
   contradicting RS<=1/16. Therefore L_j has at most three edges. With exactly
   three it is either a star K1,3 or a triangle plus an isolated vertex; the
   remaining three-edge possibility is the already excluded path P4.

2. EVERY DEGREE IS FOUR. Let a be the number of three-entry columns. Total
   row-pair incidences are60-3*a. All45 pairs intersect, so the number s of
   singleton-intersection pairs satisfies s>=90-(60-3*a)=30+3*a. But summing
   the preceding local bounds gives s<=3*(10-a)+2*a=30-a. Thus a=0. Every
   column, and hence every row, has degree four. Both inequalities are tight:
   there are30 singleton pairs and15 double pairs, with no larger intersections;
   every L_j is a star or triangle plus isolated. Applying the same reasoning
   to X* gives the corresponding classification for each row.

3. THE COLUMN TYPES CANNOT MIX. A star column has positive real entries R
   at its center and S at its three leaves. The sum and norm equations give
   R+3*S=1, R*S=k, S*S=k/2, hence S=1/5,R=2/5,k=2/25. A triangle column has
   three equal entries z and isolated entry w=1-3*z, with |z|^2=k and
   Re(z)=7*k/2. At k=2/25 these entries are

       z=(7 +/- i)/25,  w=(4 -/+ 3*i)/25,

   all nonreal. Thus at the only possible star parameter, no star entry can
   belong to a triangle-type row or conversely. Since the support is connected
   (indeed every row pair intersects), all rows and columns have the same type.
   At other parameters stars do not occur at all.

4. ALL TRIANGLES GIVE CLOSED PETERSEN SUPPORT. Each row belongs to four
   columns and has six singleton and three double row intersections: its total
   incidences are12=6+2*3. In a triangle column it contributes two singleton
   intersections if it is in the triangle and zero if isolated. Consequently
   every row is isolated in precisely one column. Relabel columns to place
   these isolated entries on the diagonal. The other support entries form a
   three-regular zero-diagonal0/1 matrix M.

   Every directed arc M_ij=1 gives a double-intersection pair {i,j}, because
   row j is isolated in column j. There are30 arcs and15 double pairs, with
   at most two opposite arcs per pair, so every double pair has both arcs.
   Hence M is symmetric. Adjacent vertices of M have no common neighbor,
   and nonadjacent vertices have exactly one: the support is I plus the
   Petersen parameters(10,3,0,1). This parameter argument does not require
   appealing to uniqueness of the Petersen graph.

   Any three vertices have at least eight vertices in the union of their
   closed neighborhoods. Indeed their induced edge count e is0,1,or2; pair
   intersections sum to3+e, so the union size is9-e+t, where t is the triple
   intersection. If e=2, the path center gives t>=1. Thus X has no3x3 zero
   submatrix. But all C events have type(3,3) by the preceding degree and
   anticommutator argument, and their A/B neighborhoods require exactly such
   a zero submatrix by triangle-freeness. Contradiction.

5. ALL STARS ARE IMPOSSIBLE EVEN FOR THE PAIR. Each row's six singleton
   intersections force it to be a star center once: if it is a center h times,
   the count is3*h+(4-h)=6. Relabel these centers onto the diagonal. Then

       X=(2*I+M)/5,
       M M*+2*(M+M*)=I+2*J,
       M* M+2*(M+M*)=I+2*J,

   where M is a zero-diagonal0/1 matrix with every row and column sum three.
   The off-diagonal identities forbid opposite arcs, and give zero common
   out-neighbors for adjacent pairs and two for nonadjacent pairs. The graph
   D with adjacency J-I-M-M* is therefore cubic on ten vertices. Each of the
   ten in-neighborhood triples of M is a triangle of D. They are distinct:
   equal triples would give a column intersection of three, while the second
   identity allows at most two. A cubic graph on ten vertices has at most
   ten triangles, since each vertex lies in at most three. Equality here
   forces every vertex's three neighbors to form a triangle, so every
   connected component is K4. Ten is not divisible by four. Contradiction.

These cases exhaust a uniform-probability isoclinic pair. Therefore a uniform
dimension19 three-context GHZ realization cannot have even one reduced
isoclinic pair. Its three pairs must all have unequal principal angles.
The general nonuniform case remains open.
"""
from fractions import Fraction as F
import itertools
import json


def petersen_pair_control():
    """An exact complex isoclinic PAIR which admits no required third context.

    X=((7+i)A+(4-3i)I)/25, where A joins disjoint two-subsets of five.
    Keep real/imaginary parts as Fractions so no complex floating point is used.
    """
    vertices = list(itertools.combinations(range(5), 2))
    A = [[int(set(u).isdisjoint(v)) for v in vertices] for u in vertices]
    assert all(sum(row) == 3 for row in A)
    R = [[F(7*A[i][j]+4*int(i == j), 25) for j in range(10)] for i in range(10)]
    I = [[F(A[i][j]-3*int(i == j), 25) for j in range(10)] for i in range(10)]
    assert all(sum(row) == 1 for row in R)
    assert all(sum(row) == 0 for row in I)
    for i in range(10):
        for j in range(10):
            real = sum(R[i][k]*R[j][k]+I[i][k]*I[j][k] for k in range(10))
            imag = sum(I[i][k]*R[j][k]-R[i][k]*I[j][k] for k in range(10))
            assert real == F(int(i == j), 5)+F(2, 25) and imag == 0
    supports = [{j for j in range(10) if R[i][j] or I[i][j]} for i in range(10)]
    assert all(len(row) == 4 for row in supports)
    unions = [len(set().union(*(supports[i] for i in triple)))
              for triple in itertools.combinations(range(10), 3)]
    assert min(unions) == 8
    real_pair = [[F(2*A[i][j]+int(i == j), 7) for j in range(10)] for i in range(10)]
    assert all(sum(row) == 1 for row in real_pair)
    for i in range(10):
        for j in range(10):
            assert sum(real_pair[i][k]*real_pair[j][k] for k in range(10)) == F(9*int(i == j)+4, 49)
    return {'squared_cosine': '1/5', 'offdiagonal_gram': '2/25',
            'exact_complex_pair_verified': True,
            'exact_real_pair_squared_cosine': '9/49',
            'row_triples_checked': len(unions), 'minimum_closed_neighborhood_union': min(unions)}


def verify():
    if not __debug__:
        raise RuntimeError('Exact checks require Python assertions; do not use -O')
    # 9*u*(1-u)>=2 iff (3*u-1)*(3*u-2)<=0.
    for u in (F(i, 120) for i in range(121)):
        assert 9*u*(1-u)-2 == -(3*u-1)*(3*u-2)
        assert (9*u*(1-u) >= 2) == (F(1, 3) <= u <= F(2, 3))
    for u in (F(1, 3), F(2, 3)):
        assert 9*(u*u-u)/2 == -1
        # Triple product of a,b,(a-b)/sqrt(2-2u).
        assert -u*(1-u)/2 == F(-1, 9)
    degrees = []
    saturated = []
    for k in range(1, 10):
        for l in range(1, 10-k):
            reciprocal = F(1, k)+F(1, l)+F(1, 10-k-l)
            if reciprocal <= 1:
                degrees.append((k, l))
                if reciprocal == 1:
                    saturated.append((k, l))
    assert degrees == [(2, 4), (3, 3), (3, 4), (4, 2), (4, 3), (4, 4)]
    assert saturated == [(2, 4), (4, 2), (4, 4)]
    assert (1+9*F(4, 9))/10 == F(1, 2)
    assert (1+9*F(1, 6))/10 == F(1, 4)
    impossible_common_neighbor_count = 16*(1-F(1, 6))/10
    assert impossible_common_neighbor_count == F(4, 3)
    assert impossible_common_neighbor_count.denominator != 1
    surviving_types = [pair for pair in degrees if pair[0] >= 3 and pair not in [(4, 2), (4, 4)]]
    assert surviving_types == [(3, 3), (3, 4), (4, 3)]
    capacities = {r: (10-r)*3+r*6 for r in range(11)}
    assert all((capacity >= 45) == (r >= 5) for r, capacity in capacities.items())
    compatible_counts = [(r, s) for r in range(11) for s in range(11)
                         if capacities[r] >= 45 and capacities[s] >= 45 and r+s <= 10]
    assert compatible_counts == [(5, 5)]
    assert 5*3+5*4 == 35
    assert 10-(3+4-1) == 4  # Outside points for the elementary incidence proof.
    for k, l, forced_x in ((2, 4, F(1, 6)), (4, 2, F(1, 6)), (4, 4, F(4, 9))):
        # Saturation gives ||E_A g||^2=1/k and ||E_B g||^2=1/l.
        assert F(1, k)+F(1, l) == (1-forced_x)*F(9, 10)
    assert (1+9*F(7, 27))/10 == F(1, 3)
    # Local singleton-intersection graph exclusions, valid at every k>0.
    assert 3*F(1, 3)**2 != 1-9*F(1, 9)
    assert 4*F(1, 4)**2 != 1-9*F(1, 16)
    balanced_forced_k = F(1, 10)
    assert F(1, 2)-4*balanced_forced_k == 1-9*balanced_forced_k
    assert balanced_forced_k > F(1, 16)
    assert [a for a in range(11) if 30+3*a <= 30-a] == [0]
    assert 30+15 == 45 and 30+2*15 == 60
    assert [h for h in range(5) if 2*(4-h) == 6] == [1]
    assert [h for h in range(5) if 3*h+(4-h) == 6] == [1]
    star_k = F(2, 25)
    center, leaf = F(2, 5), F(1, 5)
    assert center+3*leaf == 1 and center*leaf == star_k
    assert center*center+3*leaf*leaf == 1-9*star_k
    triangle_real, triangle_imag = F(7, 25), F(1, 25)
    assert triangle_real*triangle_real+triangle_imag*triangle_imag == star_k
    assert triangle_real == 7*star_k/2
    assert 3*star_k+F(4, 25)**2+F(3, 25)**2 == 1-9*star_k
    union_lower_bounds = [9, 8, 9-2+1]
    assert min(union_lower_bounds) == 8
    assert 10 % 4 != 0
    return {'general_isoclinic_cosine_interval': '(1/3,2/3)',
            'uniform_isoclinic_squared_cosine_interval': '(1/6,4/9)',
            'uniform_degree_types': degrees,
            'saturated_types': saturated,
            'isoclinic_endpoint_degree_types': surviving_types,
            'forced_degree_four_counts': compatible_counts,
            'impossible_common_neighbor_count_at_squared_cosine_1_over_6': str(impossible_common_neighbor_count),
            'one_pair_uniform_degree_types_in_every_context': surviving_types,
            'singleton_intersection_count_forces_four_regular_pair': True,
            'four_regular_support_cases': ['all stars impossible', 'all triangles require forbidden3x3zero block'],
            'petersen_pair_positive_control': petersen_pair_control(),
            'scope': 'even one isoclinic pair excluded with all30 probabilities equal1/10',
            'allowed_supergraph_warning': 'disjoint-neighborhood necessity concerns actual nonzero overlaps'}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))

"""Exact checks for regular pair-block obstructions on thirty events.

Suppose the bipartite supports A--B and A--C are q- and r-regular. Set
F(alpha,beta,rho)=1/(1/alpha+1/beta+1/rho). This harmonic mean is strictly
concave on the positive probability simplex. Since the neighbor masses sum
to q and r over A, Lemma14 and Jensen imply

    1 <= 10 F(q/10,r/10,1-(q+r)/10).

The only ordered degree pairs are (2,4),(3,3),(3,4),(4,2),(4,3),(4,4).
For (2,4) or (4,4), equality holds. Thus every A event has probability1/10,
and all sector masses equal (q/10,r/10,1-(q+r)/10). After rephasing all
state overlaps positive, equality in the three Cauchy--Schwarz bounds gives

    <g_a,g_b> = sqrt(10)*sqrt(c_b)/q  on A--B edges,

and a weighted sum L=a P_B+b P_C acts as identity on S_A orthogonal to psi.
With rho=1-(q+r)/10, the coefficients are a=1+(q/10)/rho,
b=1+(r/10)/rho. Thus a=b=3 for (4,4), and a=3/2,b=2 for (2,4).

For arbitrary projections E,F and L=aE+bF, the exact identity

    a(LE+EL) = L^2-bL+a(a+b)E

implies that compression of E to any lambda-eigenspace of L is the scalar
lambda(lambda-b)/(a(2lambda-a-b)), when the denominator is nonzero. Taking
lambda1 on S_A orthogonal to psi gives compression1/6 for (4,4), and4/9
for (2,4). Because P_B psi=psi, the common-B probability mass for distinct
A vertices is therefore (1-t)*q^2/100, where t is that compression scalar.

For (4,4) this is2/15>1/9, so all45 A-pairs need at least2 common neighbors;
but regularity supplies only10*C(4,2)=60 incidences. For (2,4) this is1/45>0,
so all45 A-pairs need a common neighbor, but only10*C(2,2)=10 are available.
Both contradictions are dimension-independent and allow complex overlaps
and nonuniform probabilities outside A. Hence if all three pair blocks are
regular, their degree triple must be (3,3,3) or a permutation of (3,3,4).

Optional orbit observation: a vertex with one neighbor in another context
has c<h(1/9)=2/27, using h(t)=t(1-t)/(1+3t), increasing on (0,1/3).
An automorphism orbit of size m>=3 cannot have cross-degree1: its total mass
would be below2m/27, but is greater than(m-1)/9 for m<10 (equal1 for m10).

Future C3 x C3 construction capsule (no implementation here): faithful
actions on10 have orbit types R=9+1, T=3+3+3+1, or D=3+3+1+1+1+1. At least
two distinct kernel lines are required among the3-orbits in T,D. R cannot
pair with D: the four fixed D vertices have only two possible nonempty
proper invariant neighborhoods in R, contradicting distinct pair rows.
For two T/D contexts the sets of3-orbit kernel lines must coincide: between
different kernel lines K,L the action on (G/K)x(G/L) is transitive, so its
block is all-or-none and cannot distinguish points in either orbit. Thus
RTT only needs shared sets of2 or3 kernel lines. A regular9 to3-orbit block
is specified by a nonempty proper subset of F3 (size1 or2). Fourier blocks
are small and may help unrestricted theta searches, but averaging a Gram
matrix need not preserve rank19. The all-R action is already excluded by
the large-orbit theorem. These reductions are necessary, not exclusions of
the remaining mixed families.
"""
from fractions import Fraction as Q
from math import comb
import json


def add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for word, coefficient in polynomial.items():
            result[word] = result.get(word, Q(0)) + coefficient
    return {word: coefficient for word, coefficient in result.items() if coefficient}


def scale(polynomial, scalar):
    return {word: scalar * coefficient for word, coefficient in polynomial.items()}


def multiply(left, right):
    result = {}
    for u, a in left.items():
        for v, b in right.items():
            word = u + v
            # The only relations used are E^2=E and F^2=F.
            word = ''.join(letter for i, letter in enumerate(word)
                           if i == 0 or word[i-1] != letter)
            result[word] = result.get(word, Q(0)) + a * b
    return result


def check_case(q, r):
    alpha, beta, rho = Q(q, 10), Q(r, 10), Q(10 - q - r, 10)
    assert 10 / (1 / alpha + 1 / beta + 1 / rho) == 1
    a, b = 1 + alpha / rho, 1 + beta / rho
    E, F = {'E': Q(1)}, {'F': Q(1)}
    L = add(scale(E, a), scale(F, b))
    lhs = scale(add(multiply(L, E), multiply(E, L)), a)
    rhs = add(multiply(L, L), scale(L, -b), scale(E, a * (a + b)))
    assert lhs == rhs
    compression = (1 - b) / (a * (2 - a - b))
    common_mass = (1 - compression) * Q(q * q, 100)
    required_common_neighbors = 2 if common_mass > Q(1, 9) else 1
    required_incidence = comb(10, 2) * required_common_neighbors
    available_incidence = 10 * comb(q, 2)
    assert common_mass > 0 and available_incidence < required_incidence
    return {'degrees': [q, r], 'projection_coefficients': [str(a), str(b)],
            'projection_compression': str(compression),
            'common_neighbor_mass': str(common_mass),
            'required_pair_incidences': required_incidence,
            'available_pair_incidences': available_incidence}


def main():
    allowed = [(q, r) for q in range(1, 10) for r in range(1, 10)
               if q + r < 10 and Q(1, q) + Q(1, r) + Q(1, 10-q-r) <= 1]
    assert allowed == [(2, 4), (3, 3), (3, 4), (4, 2), (4, 3), (4, 4)]
    triples = [(q, r, s) for q in range(1, 10) for r in range(q, 10)
               for s in range(r, 10)
               if all(pair in allowed and pair not in ((2, 4), (4, 2), (4, 4))
                      for pair in ((q, r), (q, s), (r, s)))]
    assert triples == [(3, 3, 3), (3, 3, 4)]
    print(json.dumps({'regular_block_exclusions': [check_case(4, 4), check_case(2, 4)],
                      'remaining_sorted_degree_triples': triples}, indent=2))


if __name__ == '__main__':
    main()

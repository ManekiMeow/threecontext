"""Every transitive permutation group of degree9 contains regular C9 or C3^2.

Proof. Let H stabilize a point in the transitive group G. Choose a Sylow3
subgroup Q of H and extend it to a Sylow3 subgroup P of G. Since [G:H]=9,
P intersect H has the same order as Q and [P:P intersect H]=9. Thus P is
transitive. Conjugate P into a Sylow3 subgroup W=C3 wr C3 of S9, whose base
V=F3^3 rotates within three blocks and whose top C3 cycles the blocks.

P projects onto the top C3, and B=P intersect V is nonzero, since |P|>=9.
Choose g projecting to a top generator. Conjugation by g cyclically permutes
the three base coordinates, so B is stable under the cyclic operator T.
Since (T-I)^3=0, applying the largest nonzero power of T-I to any nonzero
base vector gives a nonzero T-fixed vector. The fixed space is precisely
the diagonal line, hence B contains z=(1,1,1). This z is central in W.
If g has order9, it acts as a9-cycle. Otherwise g has order3. Then g and z
commute and generate C3^2: g moves between blocks, whereas z cycles every
block internally. Their action is transitive on9 points and hence regular.

For a context-preserving automorphism group acting as a transitive9-orbit
plus a fixed event, Proposition39 makes the action on that context faithful:
an element fixing that entire context fixes every vertex, by injectivity of
cross-context neighborhoods. Thus the group embeds in S9 and the preceding
argument applies. The quantum nonexistence corollary then follows from the
separate order9 and regular-C3^2 action exclusions; this checker only proves
and enumerates the finite permutation-group step.

Independent exact finite check. All28 linear subspaces of F3^3 are enumerated;
exactly four are T-invariant, with dimensions0,1,2,3. Every transitive subgroup
of W has nonzero invariant B and a lift g=(v,1), so it is B union Bg union Bg^2.
Enumerating all such B and v yields exactly13 transitive subgroups of W.
For each, check closure, transitivity, and an explicit regular C9 or C3^2.
No group-theory package or numerical arithmetic is used.
"""
from itertools import product
import json


ZERO = (0, 0, 0)
VECTORS = tuple(product(range(3), repeat=3))
IDENTITY = (ZERO, 0)


def add(v, w):
    return tuple((x+y) % 3 for x, y in zip(v, w))


def scalar(k, v):
    return tuple(k*x % 3 for x in v)


def shift(v, k=1):
    return tuple(v[(i-k) % 3] for i in range(3))


def multiply(g, h):
    v, t = g
    w, u = h
    return add(v, shift(w, t)), (t+u) % 3


def order(g):
    value = IDENTITY
    for n in range(1, 10):
        value = multiply(value, g)
        if value == IDENTITY:
            return n
    raise AssertionError('An element of this wreath product must have order at most9')


def image(g, point):
    v, t = g
    block, position = divmod(point, 3)
    new_block = (block+t) % 3
    return 3*new_block+(position+v[new_block]) % 3


def subspaces():
    spaces = {frozenset({ZERO})}
    todo = list(spaces)
    while todo:
        space = todo.pop()
        for v in VECTORS:
            if v in space:
                continue
            enlarged = frozenset(add(w, scalar(k, v)) for w in space for k in range(3))
            if enlarged not in spaces:
                spaces.add(enlarged)
                todo.append(enlarged)
    return spaces


def main():
    elements = [(v, t) for v in VECTORS for t in range(3)]
    assert all(image(multiply(g, h), point) == image(g, image(h, point))
               for g, h, point in product(elements, elements, range(9)))
    N = lambda v: add(shift(v), scalar(2, v))
    assert all(N(N(N(v))) == ZERO for v in VECTORS)
    assert {v for v in VECTORS if N(v) == ZERO} == {(k, k, k) for k in range(3)}
    spaces = subspaces()
    invariant = [B for B in spaces if {shift(v) for v in B} == set(B)]
    assert len(spaces) == 28 and sorted(map(len, invariant)) == [1, 3, 9, 27]
    z = ((1, 1, 1), 0)
    assert all(multiply(g, z) == multiply(z, g) for g in elements)
    subgroups = set()
    for B in invariant:
        if len(B) == 1:
            continue
        assert z[0] in B
        for v in VECTORS:
            g = (v, 1)
            g2 = multiply(g, g)
            assert multiply(g2, g)[0] in B
            subgroup = frozenset(multiply((b, 0), power)
                                 for b in B for power in (IDENTITY, g, g2))
            assert len(subgroup) == 3*len(B)
            subgroups.add(subgroup)
    assert len(subgroups) == 13
    counts = {'contains_regular_C9': 0, 'regular_C3square_without_C9': 0}
    sizes = {}
    for P in subgroups:
        assert {image(g, 0) for g in P} == set(range(9))
        assert all(multiply(g, h) in P for g, h in product(P, repeat=2))
        sizes[len(P)] = sizes.get(len(P), 0)+1
        order9 = next((g for g in P if order(g) == 9), None)
        if order9 is not None:
            power, orbit = IDENTITY, set()
            for _ in range(9):
                orbit.add(image(power, 0))
                power = multiply(power, order9)
            assert len(orbit) == 9 and power == IDENTITY
            counts['contains_regular_C9'] += 1
        else:
            g = next(g for g in P if g[1] == 1)
            assert order(g) == order(z) == 3 and multiply(g, z) == multiply(z, g)
            g2, z2 = multiply(g, g), multiply(z, z)
            regular = {multiply(u, v) for u in (IDENTITY, g, g2)
                       for v in (IDENTITY, z, z2)}
            assert len(regular) == 9 and regular <= P
            assert {image(h, 0) for h in regular} == set(range(9))
            counts['regular_C3square_without_C9'] += 1
    assert sizes == {9: 9, 27: 3, 81: 1}
    assert counts == {'contains_regular_C9': 9, 'regular_C3square_without_C9': 4}
    print(json.dumps({'linear_subspaces': len(spaces), 'invariant_subspaces': len(invariant),
                      'transitive_subgroups_of_Sylow81': len(subgroups),
                      'subgroup_order_counts': sizes, 'verified_regular_subgroups': counts}, indent=2))


if __name__ == '__main__':
    main()

"""Exact representatives for C3 x C3 RTT with opposed kernel multiplicities.

A has a regular9-orbit plus fixed a. B has3-orbits X0,X1,Y and fixed b;
C has3-orbits X,Y0,Y1 and fixed c, with distinct kernel lines X,Y.
The reduction uses Proposition39 neighborhood injectivity, the strict
probability bound, and balanced pair components. It needs no Gram averaging.

Different-kernel BC blocks are empty, because their nonempty moving-A coset
neighborhoods intersect. A-to3-orbit supports have size1 or2, since empty or
full supports repeat rows. No transitive3-orbit has cross-degree1: the
single-neighbor bound gives every probability<2/27 and orbit mass<2/9,
whereas its seven-event complement has mass<7/9.

If b sees moving A, it cannot touch moving C, so it sees only c in C. The
kernel having two C3-orbits but only one B3-orbit would then force all six
C vertices to have degree2 into the same B3-set: degree1 is excluded and
degree3 repeats rows. Only three2-subsets exist, violating row injectivity.
The symmetric mixed case is identical. Thus ab=ac=1,bc=0.

As in the three-kernel RTT reduction, a touches exactly one B3-orbit and
one C3-orbit. Same-kernel choices are impossible: one selected orbit would
have no neighbor in the opposite context, since its only same-kernel block
and opposite fixed-point edge are both forbidden by triangles. Selecting
the two singleton-kernel orbits B_Y,C_X is also impossible: BC splits into
an X-side containing6 B-events and4 C-events (including c), and a Y-side
containing4 B-events and6 C-events (including b). No union of balanced pair
components can have unequal side sizes. Hence select B_X0 and C_Y0.

B_X0 has C-neighbors only in C_X, and C_Y0 has B-neighbors only in B_Y.
Their cyclic supports must therefore have size2. Triangle-freeness forces
the corresponding four A-supports to have size1. For each remaining orbit
B_X1,C_Y1, let u be its A-support size and s its same-kernel BC support size.
The only possibilities are(u,s)=(1,2),(1,1),(2,1). Its opposite fixed-point
edge is forced: without it s1 gives degree1, while s2 repeats the three
neighborhoods of the selected orbit. Thus c--B_X1 and b--C_Y1 are present.
At least one of b--C_X,c--B_Y must be present; otherwise the same6:4 and4:6
partition contradicts balanced pair components.

Choose quotient coordinates X=x,Y=y. Independently translate fiber labels
to normalize singleton A-supports to{0} and double supports to{1,2}.
The size2 BC supports are uniquely{1,2}; a size1 support is either{1} or{2}
when both A-supports are singleton, or uniquely{0} when one is double.
The two remaining singleton phases are removed independently by x->-x
and y->-y, acting also on all fibers of the corresponding kernel.
There are therefore at most3*3*3=27 canonical representatives. The constructor and explicit
global phase isomorphisms below check this classification's final arithmetic.

Run with --generate to discover rational duals (requires the project venv);
run without arguments to replay all saved duals exactly using stdlib only.
Numerical failure would be inconclusive; only verified bounds<3 exclude.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path


MODES = ('cycle', 'matching', 'wide')
FLAGS = ((0, 1), (1, 0), (1, 1))
FOLDER = Path(__file__).with_name('event30_c3square_opposed_certificates')


def canonical_graph(mode_B, mode_C, flag_bX, flag_cY, phase_B=1, phase_C=1):
    assert mode_B in MODES and mode_C in MODES
    assert (flag_bX, flag_cY) in FLAGS and phase_B in (1, 2) and phase_C in (1, 2)
    edges = set()
    def add(u, v):
        edges.add(tuple(sorted((u, v))))
    for x, y in product(range(3), repeat=2):
        point = 3*x+y
        for k, label in enumerate((x, x, y)):
            for j in range(3):
                if (j != label if mode_B == 'wide' and k == 1 else j == label):
                    add(point, 10+3*k+j)
        for k, label in enumerate((x, y, y)):
            for j in range(3):
                if (j != label if mode_C == 'wide' and k == 2 else j == label):
                    add(point, 20+3*k+j)
    for i, j in product(range(3), repeat=2):
        if i != j:
            add(10+i, 20+j)
            add(16+i, 23+j)
        condition_B = (i != j) if mode_B == 'cycle' else ((j-i) % 3 == phase_B if mode_B == 'matching' else i == j)
        condition_C = (i != j) if mode_C == 'cycle' else ((j-i) % 3 == phase_C if mode_C == 'matching' else i == j)
        if condition_B:
            add(13+i, 20+j)
        if condition_C:
            add(16+i, 26+j)
    add(9, 19)
    add(9, 29)
    for j in range(3):
        for u, v in ((9, 10+j), (9, 23+j), (29, 13+j), (19, 26+j)):
            add(u, v)
        if flag_bX:
            add(19, 20+j)
        if flag_cY:
            add(29, 16+j)
    return {'n': 30, 'parts': [list(range(10*k, 10*k+10)) for k in range(3)],
            'edges': [list(edge) for edge in sorted(edges)]}


def edge_image(graph, permutation):
    return {tuple(sorted((permutation[u], permutation[v]))) for u, v in graph['edges']}


def phase_checks():
    for axis in (0, 1):
        permutation = list(range(30))
        for x, y in product(range(3), repeat=2):
            permutation[3*x+y] = 3*((-x if axis == 0 else x) % 3)+(-y if axis == 1 else y) % 3
        starts = (10, 13, 20) if axis == 0 else (16, 23, 26)
        for start in starts:
            for j in range(3):
                permutation[start+j] = start+(-j) % 3
        assert sorted(permutation) == list(range(30))
        for mode_B, mode_C, flags in product(MODES, MODES, FLAGS):
            before = canonical_graph(mode_B, mode_C, *flags)
            after = canonical_graph(mode_B, mode_C, *flags,
                                    phase_B=2 if axis == 0 else 1,
                                    phase_C=2 if axis == 1 else 1)
            assert edge_image(before, permutation) == {tuple(edge) for edge in after['edges']}
    return 54


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    args = parser.parse_args()
    from overnight30_screen import base_validation, verify_certificate
    if args.generate:
        from overnight30_screen import exact_dual_certificate
        FOLDER.mkdir(exist_ok=True)
    results = []
    for mode_B, mode_C, flags in product(MODES, MODES, FLAGS):
        graph = canonical_graph(mode_B, mode_C, *flags)
        base_validation(graph)
        name = f'{mode_B}_{mode_C}_{flags[0]}{flags[1]}'
        path = FOLDER/(name+'.json')
        if args.generate:
            certificate, optimum = exact_dual_certificate(graph)
            if certificate is None:
                raise RuntimeError(f'No strict exact dual found for {name}; numeric={optimum}')
            path.write_text(json.dumps(certificate, indent=2)+'\n')
        else:
            certificate = json.loads(path.read_text())
        assert certificate['graph'] == graph and certificate['method'] == 'theta_dual'
        bound = F(certificate['witness']['bound'])
        assert bound < 3
        verify_certificate(certificate)
        results.append({'name': name, 'edges': len(graph['edges']), 'exact_bound': str(bound)})
    print(json.dumps({'verified_canonical_supports': len(results),
                      'verified_independent_phase_isomorphisms': phase_checks(),
                      'largest_exact_bound': str(max(F(row['exact_bound']) for row in results)),
                      'certificates': results}, indent=2))


if __name__ == '__main__':
    main()

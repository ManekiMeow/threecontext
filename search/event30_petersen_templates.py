"""Finite maximal-neighborhood reduction for a closed-Petersen context pair.

Fix A--B support I+P, where P=KG(5,2). A third-context event has two nonempty
neighbor sets S,T with no edge between them. Enlarging S,T to a maximal zero
rectangle preserves triangle-freeness and relaxes overlap constraints. The
110 maximal rectangles have types (1,6),(6,1):10 each; (2,4),(4,2):15 each;
and (2,3),(3,2):30 each. Hence every putative ten-event C context extends
to ten such templates. Repeated templates give same-context false twins and
are impossible by contraction, so it suffices to select ten distinct templates.

This is a finite reduction, NOT an exclusion of the family. The bounded SAT
pilot enforces only exact necessary support conditions and applies previously
verified graph obstructions. Surviving graphs are not quantum realizations.

Important scalar correction: a (2,4) event has upper bound 14/135, not 4/39.
The relaxed masses (2/9,7/18,7/18) give the former. Thus a single (1,6)
template is not excluded by the proposed budget: 9*(14/135)+2/27=136/135>1.
Two singleton templates ARE excluded: 8*(14/135)+2*(2/27)=44/45<1.
The pilot may therefore select at most one singleton template.
"""
if not __debug__:
    raise RuntimeError('Exact template verification requires assertions; do not use -O.')

from collections import Counter
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import time

PARTS = [list(range(10*k, 10*k+10)) for k in range(3)]


def petersen():
    vertices = list(itertools.combinations(range(5), 2))
    neighbors = [{j for j, b in enumerate(vertices) if set(a).isdisjoint(b)} | {i}
                 for i, a in enumerate(vertices)]
    return vertices, neighbors


def templates():
    _, neighbors = petersen()
    result = []
    for mask in range(1, 1024):
        left = {i for i in range(10) if mask >> i & 1}
        right = set(range(10))-set.union(*(neighbors[i] for i in left))
        if not right:
            continue
        closure = set(range(10))-set.union(*(neighbors[j] for j in right))
        if closure == left:
            result.append((tuple(sorted(left)), tuple(sorted(right))))
    return result


def action():
    """240 explicitly verified template permutations: S5 and interchange A,B."""
    vertices, neighbors = petersen()
    choices = templates()
    index = {item: i for i, item in enumerate(choices)}
    result = []
    for sigma in itertools.permutations(range(5)):
        permutation = [vertices.index(tuple(sorted((sigma[a], sigma[b])))) for a, b in vertices]
        assert all({permutation[j] for j in neighbors[i]} == neighbors[permutation[i]]
                   for i in range(10))
        for swap in (False, True):
            transformed = []
            for left, right in choices:
                pair = (tuple(sorted(permutation[v] for v in left)),
                        tuple(sorted(permutation[v] for v in right)))
                transformed.append(index[pair[::-1] if swap else pair])
            assert sorted(transformed) == list(range(110))
            result.append(transformed)
    assert len({tuple(p) for p in result}) == 240
    return result


def graph(selected):
    from overnight30_screen import normalize_graph, base_validation
    assert selected == sorted(set(selected)) and len(selected) == 10
    choices = templates()
    assert all(type(i) is int and 0 <= i < len(choices) for i in selected)
    _, neighbors = petersen()
    edges = [[i, 10+j] for i, row in enumerate(neighbors) for j in sorted(row)]
    for c, index in enumerate(selected, 20):
        left, right = choices[index]
        edges.extend([a, c] for a in left)
        edges.extend([10+b, c] for b in right)
    record = normalize_graph({'parts': PARTS, 'edges': edges})
    base_validation(record)
    return record


def support_constraints(selected):
    """Independent combinatorial verification of the SAT model constraints."""
    from overnight30_screen import base_validation
    record = graph(selected)
    adj, _ = base_validation(record)
    assert all(1 <= len(adj[v] & set(PARTS[2])) <= 6 for v in range(20))
    for part in PARTS:
        for u, v in itertools.combinations(part, 2):
            assert any(adj[a] & adj[v] for a in adj[u])
    assert sum(min(map(len, templates()[i])) == 1 for i in selected) <= 1
    return record


def build(matchings=False):
    from pysat.card import CardEnc, EncType
    from pysat.formula import IDPool
    choices = templates()
    _, neighbors = petersen()
    pool = IDPool(start_from=111)
    clauses = CardEnc.equals(list(range(1, 111)), bound=10, vpool=pool,
                             encoding=EncType.seqcounter).clauses
    singleton = [i+1 for i, (s, t) in enumerate(choices) if min(len(s), len(t)) == 1]
    clauses += CardEnc.atmost(singleton, bound=1, vpool=pool, encoding=EncType.seqcounter).clauses
    for side in (0, 1):
        for vertex in range(10):
            incident = [i+1 for i, pair in enumerate(choices) if vertex in pair[side]]
            clauses.append(incident)
            clauses += CardEnc.atmost(incident, bound=6, vpool=pool,
                                     encoding=EncType.seqcounter).clauses
        for u, v in itertools.combinations(range(10), 2):
            eligible = [i+1 for i, pair in enumerate(choices)
                        if (u in pair[side] and neighbors[v].intersection(pair[1-side]))
                        or (v in pair[side] and neighbors[u].intersection(pair[1-side]))]
            clauses.append(eligible)
    incompatible = []
    for i, j in itertools.combinations(range(110), 2):
        s, t = choices[i]; u, v = choices[j]
        if not any(neighbors[a].intersection(v) for a in s) and not any(neighbors[a].intersection(t) for a in u):
            incompatible.append([i, j]); clauses.append([-i-1, -j-1])
    matching_variables = []
    if matchings:
        # Ten matched vertices and ten selected templates, with injective
        # matching and match=>selection, cover every selected template too.
        for side in (0, 1):
            variables = {(v, i): pool.id() for v in range(10)
                         for i, pair in enumerate(choices) if v in pair[side]}
            for v in range(10):
                incident = [var for (vertex, i), var in variables.items() if vertex == v]
                clauses += CardEnc.equals(incident, bound=1, vpool=pool,
                                          encoding=EncType.seqcounter).clauses
            for i in range(110):
                incident = [var for (v, template), var in variables.items() if template == i]
                clauses += CardEnc.atmost(incident, bound=1, vpool=pool,
                                          encoding=EncType.seqcounter).clauses
                clauses.extend([-var, i+1] for var in incident)
            matching_variables.append(variables)
    return clauses, {'variables': pool.top, 'clauses': len(clauses),
                     'templates': len(choices), 'incompatible_template_pairs': len(incompatible),
                     'perfect_matching_constraints': matchings}, matching_variables


def self_test():
    choices = templates()
    assert len(choices) == 110
    counts = Counter(tuple(map(len, item)) for item in choices)
    assert counts == {(1, 6): 10, (6, 1): 10, (2, 4): 15, (4, 2): 15,
                      (2, 3): 30, (3, 2): 30}
    _, neighbors = petersen()
    for mask in range(1, 1024):
        left = {i for i in range(10) if mask >> i & 1}
        right = set(range(10))-set.union(*(neighbors[i] for i in left))
        if right:
            closure = set(range(10))-set.union(*(neighbors[j] for j in right))
            assert (tuple(sorted(closure)), tuple(sorted(right))) in choices
            assert left <= closure
    assert F(1, 1)/(F(9, 2)+F(18, 7)+F(18, 7)) == F(14, 135)
    assert 9*F(14, 135)+F(2, 27) == F(136, 135) > 1
    assert 8*F(14, 135)+2*F(2, 27) == F(44, 45) < 1
    return {'maximal_templates': 110, 'type_counts': {str(k): v for k, v in counts.items()},
            'explicit_template_automorphisms': len(action()),
            'single_singleton_budget': '136/135', 'two_singleton_budget': '44/45'}


def pilot(output, seconds, matchings=False):
    from pysat.solvers import Cadical153
    from overnight30_screen import fast_screen, graph_id, verify_certificate, unbalanced_coloring, certificate
    output.mkdir(parents=True, exist_ok=False)
    controls = self_test()
    start = time.monotonic()
    clauses, metadata, matching_variables = build(matchings)
    permutations = action()
    count, methods, surviving, blocked, unknown = 0, Counter(), 0, set(), False
    with Cadical153(bootstrap_with=clauses) as solver, (output/'records.jsonl').open('w') as stream:
        while time.monotonic()-start < seconds:
            solver.conf_budget(2000)
            solved = solver.solve_limited(expect_interrupt=True)
            if solved is None:
                unknown = True
                break
            if not solved:
                # Exhaustion would require a separately published exact proof.
                metadata['solver_reports_exhaustion'] = True
                break
            full_model = {i for i in solver.get_model() if i > 0}
            selected = sorted(i-1 for i in full_model if 1 <= i <= 110)
            record = support_constraints(selected)
            matching_witness = []
            for side, variables in enumerate(matching_variables):
                matched = [[v, i] for (v, i), var in variables.items() if var in full_model]
                assert sorted(v for v, i in matched) == list(range(10))
                assert sorted(i for v, i in matched) == selected
                assert all(v in templates()[i][side] for v, i in matched)
                matching_witness.append(matched)
            orbit = {tuple(sorted(p[i] for i in selected)) for p in permutations}
            assert not orbit.intersection(blocked)
            blocked.update(orbit)
            for equivalent in orbit:
                solver.add_clause([-i-1 for i in equivalent])
            screened = fast_screen(record)
            if screened['status'] != 'excluded' and time.monotonic()-start < seconds-0.2:
                coloring, status = unbalanced_coloring(record, conflict_budget=1000)
                screened['coloring_search'] = status
                if coloring is not None:
                    screened.update(status='excluded', certificate=certificate(
                        record, 'unbalanced_coloring', coloring=coloring, discovery='SAT witness'))
            proof = screened.get('certificate')
            if proof is not None:
                verify_certificate(proof)
                methods[proof['method']] += 1
            else:
                surviving += 1
            stream.write(json.dumps({'selected': selected, 'graph_id': graph_id(record),
                                     'graph': record, 'status': screened['status'],
                                     'certificate': proof,
                                     'matchings': matching_witness,
                                     'coloring_search': screened.get('coloring_search')}, separators=(',', ':'))+'\n')
            stream.flush()
            count += 1
    summary = dict(metadata, controls=controls, seconds_requested=seconds,
                   seconds_elapsed=time.monotonic()-start, representatives=count,
                   selected_subsets_blocked=len(blocked), survivors=surviving,
                   exclusion_methods=dict(methods), limited_sat_unknown=unknown,
                   numerical_solves=0,
                   scope='Bounded representatives under the explicit240-element action, not full graph-isomorphism classes or a family exclusion.')
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot', type=Path)
    parser.add_argument('--seconds', type=float, default=18.0)
    parser.add_argument('--matchings', action='store_true')
    args = parser.parse_args()
    print(json.dumps(pilot(args.pilot, args.seconds, args.matchings) if args.pilot else self_test(), indent=2))

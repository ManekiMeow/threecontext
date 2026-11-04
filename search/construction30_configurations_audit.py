"""Independent audit of the linear 10_3 configuration certificates.

Uses a different cubic-graph traversal without WL hashing, a SAT exact-cover
encoding, and integer Bareiss elimination instead of rational Gaussian
elimination. Does not import construction30_configurations.py or use an SDP.

Requires NetworkX and python-sat:
    python search/construction30_configurations_audit.py
"""
import itertools
import json
from pathlib import Path

import networkx as nx
from pysat.solvers import Cadical153


def cubic_graphs():
    seen, final = {}, []

    def visit(graph):
        key = tuple(sorted(dict(graph.degree()).values()))
        bucket = seen.setdefault(key, [])
        if any(nx.is_isomorphic(graph, prior) for prior in bucket):
            return
        bucket.append(graph.copy())
        unfinished = [v for v in graph if graph.degree(v) < 3]
        if not unfinished:
            final.append(graph.copy())
            return
        # Different from the main checker: last vertex, not maximum degree.
        v = unfinished[-1]
        possible = [u for u in unfinished if u != v and not graph.has_edge(u, v)]
        for neighbours in itertools.combinations(possible, 3-graph.degree(v)):
            extended = graph.copy()
            extended.add_edges_from((v, u) for u in neighbours)
            visit(extended)

    visit(nx.empty_graph(10))
    assert len(final) == 21 and sum(nx.is_connected(g) for g in final) == 19
    print("independent cubic enumeration: 21 classes;", sum(map(len, seen.values())), "partial states")
    return final


def covers(graph):
    triples = [t for t in itertools.combinations(range(10), 3)
               if all(graph.has_edge(*e) for e in itertools.combinations(t, 2))]
    clauses = []
    for a, b in graph.edges():
        variables = [i+1 for i, t in enumerate(triples) if a in t and b in t]
        clauses.append(variables)
        clauses.extend([-u, -v] for u, v in itertools.combinations(variables, 2))
    with Cadical153(bootstrap_with=clauses) as solver:
        while solver.solve():
            chosen = [v for v in solver.get_model() if 0 < v <= len(triples)]
            assert len(chosen) == 10
            yield [triples[v-1] for v in chosen]
            solver.add_clause([-v for v in chosen])


def incidence(lines):
    graph = nx.Graph()
    graph.add_nodes_from((v, {"part": v//10}) for v in range(20))
    graph.add_edges_from((v, 10+j) for j, line in enumerate(lines) for v in line)
    return graph


def isomorphic(a, b):
    return nx.is_isomorphic(a, b, node_match=lambda x, y: x["part"] == y["part"])


def bareiss_certificate(case):
    lines = list(map(set, case["lines"]))
    adjacency = [set() for _ in range(30)]
    for j, line in enumerate(lines):
        for p in line:
            adjacency[p].update((10+j, 20+j))
            adjacency[10+j].add(p)
            adjacency[20+j].add(p)
    for j, k in itertools.product(range(10), repeat=2):
        if not lines[j] & lines[k]:
            adjacency[10+j].add(20+k)
            adjacency[20+k].add(10+j)
    free = iter(case["free_numerators"])
    matrix = [[0]*30 for _ in range(30)]
    for i in range(30):
        matrix[i][i] = 17500  # 10000 * (11/4 - 1)
        for j in range(i+1, 30):
            matrix[i][j] = matrix[j][i] = -10000 if j in adjacency[i] else -next(free)
    sentinel = object()
    assert next(free, sentinel) is sentinel
    previous = 1
    for i in range(30):
        pivot = matrix[i][i]
        assert pivot > 0  # successive leading principal minors (Sylvester)
        for j in range(i+1, 30):
            for k in range(i+1, 30):
                numerator = pivot*matrix[j][k]-matrix[j][i]*matrix[i][k]
                assert numerator % previous == 0
                matrix[j][k] = numerator//previous
        previous = pivot


def main():
    path = Path(__file__).with_name("construction30_configurations_certificates.json")
    data = json.loads(path.read_text())
    assert data["denominator"] == 10000 and data["bound"] == "11/4"
    saved = [incidence(case["lines"]) for case in data["cases"]]
    assert len(saved) == 10
    reached, count = set(), 0
    for leave in cubic_graphs():
        for lines in covers(nx.complement(leave)):
            count += 1
            graph = incidence(lines)
            matches = [i for i, candidate in enumerate(saved) if isomorphic(graph, candidate)]
            assert len(matches) == 1
            reached.update(matches)
    assert count == 20 and reached == set(range(10))
    print("independent SAT exact covers: 20 decompositions; all 10 saved classes covered")
    for case in data["cases"]:
        bareiss_certificate(case)
    print("independent integer Bareiss check: all 300 leading principal minors positive")
    print("Verified theta(complement(H)) <= 11/4 for every member of this family.")


if __name__ == "__main__":
    main()

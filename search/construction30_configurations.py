"""Exhaust the balanced linear-incidence generalization of H33.

A linear 10_3 configuration has ten points and ten three-point lines, each point
on three lines, with no pair of points on two lines. Construct H on three parts:
the points and two copies of the lines. Point-line edges are incidences; edges
between the line copies join disjoint lines. H is triangle-free and 6-regular.

All ten incidence-isomorphism classes are generated without a graph catalogue:
enumerate cubic leave graphs on ten points, then triangle-decompose each leave's
complement. Exact rational Lovasz dual certificates exclude theta(Hbar)=3 for
every class, allowing arbitrary probabilities, edge weights, and complex vectors.

  python search/construction30_configurations.py generate
  python search/construction30_configurations.py verify

Verification needs NetworkX and standard Python only, not an SDP solver.
Generation additionally needs the dependencies of search/dualcert.py.
"""
import argparse
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path


def cubic_leaves():
    """All unlabeled simple cubic graphs on ten vertices by partial completion.

    Every completion of a partial graph must supply the missing neighbours of
    the chosen unfinished vertex. We try all of them. Isomorphic partial graphs
    have identical sets of completed isomorphism classes, so keeping one is safe.
    WL hashes select buckets only: every deduplication uses exact isomorphism.
    """
    import networkx as nx
    seen, final = {}, []
    count = 0

    def visit(graph):
        nonlocal count
        key = nx.weisfeiler_lehman_graph_hash(graph, node_attr="color", iterations=3)
        bucket = seen.setdefault(key, [])
        if any(nx.is_isomorphic(graph, prior) for prior in bucket):
            return
        bucket.append(graph.copy())
        count += 1
        todo = [v for v in graph if graph.degree(v) < 3]
        if not todo:
            final.append(graph.copy())
            return
        v = max(todo, key=lambda x: graph.degree(x))
        candidates = [u for u in todo if u != v and not graph.has_edge(u, v)]
        for neighbours in itertools.combinations(candidates, 3-graph.degree(v)):
            graph.add_edges_from((v, u) for u in neighbours)
            visit(graph)
            graph.remove_edges_from((v, u) for u in neighbours)

    graph = nx.empty_graph(10)
    nx.set_node_attributes(graph, 0, "color")
    visit(graph)
    return final, count


def triangle_decompositions(graph):
    """All exact covers of the edges by triangles, with deterministic pivots."""
    edges = [tuple(sorted(e)) for e in graph.edges()]
    index = {edge: i for i, edge in enumerate(edges)}
    triples = [t for t in itertools.combinations(range(10), 3)
               if all(graph.has_edge(*e) for e in itertools.combinations(t, 2))]
    masks = [sum(1 << index[e] for e in itertools.combinations(t, 2)) for t in triples]
    containing = [[i for i, mask in enumerate(masks) if mask >> e & 1] for e in range(len(edges))]

    def visit(remaining, selected):
        if not remaining:
            yield sorted(triples[i] for i in selected)
            return
        options = min(([i for i in containing[e] if masks[i] & remaining == masks[i]]
                       for e in range(len(edges)) if remaining >> e & 1), key=len)
        for i in options:
            yield from visit(remaining ^ masks[i], selected+[i])

    yield from visit((1 << len(edges))-1, [])


def incidence_graph(lines):
    import networkx as nx
    graph = nx.Graph()
    graph.add_nodes_from((v, {"part": v//10}) for v in range(20))
    graph.add_edges_from((v, 10+j) for j, line in enumerate(lines) for v in line)
    return graph


def same_configuration(first, second):
    import networkx as nx
    return nx.is_isomorphic(first, second, node_match=lambda a, b: a["part"] == b["part"])


def enumerate_configurations():
    import networkx as nx
    leaves, states = cubic_leaves()
    configurations, incidences = [], []
    covers = 0
    for leave in leaves:
        for lines in triangle_decompositions(nx.complement(leave)):
            covers += 1
            incidence = incidence_graph(lines)
            if not any(same_configuration(incidence, prior) for prior in incidences):
                configurations.append(lines)
                incidences.append(incidence)
    stats = {"partial_cubic_states": states, "cubic_leaves": len(leaves),
             "connected_cubic_leaves": sum(nx.is_connected(g) for g in leaves),
             "triangle_decompositions": covers, "configuration_classes": len(configurations)}
    assert stats == {"partial_cubic_states": 308, "cubic_leaves": 21,
                     "connected_cubic_leaves": 19, "triangle_decompositions": 20,
                     "configuration_classes": 10}, stats
    return stats, configurations


def support(lines):
    """Rebuild and check a configuration and its 30-vertex nonorthogonality graph."""
    lines = [frozenset(line) for line in lines]
    assert len(lines) == len(set(lines)) == 10
    assert all(len(line) == 3 and line <= set(range(10)) for line in lines)
    assert all(sum(p in line for line in lines) == 3 for p in range(10))
    assert all(len(a & b) <= 1 for a, b in itertools.combinations(lines, 2))
    adj = [set() for _ in range(30)]

    def edge(a, b):
        adj[a].add(b)
        adj[b].add(a)

    for j, line in enumerate(lines):
        for p in line:
            edge(p, 10+j)
            edge(p, 20+j)
    for j, k in itertools.product(range(10), repeat=2):
        if lines[j].isdisjoint(lines[k]):
            edge(10+j, 20+k)
    assert all(len(ns) == 6 for ns in adj)
    assert all(v//10 != u//10 and not (adj[v] & adj[u]) for u in range(30) for v in adj[u])
    return adj


def positive_definite(matrix):
    """Exact symmetric elimination over Q; all certificates use strict slack."""
    matrix = [row[:] for row in matrix]
    for i in range(len(matrix)):
        pivot = matrix[i][i]
        if pivot <= 0:
            return False
        for j in range(i+1, len(matrix)):
            multiplier = matrix[j][i]/pivot
            for k in range(i+1, len(matrix)):
                matrix[j][k] -= multiplier*matrix[i][k]
    return True


def verify_case(case, denominator, bound):
    adj = support(case["lines"])
    free = iter(case["free_numerators"])
    matrix = [[Q(1) for _ in range(30)] for _ in range(30)]
    for i in range(30):
        for j in range(i+1, 30):
            if j not in adj[i]:
                matrix[i][j] = matrix[j][i] = Q(next(free), denominator)
    sentinel = object()
    assert next(free, sentinel) is sentinel
    slack = [[(bound if i == j else 0)-matrix[i][j] for j in range(30)] for i in range(30)]
    assert positive_definite(slack), "rational Lovasz dual certificate failed"


def generate(path):
    import networkx as nx
    from dualcert import dual_upper
    stats, configurations = enumerate_configurations()
    denominator, bound = 10000, Q(11, 4)
    cases = []
    for i, lines in enumerate(configurations):
        adj = support(lines)
        graph = nx.Graph()
        graph.add_nodes_from(range(30))
        graph.add_edges_from((u, v) for u in range(30) for v in adj[u] if u < v)
        nodes, _, matrix, optimum = dual_upper(graph)
        assert nodes == list(range(30))
        case = {"lines": lines, "theta_numeric": optimum,
                "free_numerators": [round(float(matrix[u, v])*denominator)
                                    for u in range(30) for v in range(u+1, 30) if v not in adj[u]]}
        verify_case(case, denominator, bound)
        cases.append(case)
        print(f"class {i+1}: numerical theta {optimum:.9f}, exact upper bound {bound}", flush=True)
    path.write_text(json.dumps({"scope": "all linear 10_3 configurations, identical line systems in B and C",
                               "enumeration": stats, "bound": str(bound), "denominator": denominator,
                               "cases": cases}, indent=2)+"\n")


def verify(path):
    data = json.loads(path.read_text())
    stats, configurations = enumerate_configurations()
    assert data["enumeration"] == stats
    assert len(data["cases"]) == len(configurations)
    saved = [incidence_graph(case["lines"]) for case in data["cases"]]
    for lines in configurations:
        incidence = incidence_graph(lines)
        assert sum(same_configuration(incidence, candidate) for candidate in saved) == 1
    bound = Q(data["bound"])
    assert bound < 3
    for case in data["cases"]:
        verify_case(case, data["denominator"], bound)
    print("exhaustive enumeration:", stats)
    print(f"all {len(configurations)} exact rational duals verified: theta(Hbar) <= {bound} < 3")
    print("No member of this family is a three-context GHZ graph, in any dimension or over C.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("generate", "verify"))
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("construction30_configurations_certificates.json"))
    args = parser.parse_args()
    if args.mode == "generate":
        generate(args.output)
    else:
        verify(args.output)

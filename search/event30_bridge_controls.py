"""Exact controls separating the two pair-component bridge obstructions.

Both controls are triangle-free, have three size-ten color classes, and every
vertex has a neighbor in each other class.  For the pair of contexts (0, 1):

* nonempty_shortage has bridge counts (4, 3, 3) for component side sizes
  (3, 3, 4), while every third-context vertex bridges some component;
* unbridged_only has bridge counts (3, 3, 4), meeting all three cardinality
  bounds, but vertex 20 bridges no component.

These are controls for the specified pair, not otherwise feasible candidates.
No pair component in either control is complete bipartite.
Run directly to verify the controls and the screen's exact certificates.
"""
import json


def controls():
    parts = [list(range(10 * k, 10 * k + 10)) for k in range(3)]
    pair_edges = []
    # Two C6 components, each represented as K3,3 minus its diagonal.
    for offset in (0, 3):
        pair_edges.extend((offset + i, 10 + offset + j)
                          for i in range(3) for j in range(3) if i != j)
    # One C8 component, with diagonal and next-diagonal pair edges.
    pair_edges.extend((6 + i, 16 + (i + shift) % 4)
                      for i in range(4) for shift in (0, 1))

    def graph(neighborhoods):
        edges = pair_edges + [(u, 20 + k) for k, row in enumerate(neighborhoods)
                              for u in row]
        return {"n": 30, "parts": parts,
                "edges": [list(edge) for edge in sorted(set(edges))]}

    # Three anticliques cover both sides of C8, leaving only three bridges
    # for a balanced component with four vertices on each side.
    shortage = [[i, 10 + i] for i in range(6)]
    shortage += [[6, 7, 19], [8, 16, 17], [9, 18], [0, 10]]

    # Initially one bridge per vertex on either side. Move vertex 20 across
    # components and let vertex 23 bridge both the first and second components.
    unbridged = [[i, 10 + i] for i in range(6)]
    unbridged += [[6 + i, 16 + (i + 2) % 4] for i in range(4)]
    unbridged[0] = [0, 13]
    unbridged[3] = [0, 3, 10, 13]
    return {"nonempty_shortage": graph(shortage), "unbridged_only": graph(unbridged)}


def independent_pair_data(graph):
    """Compute the relevant pair components independently of the main screen."""
    adjacency = [set() for _ in range(30)]
    for u, v in graph["edges"]:
        adjacency[u].add(v)
        adjacency[v].add(u)
    assert all(not (adjacency[u] & adjacency[v]) for u, v in graph["edges"])
    for part in graph["parts"]:
        for u in part:
            assert not adjacency[u].intersection(part)
            assert all(adjacency[u].intersection(other) for other in graph["parts"]
                       if other is not part)
    unseen, result = set(range(20)), []
    while unseen:
        todo = [min(unseen)]
        found = set(todo)
        while todo:
            new = (adjacency[todo.pop()] & set(range(20))) - found
            found.update(new)
            todo.extend(new)
        unseen.difference_update(found)
        left, right = found & set(range(10)), found & set(range(10, 20))
        bridges = [v for v in range(20, 30) if adjacency[v] & left and adjacency[v] & right]
        assert any(v not in adjacency[u] for u in left for v in right)
        result.append({"component": sorted(found), "left": sorted(left),
                       "right": sorted(right), "bridges": bridges})
    return result


def main():
    from overnight30_screen import certificate, pair_bridge_data, verify_certificate

    results = {}
    for name, graph in controls().items():
        data = independent_pair_data(graph)
        assert pair_bridge_data(graph, [0, 1]) == data
        covered = set().union(*(set(item["bridges"]) for item in data))
        shortage = [item for item in data if len(item["bridges"]) < len(item["left"])]
        if name == "nonempty_shortage":
            assert covered == set(range(20, 30))
            assert len(shortage) == 1 and len(shortage[0]["bridges"]) == 3
            cert = certificate(graph, "bridge_shortage", pair=[0, 1],
                               component=shortage[0]["component"], bridges=shortage[0]["bridges"])
        else:
            assert not shortage and set(range(20, 30)) - covered == {20}
            cert = certificate(graph, "unbridged_vertex", pair=[0, 1], vertex=20)
        verify_certificate(cert)
        results[name] = {"pair": [0, 1], "bridge_data": data, "certificate": cert}
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()

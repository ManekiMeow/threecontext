"""Reproduce the four colored automorphism orders using exhaustive VF2."""
import json
from pathlib import Path

import networkx as nx

DIRECTORY = Path(__file__).resolve().parent
records = {r['graph_id']: r for r in map(json.loads,
           (DIRECTORY/'strongest_candidates.jsonl').read_text().splitlines())}
audit = json.loads((DIRECTORY/'automorphism_audit.json').read_text())
for row in audit:
    graph = nx.from_graph6_bytes(records[row['graph_id']]['g6'].encode('ascii'))
    nx.set_node_attributes(graph, {v: v//10 for v in graph}, 'context')
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph, graph, node_match=lambda a, b: a['context'] == b['context'])
    count = sum(1 for _ in matcher.isomorphisms_iter())
    assert count == row['context_preserving_automorphism_order']
    print(row['graph_id'], count)

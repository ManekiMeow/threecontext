"""Exact full and context-preserving automorphism counts for two saved guides."""
import json
from pathlib import Path
import networkx as nx
if not __debug__:
    raise RuntimeError('Exact checks require assertions; do not use -O')
root=Path(__file__).resolve().parent
audit=json.loads((root/'automorphism_audit.json').read_text())
for role,row in audit.items():
    graph=nx.from_graph6_bytes(row['g6'].encode('ascii'))
    full=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter())
    nx.set_node_attributes(graph,{v:v//10 for v in graph},'context')
    preserving=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(graph,graph,
        node_match=lambda a,b:a['context']==b['context']).isomorphisms_iter())
    assert full==row['automorphisms_uncolored']
    assert preserving==row['automorphisms_context_preserving']
    print(role,full,preserving)

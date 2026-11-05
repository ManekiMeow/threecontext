"""Exact full and context-preserving automorphism counts for the new leader."""
if not __debug__:
    raise RuntimeError('Exact checks require assertions; do not use -O')
import json
from pathlib import Path
import networkx as nx
row=json.loads((Path(__file__).resolve().parent/'automorphism_audit.json').read_text())
graph=nx.from_graph6_bytes(row['g6'].encode())
full=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(graph,graph).isomorphisms_iter())
nx.set_node_attributes(graph,{v:v//10 for v in graph},'context')
preserving=sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(graph,graph,node_match=lambda a,b:a['context']==b['context']).isomorphisms_iter())
assert full==row['full_automorphism_group_order']==1
assert preserving==row['context_preserving_automorphism_group_order']==1
print('Full and context-preserving automorphism groups are both trivial.')

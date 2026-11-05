"""Independent exact nonisomorphism checks within the bounded selected set."""
if not __debug__:
    raise RuntimeError('Exact verification requires assertions; do not use -O')
import collections
import json
from pathlib import Path
import networkx as nx
HERE=Path(__file__).resolve().parent
buckets=collections.defaultdict(list)
prior_dates=('2026-11-04T2200','2026-11-04T2300','2026-11-05T0000','2026-11-05T0100')
for date in prior_dates:
    for row in map(json.loads,(HERE.parent/date/'strongest_candidates.jsonl').read_text().splitlines()):
        graph=nx.from_graph6_bytes(row['g6'].encode())
        buckets[nx.weisfeiler_lehman_graph_hash(graph)].append(graph)
count=0; comparisons=0
for row in map(json.loads,(HERE/'selected_representatives.jsonl').read_text().splitlines()):
    graph=nx.from_graph6_bytes(row['g6'].encode()); key=nx.weisfeiler_lehman_graph_hash(graph)
    for previous in buckets[key]:
        comparisons+=1
        assert not nx.is_isomorphic(graph,previous), row['graph_id']
    buckets[key].append(graph);count+=1
assert count==124
print(f'Verified {count} selected new representatives against one another and prior reviewed representatives; {comparisons} exact isomorphism comparisons after invariant hashing. No full-pool census.')

# Closed-Petersen pair: finite reduction and bounded pilots

Fix the allowed A--B support to `I+P`, where `P=KG(5,2)`. Every third-context
neighborhood is a zero rectangle of this support and can be enlarged to a
maximal rectangle. This preserves triangle-freeness and relaxes overlap
constraints. There are exactly **110 maximal templates**:

| A/B degrees | Count |
|---|---:|
| (1,6), (6,1) | 10 each |
| (2,4), (4,2) | 15 each |
| (2,3), (3,2) | 30 each |

Selecting ten distinct templates therefore covers all possible ten-event
third-context extensions. Repeated templates would give forbidden false twins.
Distinct maximal rectangles always give a length-three path between their C
vertices: otherwise their union would be a larger zero rectangle.

The exact scalar bound for a (2,4) event is **14/135**, not 4/39. The relaxed
masses `(2/9,7/18,7/18)` show why the latter proposed bound is invalid.
A singleton side gives the valid upper bound 2/27, so two singleton templates
are impossible by the context budget **8·14/135+2·2/27=44/45<1**. One singleton
is not excluded this way: the corresponding budget is **136/135>1**.

Two separate bounded SAT pilots use this reduction. They enforce ten selected
templates, at most one singleton template, original A/B-to-C degrees from one
through six, and every same-context contraction condition. The second pilot
also encodes perfect matchings to all ten selected templates on both sides.

| Pilot | Seconds | Representatives | Exact individual exclusions |
|---|---:|---:|---|
| Support constraints | 18.04 | 469 | 302 Hall; 167 probability propagation |
| Added perfect matchings | 28.09 | 280 | 216 probability propagation; 64 bridge shortage |

Every generated graph is exactly excluded. There were **zero numerical solves**.
The 749 representatives are distinct under the explicitly verified 240-element
action of the Petersen S5 permutations and A/B interchange, including across
the two pilots. They need not be pairwise nonisomorphic under arbitrary graph
isomorphisms and were not deduplicated against the entire prior search pool.
Neither pilot exhausted the family. **Closed-Petersen support without
isoclinicity remains unresolved.**

`summary.json` binds immutable raw archives in experiment chunk 104. Each
archive includes its discovery helper snapshot, selected templates, reconstructed
graphs, exact exclusion certificates, and, for the second pilot, both matching
witnesses. The verifier independently enumerates maximal rectangles from right
subsets, checks the declared action, reconstructs every graph, replays every
certificate, and checks the inter-pilot overlap count:

```sh
python search/event30_petersen_templates.py
python search/overnight30_petersen_reviews/2026-11-05T0200/verify.py
```

Use `--data-repo PATH` on the second command when the experiment data are in a
different checkout. Raw pilot records are auxiliary data and are not part of
the main supervisor's graph/score census.

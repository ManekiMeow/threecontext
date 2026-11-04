# Complete aligned two-kernel RTT exclusion

The group is C3×C3. Context A is its regular nine-point action plus a fixed
point; B and C each have three three-point orbits with quotient coordinates
`(x,x,y)` and one fixed point. All complex realizations and nonuniform event
probabilities are allowed. This is the aligned two-kernel case of Theorem 50.

The finite proof covers 364 SAT assignments. Of these, 299 have individual
exact structural certificates and yield positive absent-orbit clauses: a
realization supported in any subgraph of an excluded allowed support would
also realize that allowed support. The other 65 assignments have explicit
context-preserving maps to fourteen representatives, each with an exact
rational theta upper bound below three. The largest of those fourteen bounds
is `1714/625 = 2.7424`; this is not asserted as a uniform quantitative bound
on the structurally rejected assignments.

The final CNF contains the base constraints and exactly the 364 bound clauses.
Its LRAT proof derives the empty clause in 222 verified reverse-unit-propagation
steps. Thus every support satisfying the necessary base constraints is covered
by an exact exclusion. `coverage_manifest.json` binds the clauses, graphs,
certificates, vertex maps and proof files by SHA-256.

## Mathematical coverage of the base CNF

Two commuting translations generate all nine group elements. Their 54 complete
edge orbits partition all 300 possible cross-context edges; no orbit is omitted
or joined to a different orbit. The base CNF imposes triangle-freeness,
nonempty proper cross-neighborhoods, total degree at most ten, and distinct
neighborhoods in each directed pair of contexts (Proposition 39).

Each regular-to-three-orbit support is a singleton or its complement: an empty
or full block would give identical pair-neighborhood rows. Independent
translations of the six fiber coordinates normalize these sets to `{0}` or
`{1,2}` while leaving every B–C orbit variable free. This is a complete
relabeling reduction and does not constrain probabilities or overlap weights.

The additional orbit-degree bounds follow from strict probability budgets.
Moving three-orbits have cross-degrees between two and six and total degree
at most nine. The regular nine-orbit has cross-degrees between two and five
and total degree at most nine. Degree one bounds each probability below
`2/27`; the large opposite degrees and remaining total-degree-ten cases are
excluded by the harmonic bounds in `event30_orbit_degree.py`.

If both fixed edges `a-b` and `a-c` occur, singleton-neighborhood inequalities
force a to touch a moving orbit in each other context. The total-degree bound
then permits exactly one per context. These conditional clauses impose no
additional condition on either mixed fixed-point case. In particular, the
unproved C-side condition considered during research was never included.

Cardinality constraints use the reviewed standard sequential-counter encoding.
Repeated orbit literals count separate incident edges: substituting a common
variable into the ordinary occurrence-variable encoding preserves its
existential-extension property. XOR auxiliaries exactly encode row differences;
all auxiliary IDs are disjoint. Independent mathematical and code review found
no extra ansatz in these constraints.

## Reproduction

Standard Python checks the actual group action and complete edge-orbit partition,
all 313 exact certificates, the 65 vertex maps, all clause bindings, hashes,
and the full LRAT proof:

```
python search/overnight30_c3square_search.py --verify-coverage search/event30_c3square_aligned_certificates
```

Also regenerate the base CNF exactly using `python-sat==1.9.dev20`:

```
/workspace/threecontext-venv/bin/python search/overnight30_c3square_search.py --verify-coverage search/event30_c3square_aligned_certificates --rebuild-base
```

The recorded base has 9,855 variables and 24,582 clauses; the final CNF has
24,946 clauses. The mathematical interpretation above is separate from the
propositional UNSAT check. All proof assertions must be enabled; the LRAT
checker explicitly refuses optimized Python.

The original DRAT proof was independently verified with `drat-trim` commit
`2e3b2dc0ecf938addbd779d42877b6ed69d9a985`, then converted to the saved text
LRAT proof. That external tool is not required for replay. The small independent
checker `event30_lrat.py` accepts only positive-hint RUP steps and was reviewed
against 61,440 exhaustive two-variable applications, tiny satisfiable and
unsatisfiable controls, and tampered versions of the actual proof.

The broader bounded construction pool remains experiment data, separately
archived by the overnight publisher. This directory contains the complete
finite exclusion proof.

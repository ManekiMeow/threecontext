# Bounded expansion of the published C5 quotient review

The published checkpoint `bb92744f2b4fa70cb1a11b12ca4290bee9546f9e` contains
31,021 scored C5-double records in four completed pools. Their orbit counts
give 5,503 distinct quotient classes under context permutations and orbit
swaps. This expands the earlier three-pool census of 5,273 classes.

A 60-second analytic review skipped the previously selected 150 classes and
examined the next 3,560 classes in descending numerical score. It produced
**3,378 further exact quotient exclusions**:

- 502 rational sector-probability LP certificates;
- 1,299 exact affine/Cauchy saturation certificates;
- 1,577 probability LP certificates using disjoint-neighborhood moment bounds.

The remaining 182 reviewed classes were not excluded by these methods; another
1,793 published classes were not reached within the time limit. These are
counts within the published pool, not an exhaustive classification of the
C5-double family. All conditional exclusions require an actual automorphism
group with the stated orbits, and cover every cyclic edge phase of that
quotient.

The highest-scoring guide outside the earlier 150 classes has numerical theta
about **2.8894779990**, from chunk55, record index32565. Its quotient passed the
analytic attempts and its representative passed the structural screens used
in this review. The particular support already has an exact theta upper bound
`1445246758047/500000000000 < 3`; it is a mutation guide, not a GHZ realization.
Its individual dual and graph binding are preserved and verified here.

Proofs and per-class metadata are stored in compact gzip JSONL archives. Run:

```sh
python search/overnight30_quotient_expansion/2026-11-05T0000/verify.py
```

This uses only the standard library and replays every exact proof, checks each
saved representative's actual simultaneous C5 automorphism and quotient, and
checks the guide's individual exclusion. Numerical optimization was used only
for certificate discovery. `discover.py` reproduces a bounded review into a
new output directory and refuses to overwrite a completed catalog.

Worker integration replays every new proof before generating conditional C5
clauses. The new proofs took 4.01 seconds to verify; the remaining catalog took
4.57 seconds. All catalogs together generate 163,348 patterns and 165,652 SAT
clauses with 72 auxiliary variables. CNF construction and solver bootstrap took
0.72 seconds, with peak resident memory about 219 MiB. New examples from all
three proof methods were rejected, while the mutation guide remained admitted.

A fresh 20-second worker control generated 336 models and completed 81 scores.
Its sole reported `solver_error`, record 335, reads `deadline before numerical
solve`: this is a budget skip, with no substantive solver failure. A separate
resume control successfully replayed 47,691 historical models and 11,070
structural cuts, but its 14-second budget expired during the 68.55-second setup;
it added no new models or scores. Its cumulative three errors came from the
inherited checkpoint. These controls show both functional compatibility and
the growing cost of replaying accumulated historical certificates.

The benchmark JSON files record these measurements. Raw control outputs and
logs are preserved separately on the experiment branch, without changing its
campaign indexes, under
`search/overnight30_runs/2026-11-04/chunks/00076_asymmetric/auxiliary_quotient_expansion_benchmarks/`.
The auxiliary manifest hashes all 27 archived files and the compact archive.

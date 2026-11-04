# Multiplier-family classification, night of 30–31 Oct 2026 (RESULTS.md §3.13)

* `classify_eE.log` (`multiplier_classify.py E`): K4-free single-class patterns for every prime p ≡ 1 mod E up to the
  Weil threshold p₀(E) (complete for E = 2–6).
* `single_full_eE.log` (`single_full_run.py`): (C2), (C4), (TRI), (DEG) and ϑ = 3 for those patterns.
* `unions_e*.log` (`multiplier_classify.py unions`): union patterns as boxes in the K4-free set, with the survivors
  of the exact checks.
* `union_theta_e*.log` (`union_theta_run_e.py`): ϑ = 3 for the multi-class survivors, via the orbit-reduced
  `theta_orbits`.
* `iso_unions.py`: isomorphism classes of the PASS lines.
* `tmin_e*.log` (`tmin_e.py`): the flat bound t_min ≥ √h for the e = 5, 6 single-class candidates.

The helper scripts were run from the session's scratch directory, and their paths point there. To rerun, adjust the
SC paths.

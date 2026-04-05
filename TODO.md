# TODO

## [HOUSEKEEPING] 1. Flag duplicate scripts (NON-BLOCKING)

Rule 11 says no duplicate scripts. Known suspects:
- `experiments/analytic/weil_optimized.py` vs `weil_optimized_v2.py` vs `weil_optimized_v3.py`
- `experiments/sieve/ht_transfer_attempt.py` vs `ht_signed_transfer_v2.py`

Action: List duplicates in this file below. Human will review before deletion.
**This is non-blocking — proceed to research after listing them.**

Duplicates found (Session 39):
- [ ] `experiments/analytic/weil_optimized.py` vs `weil_optimized_v2.py` vs `weil_optimized_v3.py` (v3 same size as v1, possible revert)
- [ ] `experiments/sieve/ht_transfer_attempt.py` vs `ht_signed_transfer_v2.py`
- [ ] `experiments/circuit_complexity/k_party_nof.py` vs `k_party_nof_v2.py`
- [ ] `experiments/circuit_complexity/approx_degree.py` + `approx_degree_quick.py` + `approx_degree_small.py` + `approx_degree_prime.py` + `approx_degree_counting.py` (5 scripts; may be distinct experiments but violates one-script-per-experiment rule)
- [ ] `experiments/information_theory/information_shortcut.py` vs `information_shortcut_v2.py`
- [ ] `experiments/other/breakthrough_attempt_v2.py` (no v1 exists)
- [ ] `experiments/other/self_correcting_v2.py` (no v1 exists)
Human review needed before merging/deleting. All have companion _results.md files.

---

## [RESEARCH] 2. Benchmark Helfgott-Thompson O(x^{3/5}) for M(x)

`experiments/sieve/mertens_speedup.py` references H-T (2021). Transfer to pi(x) fails
(signed->unsigned barrier, S16), but the M(x) algorithm itself was never benchmarked.
Won't help pi(x), but establishes the constant-factor landscape.

---

## [RESEARCH] 3. Benchmark inversion_search.py fixed-point iteration

`experiments/sieve/inversion_search.py` implements p_{k+1} = R^{-1}(n + sum_rho R(p_k^rho)).
Fix import issues (depends on `riemann_explicit` module in experiments/analytic/) and
measure convergence rate for n = 10..10000.

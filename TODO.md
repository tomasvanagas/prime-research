# TODO

Priority: HOUSEKEEPING FIRST, then research. Do not start new experiments until
all items below marked [HOUSEKEEPING] are complete.

---

## [HOUSEKEEPING] 1. Close "Zeta Zero Compressibility" (OPEN_PROBLEMS.md #3)

Session 36 proved this is closed:
- 82% of Fourier modes needed for reconstruction (spectral_algebraic_structure_results.md)
- Kt grows linearly ~5.58*N (experiments/information_theory/kt_complexity/SYNTHESIS.md)
- No algebraic structure in spectrum

**Action:**
- In `status/OPEN_PROBLEMS.md`, move Problem #3 to CLOSED status with citation:
  "CLOSED Session 36: Kt(delta) grows linearly (5.58*N), 82% of spectral modes
  required, no algebraic structure. See kt_complexity/SYNTHESIS.md."
- Add entry to `status/CLOSED_PATHS.md` under Analytic section.

---

## [HOUSEKEEPING] 2. Backfill missing _results.md files

~200 Python scripts lack companion results files. This is a rule violation (Rule 10).
Run each script (or read its code + output comments) and write a short _results.md
with: (a) what it tested, (b) verdict (PASS/FAIL/PARTIAL), (c) key numbers.

Priority order (by gap severity):
1. `experiments/sieve/` — 19 scripts missing results (only 1 summary exists)
2. `experiments/proposals/` — 38 scripts missing results
3. `experiments/ml/` — 7 scripts, 0 results
4. `experiments/dynamical/` — 7 scripts, 0 results
5. `experiments/quantum/` — 4 scripts, 0 results
6. `experiments/topological/` — 5 scripts, 0 results
7. `experiments/analytic/` — ~45 scripts missing results
8. `experiments/circuit_complexity/` — ~81 scripts missing results
9. `experiments/wildcard/` — ~47 scripts missing results

**Strategy:** Use sub-agents. For each directory, spawn an agent that:
- Lists all .py files without a companion _results.md
- Reads each .py (docstring + any printed output/comments near __main__)
- If the script is runnable in <30s, run it and capture output
- If not runnable (missing deps, too slow), summarize from code comments
- Write a _results.md per script with: purpose, verdict, key findings
- If multiple scripts are trivial variants of the same idea, note that in results
  and flag for cleanup (Rule 11: no duplicate scripts)

**DO NOT skip this.** It is the single biggest gap in the project.

---

## [HOUSEKEEPING] 3. Write unified pseudorandomness summary

Sessions 28, 31, 35 collectively proved pi(x) mod 2 is random-like under 20+
independent structural measures. No single document lists them all.

**Action:** Create `novel/pseudorandomness_of_pi.md` listing ALL measures with
their values, the "random baseline" comparison, and session citations:

Measures to include (at minimum):
- Approximate degree = N/2 (S28)
- Communication rank = 2^{N/2-1}+2 (S17)
- GF(2) ANF sparsity = 0.50 (S35)
- SLP complexity = random (S35)
- Boolean Fourier: noise sensitivity ~0.9, total influence ~0.92 (S17)
- BDD size ~ 2^{0.73*N} (S20-28)
- PTF degree grows as N/2 (S35)
- SAT-based circuit size ~ N^9.8, matching random (S35)
- Linear complexity over GF(2..23) = 0.5000 (S24)
- Not k-automatic for any k (S24)
- LFSR complexity = N/2 (S28)
- k-party NOF rank = 2^{ceil(N/k)} for k>=3 (S23)
- Spectral flatness of zero contributions = 0.91 (S25)
- Entropy rate of delta(n) = 5.8 bits (S36)
- 1/f^1.69 spectrum, no algebraic structure (S20)
- Not holonomic / not D-finite (S23)
- Full tensor rank (S23)
- Per-bit influence: LSB-half ~2x MSB-half (S28)
- SVD spectral decay is power-law ~i^{-1} (S19)
- Compression ratio only 2-2.6x vs random (S35)

This is arguably the project's strongest collective finding.

---

## [HOUSEKEEPING] 4. Elevate completed Session 35 experiments

SAT-based circuit minimization and threshold circuit construction are COMPLETE
with results files, but still labeled "pending" in archive/ephemeral/.

**Action:**
- Read `experiments/circuit_complexity/sat_circuit_minimization_results.md`
- Read `experiments/circuit_complexity/threshold_circuit_construction_results.md`
- Add their findings to `status/SESSION_INSIGHTS.md` under Session 35
- Update `archive/ephemeral/critique_latest.md` to mark them as completed
- Add to CLOSED_PATHS if they close any approach

---

## [HOUSEKEEPING] 5. Document Lambert W approximation formula

`algorithms/v1_pade_approximation.py` contains a "Lambert W Prime Theorem" —
an O(1) approximation with 0.019% error via Pade [2,2] in W(n) basis.
Not mentioned anywhere in status files.

**Action:**
- Add to `status/BEST_ALGORITHMS.md` under a new "Approximation Methods" section
- Benchmark it against R^{-1}(n) for speed and accuracy
- Note: this is an APPROXIMATION (0.019% error), not exact

---

## [HOUSEKEEPING] 6. Run and document Proposal 21

`experiments/proposals/proposal21_zero_clustering_truncation.py` exists but was
never run and has no results file.

**Action:**
- Run it (should be <30s with 200 zeros)
- Write `proposal21_zero_clustering_truncation_results.md`
- Add verdict to CLOSED_PATHS if it fails (likely: clustering cannot reduce
  zero count to polylog because phases decohere at large x)

---

## [HOUSEKEEPING] 7. Flag duplicate/orphan scripts for cleanup

Rule 11 says no duplicate scripts. During the results backfill (item 2),
flag any scripts that are trivial variants of each other. List them here
for manual review before deletion.

Known suspects:
- `experiments/sieve/ht_transfer_attempt.py` vs `ht_signed_transfer_v2.py`
- Multiple proposal scripts that may duplicate closed paths
- Any `*_v2.py`, `*_quick.py`, `*_small.py` patterns

---

## [RESEARCH] 8. Investigate Helfgott-Thompson O(x^{3/5}) benchmark

`experiments/sieve/mertens_speedup.py` references H-T (2021) O(x^{3/5+eps})
for M(x). The transfer to pi(x) was proven to fail (signed->unsigned barrier),
but the H-T algorithm for M(x) itself was never benchmarked.

**Action:**
- Implement H-T algorithm for M(x) and benchmark against Lucy DP
- Document the signed->unsigned transfer barrier clearly
- This won't help pi(x) but establishes the constant-factor landscape

---

## [RESEARCH] 9. Benchmark inversion_search.py fixed-point iteration

`experiments/sieve/inversion_search.py` implements p_{k+1} = R^{-1}(n + sum_rho R(p_k^rho)).
It was never run to completion or documented.

**Action:**
- Fix any import issues (depends on `riemann_explicit` module)
- Run for n = 10..10000, measure convergence rate
- Document whether iteration converges, how many steps, and total cost
- Write results file

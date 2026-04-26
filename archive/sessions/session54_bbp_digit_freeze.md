# Session 54 — Fresh-perspective wildcard: BBP-style digit spigot

**Date:** 2026-04-26
**Mode:** Fresh thinking, no reading of CLOSED_PATHS / CLAUDE.md
**Outcome:** Two clean nulls. No breakthrough.

## Premise
Computing pi(n) for huge n is information-theoretically tied to the
behavior of zeta zero contributions. The BBP formula isolates the k-th
hex digit of pi *without* computing earlier digits — could pi(x) (or
psi(x)) admit anything similar?

## Brainstorm (5 ideas, ranked by testability)
1. **BBP-style digit-freeze** for psi(x) partial sums — TESTED.
2. **Random-subset compressed sensing** of zeros — TESTED.
3. **Hierarchical multipole sieve** (FMM analog) — not coded; sieve has
   no quadratic kernel to compress.
4. **Galois descent for pi(x) mod q** — already covered by S35
   pseudorandomness findings.
5. **Spectral low-pass smoothing** pi*K_T — equivalent to (1).

## Experiment 1: bbp_digit_freeze.py
Compute psi(x) approximations from K=10..2000 zeros for x=1e4..1e7.
Measure residuals; track digit stabilization; hunt for anomalous x.

**Result:** Residuals are Gaussian with std proportional to sqrt(x).
Going from K=10 to K=2000 reduces error only ~4-300x, never enough
to add a digit per *constant factor* of K. Anomaly density (1 in 4001
samples below |res|=0.012) matches Gaussian-tail prediction — no
structural shortcut.

**Implication:** Zeros-per-digit grows ~10^d (geometric in digits),
not O(d) as a BBP-style spigot would require. Standard `sqrt(x)`
information bound holds.

## Experiment 2: random_zero_subset.py
Compare first-K=50 vs random-K=50 from a pool of M=2000 zeros.

**Result:** Random subsets are 4x WORSE (rms 404 vs 97). Even the
best of 50 random trials (rms 310) is 3x worse than first-K. Reason:
amplitudes are 1/|rho_k|-weighted, concentrating mass on low-frequency
zeros that random sampling misses.

**Implication:** Compressed-sensing-style "K random measurements
recover the signal" doesn't apply — the signal is dense, not sparse,
in the zero basis.

## Status updates
- Appended both null results to `status/CLOSED_PATHS.md`.
- Appended Session 54 entry to `status/SESSION_INSIGHTS.md`.
- Two scripts + companion `_results.md` saved in `experiments/wildcard/`.

## Literature check (2026)
WebSearch returned no new prime-counting algorithm publications.
State-of-art remains Deléglise-Rivat-Gourdon at O(x^{2/3}). FKBJ
analytic method at O(x^{1/2+eps}). No polylog candidate.
EOF
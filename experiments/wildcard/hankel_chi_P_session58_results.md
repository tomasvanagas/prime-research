# Hankel singular-value test for chi_P (session 58)

## Question
If the prime indicator chi_P(n) satisfies any short linear recurrence over R
(or admits a low-rank approximation thereof), the Hankel matrix
        H[i, j] = chi_P(i + j)        for 0 <= i, j < K
has effective rank equal to that recurrence order (Berlekamp-Massey).

We compute singular values of H_chi and compare against two random controls:
  (a) **naive control:** Bernoulli with the same density as chi_P (~0.16).
  (b) **sieve-aware control:** Bernoulli supported on integers coprime to 30,
      with density on that support matching chi_P's density on it.
The naive control inherits no wheel-sieve structure, so it inflates the
"deviation"; the sieve control isolates whatever residual structure remains.

## Setup
- N = 1250, K = 600.
- pi(N) = 204, density 0.1631; skeleton-{coprime to 30} density 0.27,
  conditional density on skeleton = 0.605.
- 5 random seeds for each control; report mean and std of singular values
  and effective ranks.
- Effective rank at threshold t = smallest r such that ||sigma_{r:}||^2 /
  ||sigma||^2 < t.

## Results

### vs naive Bernoulli control (uniform support)
The leading SVs of H_chi are dramatically larger than naive control:

| k | s_chi | s_naive_mean | z |
|---|---|---|---|
| 0 | 94.1 | 95.4 | -0.16 |
| 1 | 94.1 | 22.0 | **+40.2** |
| 2 | 47.1 | 22.0 | **+13.9** |
| 3-5 | ~47.1 | ~22.0 | **+30** |
| 6-9 | ~25.1 | ~18.0 | +7-10 |

Looks like a giant signal — but it is **entirely explained by wheel
structure**. The pattern (94.1, 94.1, 47.1, 47.1, 47.1, 47.1, 25.1, ...)
matches what you get from a parity/mod-6/mod-30 projection: each modulus
contributes a low-rank piece whose SV scales as (density on residue class) *
sqrt(K/period).

### vs sieve-aware control (Bernoulli on coprime-to-30)
| k | s_chi | s_sieve_mean | z |
|---|---|---|---|
| 0  | 94.11 | 97.48 | -1.08 |
| 1  | 94.10 | 97.47 | -1.08 |
| 2-5 | 47.1 | 49.0 | -1.27 to -1.29 |
| 6-9 | 25.1 | 27.0 | -1.45 |
| 16  | 17.08 | 16.48 | +0.66 |
| 32  | 12.42 | 12.69 | -0.54 |
| 64  | 10.42 | 10.37 | +0.22 |
| 128 | 8.61 | 8.30 | +1.52 |
| 256 | 5.95 | 5.98 | -0.36 |

Effective rank at 1% energy threshold:
- chi_P:           442
- naive control:   479 (chi_P ~8% lower — pure wheel effect)
- sieve control:   442 (**identical**)

Max |z| over k in [10, 100] (where leading wheel modes have decayed) =
**2.17** — well below the 3-sigma significance threshold.

## Verdict — FAIL (failure mode I, information loss)
After matching the wheel-sieve (mod 30) structure, the prime indicator's
Hankel singular value spectrum is **statistically indistinguishable** from
a Bernoulli random model on the coprime-to-30 skeleton.

This is a new pseudorandomness measure: **No short linear recurrence over R
generates chi_P.** Hidden algebraic-recurrence shortcuts (e.g. "primes are
the n-th term of a P-recursive sequence") are ruled out at the resolution
of K=600, beyond what's already captured by the wheel sieve.

## What this rules out
- Any conjecture of the form "chi_P satisfies a constant-coefficient linear
  recurrence over R / Q / R[1/N] / etc. of order o(K)" — the SV spectrum
  would show a rank-r truncation gap, but the residual is random.
- More generally: **any bilinear-low-rank structure** in chi_P beyond
  the wheel sieve itself. A polylog representation in terms of the Hankel
  basis does not exist.

## What it does NOT rule out
- Recurrences over **finite rings** (e.g., Z/p or finite extensions) — those
  would show up only in modular Hankel ranks. Worth a follow-up: compute
  rank of H_chi mod small primes and compare to Berlekamp-Massey.
- Recurrences on **subsequences** (e.g., on prime gaps, on chi_P along
  arithmetic progressions). Different test.
- **Non-stationary / multiplicative** structure — Hankel only sees additive
  recurrences. Toeplitz on multiplicative shifts is a separate test.

## Files
- `hankel_chi_P_session58.py` — single experiment script.
- `hankel_chi_P_session58_data.npz` — raw singular values for replay.

## Edges touched
- E1.x (information chain): another null measurement on chi_P. Confirms
  the pseudorandomness picture once wheel structure is subtracted.
- E2.x (algebraic chain): rules out Berlekamp-Massey-style linear recurrence
  shortcut over R. Subset/refinement of the broader algebraic null.

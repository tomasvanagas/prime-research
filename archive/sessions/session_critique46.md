# Session Synthesis — Critique-46 (2026-04-26)

**Mode:** critique. Verifies the four S46-fresh proposals
(`archive/ephemeral/proposals_session.md`) against CLOSED_PATHS.md.

## Summary

The proposer drafted four concrete proposals with runnable code, ran
each, and self-reported all four as CLOSED with the right failure mode.
The critique confirms each proposer self-verdict and adds DUPLICATE-PLUS
verdicts citing the closest CLOSED_PATHS parents.

## Verdicts

| # | Proposal | Critic verdict | Mode | Closest prior |
|---|----------|----------------|------|----------------|
| A | Wynn-epsilon (= Shanks transform) on psi(x) zero-sum | DUPLICATE-PLUS | I | lines 26, 49, 657, 697 |
| B | Borel-Pade resummation of Cipolla asymptotic for p(n) | DUPLICATE-PLUS | E | lines 40, 43, 49 |
| C | Depth-5 multiplicative fingerprint primality oracle | DUPLICATE-PLUS | C | line 689 |
| D | Zero-aware control variate for Monte-Carlo pi(x) | DUPLICATE-PLUS | E+I | lines 256, 257, 688 |

## Key empirical findings

### Proposal A — Wynn-epsilon catastrophically diverges

Geometric ladder T_k = 5*2^k, k=0..5, K_max=160 zeros, x in {100, 1000,
10000}.

| x | best partial err | Wynn extrap err | Wynn worse by |
|---|---|---|---|
| 100 | 0.024 | 78.18 | 3000x |
| 1000 | 0.525 | 995.42 | 1900x |
| 10000 | 2.43 | 10013.88 | 4100x |

Empirical |psi - psi_T| log-log slope b ≈ -0.25 average (RH predicts -1).
Wynn-eps_2 = Shanks transform (already tested by name in line 26);
P-A is a strict subcase that confirms the closure under one more
specific transform.

### Proposal B — Borel-Pade strictly worse than raw Cipolla

Cipolla truncation K=1..7 + Borel-Pade resummation at L = log n,
M = log log n.

| n | best Cipolla err | Borel-Pade err | penalty |
|---|---|---|---|
| 10 | 5.97 (K=1) | 85.67 | 14x |
| 100 | 27.77 (K=2) | 204.62 | 7.4x |
| 1000 | 78.6 (K=2) | 534.05 | 6.8x |
| 10000 | 183.1 (K=3) | 260.33 | 1.4x |

Cipolla relative error 0.21 → 0.0017 from n=10 to 10⁴ (asymptotic-series
character preserved), but absolute error grows. Diagnosis: Stokes-line
or multi-instanton structure that Borel-Pade in single direction can't
capture; OR asymptotic series is past optimal truncation point at small K.

### Proposal C — depth-5 fingerprint achieves perfect separation in [2, 10⁴]

Fingerprint = (tau(n)-n^11-1 mod 691, tau(n) mod 5, sigma_1(n)-1-n,
n mod 30, n^6 mod 7).

| depth | unique prime fps | composite FPs | FP rate | prime collide rate |
|---|---|---|---|---|
| 1 | 1 | 10 | 0.11% | 100% |
| 2 | 3 | 9 | 0.10% | 74.86% |
| 3 | 3 | 0 | 0.00% | 0.00% |
| 4 | 11 | 0 | 0.00% | 0.00% |
| 5 | 12 | 0 | 0.00% | 0.00% |

USEFUL EMPIRICAL ANCHOR but not novel: structural barriers are
identical to line 689's parent closure — tau(n) for composite n requires
factoring n via multiplicativity (subexp via GNFS, not polylog), and
even with a polylog primality oracle, counting primes in [1,x] is
pi(x) again (circular at meta level). Depth-5 = 26 bits per n; finite-
domain perfect separation collapses asymptotically.

### Proposal D — PNT control variate variance reduction is 1.006x

| T | err@x=100 | err@x=1000 | err@x=10000 |
|---|---|---|---|
| 5 | 0.77 | 1.74 | 12.22 |
| 50 | 0.40 | 2.42 | 9.85 |
| 500 | 0.05 | 0.04 | 5.62 |

x=10000, T=200 zeros: PNT-control-variate var reduction 1.006x;
naive MC err 13.2 vs PNT-controlled 14.4 (slightly WORSE). Log-log slope
of |psi-psi_T| vs T: average b ≈ -0.29 (small-T regime; at T >= sqrt(x)
slope approaches -1 = info-theoretic limit).

## Pattern observation

The convergence-acceleration / variance-reduction family of interventions
on already-closed analytic primitives is now **systematically exhausted**
across sessions 5, 6, 10, 15, 25, 32, 43-46, 48, 51, 63, and critique-46.
Specific transforms catalogued: Pade, Wynn (= Shanks), Shanks (named),
Richardson, Aitken Δ², Cesaro, Fejer, Borel (single-direction and Borel-
Pade), Mellin-Barnes, Hermite mollification, Gaussian mollification,
Riesz mollification, Selberg-style Dirichlet-polynomial mollification,
PNT control variate, Harper random-multiplicative, randomized zeta
zero sampling, randomized I-E. Future proposals in this family should
be pre-disqualified via a single CLOSED_PATHS check.

The bar for any future acceleration-style proposal: must invoke a
NONLINEAR operation on zeros (line 693 explicit ask) AND must not be
mathematically equivalent to a transform already on the master list.

## Project status

Unchanged. `OPEN_PROBLEMS.md` still has Circuit Complexity of pi(x) as
the only genuinely viable direction. Within Chain E:

- Sub-attack 1 (Bernstein 2003 strengthened-r AKS) — CLOSED S66
- Sub-attack 2 (non-cyclotomic ring AKS) — CLOSED S61
- Sub-attack 3 (Healy-Viola Frobenius transplant) — CLOSED S64
- Computationally cornered.

Remaining levers: Brandt MKtP (FOCUS-3, un-engaged) and pure
new-technique work (open problem class).

## Files modified

- `archive/ephemeral/critique_latest.md` (full critique, ~5800 chars)
- `status/CLOSED_PATHS.md` (4 new entries appended, line 723-726;
  header updated to "541+ approaches", "Sessions 1-66 + critique-46")
- `status/SESSION_INSIGHTS.md` (critique-46 section appended)
- `archive/sessions/session_critique46.md` (this file)
- `.run_state` → 46

# Trace Formula Approach: Computing pi(x) Without Enumerating Zeros

**Date:** 2026-04-05
**Verdict:** CLOSED -- All four trace/moment strategies fail fundamentally.

## Motivation

The explicit formula pi(x) = li(x) - sum_rho li(x^rho) + ... looks like a trace:
Tr(f(A)) where A has eigenvalues {rho}. In random matrix theory, Tr(f(A)) can
sometimes be computed from matrix entries via moments Tr(A^k), avoiding eigenvalue
enumeration. Can this shortcut work for zeta zeros?

## Experiment 1 & 5: Moments Method

**Idea:** Expand S(x) = sum_j x^{rho_j}/rho_j as a power series using spectral
moments M_k = sum_j rho_j^k.

**Result: DIVERGES.**

| k | |M_k| | k! | |M_k|/k! |
|---|-------|-----|----------|
| 0 | 1.00e+03 | 1 | 1.00e+03 |
| 5 | 1.14e+18 | 120 | 9.53e+15 |
| 10 | 3.65e+33 | 3.63e+06 | 1.01e+27 |
| 20 | 6.43e+64 | 2.43e+18 | 2.64e+46 |
| 30 | 1.46e+96 | 2.65e+32 | 5.50e+63 |

Moments grow as |M_k| ~ gamma_max^k ~ 1419^k, vastly outpacing k!.
The series radius of convergence is x < e^{1/gamma_max} ~ 1.0007, i.e., useless.

For all test x in {100, 1000, 10000}, the moment expansion diverges violently:
- At K=5: error ~ 10^10
- At K=20: error ~ 10^54
- At K=80: error ~ 10^182

**Why it fails:** Zeta zeros are unbounded (gamma_j -> infinity), so spectral
moments grow without bound. GUE matrices have bounded eigenvalues, making their
moment expansion convergent. This is a fundamental structural difference.

## Experiment 2: GUE Random Matrix Trace

**Idea:** For GUE(N), Tr(e^{itH}) converges nicely via matrix power series.
Does this carry over to the zeta case?

**Result: GUE ANALOGY IS MISLEADING.**

| Matrix | t | K needed for 10^{-6} accuracy |
|--------|---|-------------------------------|
| GUE(50) | 1.0 | ~10 |
| GUE(50) | 5.0 | ~25 |
| GUE(50) | 10.0 | ~40 |
| GUE(100) | 1.0 | ~10 |
| GUE(200) | 1.0 | ~10 |

The GUE moment expansion converges beautifully (to machine precision by K~50)
because eigenvalues are in [-1, 1]. But zeta zeros have gamma_j up to 1419
(with 1000 zeros) and growing without bound. The bounded-eigenvalue assumption
that makes GUE traces work is violated.

## Experiment 3: Weil Explicit Formula (Geometric Side)

**Idea:** The Weil formula relates sum_rho h(rho) to sums over primes.
Can we evaluate the prime side using only partial information?

**Result: CIRCULAR.**

| y/x fraction | psi(y)/psi(x) at x=10000 |
|--------------|--------------------------|
| 1% | 0.020 |
| 10% | 0.107 |
| 50% | 0.502 |
| 80% | 0.798 |
| 100% | 1.000 |

To compute psi(x) within +/-0.5 from the geometric side requires knowing
essentially ALL primes up to x. The explicit formula trades spectral zeros
for primes -- neither side avoids enumerating the other. This is a fundamental
duality, not a computational shortcut.

## Experiment 4: Effective Rank

**Idea:** Build matrix M[i,j] = 2Re[x_i^{rho_j}/rho_j] and check if it has
low-rank structure exploitable for fast computation.

**Result: NO USEFUL COMPRESSION.**

For x in [1000, 10000] with 200 zeros:

| Rank k | Max |error| | Relative error |
|--------|-------------|----------------|
| 1 | 36.3 | 92% |
| 10 | 28.8 | 73% |
| 50 | 12.2 | 31% |
| 100 | 4.9 | 12% |
| 150 | 0.9 | 2.3% |
| 200 | ~0 | ~0 |

Energy concentration:
- 99% energy: rank ~107 (of 200)
- 99.9% energy: rank ~139
- 99.99% energy: rank ~159

The operator is essentially full-rank. To get S(x) within O(1), we need
rank comparable to the number of significant zeros -- no compression advantage.

## Experiment 6: Smoothed Trace (Gaussian Damping)

**Idea:** Damp high zeros with e^{-gamma^2/T^2} to make moments converge,
then remove the damping.

**Result: UNFAVORABLE TRADE-OFF.**

At x = 100: T = 250 (295 effective zeros) needed for |error| < 0.5.
At x = 1000: No T in [10, 2000] achieved |error| < 0.5 with 1000 zeros.
At x = 10000: Error ~ 37 even with all 1000 zeros undamped.

The smoothing reduces the number of effective zeros, but the remaining zeros
are insufficient for exact pi(x). Making T large enough for accuracy means
including essentially all zeros, defeating the purpose.

## Overall Verdict

**CLOSED.** The trace formula approach to computing pi(x) without enumerating
individual zeta zeros fails for four independent reasons:

1. **Moment divergence:** Spectral moments M_k ~ gamma_max^k grow too fast for
   the moment expansion to converge (radius of convergence ~ 1.0007).

2. **Unbounded spectrum:** The GUE analogy breaks down because GUE eigenvalues
   are bounded while zeta zeros are unbounded. The trace-from-entries method
   relies fundamentally on bounded eigenvalues.

3. **Weil duality is circular:** The "geometric side" of the explicit formula
   requires knowing all primes up to x, which is what we are trying to compute.

4. **Full-rank operator:** The zero-sum operator has no useful low-rank structure.
   Approximating S(x) to O(1) accuracy requires rank ~ number of contributing zeros.

These four failures are not independent workarounds -- they reflect the same
underlying reality: the oscillatory part of pi(x) has intrinsic information content
that cannot be bypassed by algebraic rearrangement of the explicit formula.

**Add to CLOSED_PATHS.md:** Trace formula / moment method / GUE trace shortcut
for zero sums -- divergent moments, circular duality, full-rank operator.

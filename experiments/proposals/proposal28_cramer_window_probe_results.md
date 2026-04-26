# Proposal 28 Results — Empirical Cramér-Window Probe

## Question
Is the window |p_n − R⁻¹(n)| narrow enough (polylog-scale) that, given an
anchor count π(R⁻¹(n) − win), we could finish localising p_n in a sieved
window of polylog work?

## Setup
- Sieved primes up to 1.6 × 10⁷ (covers p(10⁶) = 15 485 863).
- 80 log-spaced samples of n in [10, 10⁶].
- Computed R⁻¹(n) by Newton's method on R(x) = Σ_{k=1..30} μ(k)/k · Li(x^{1/k}),
  Li via series (x > 1).
- Tabulated δ_n = p_n − R⁻¹(n), |δ_n|/log²(p_n), |δ_n|/√(p_n), and a
  log-linear fit |δ| ~ p^α.

## Numerical findings (n ≤ 10⁶)

| metric                      | value     |
|-----------------------------|-----------|
| max |δ_n|                   | 1823.0    |
| max |δ_n| / log²(p_n)       | **6.65**  |
| mean |δ_n| / log²(p_n)      | 0.66      |
| max |δ_n| / √(p_n)          | 0.59      |
| empirical scaling α (|δ| ~ p^α) | **0.505** |

The fit α ≈ 0.505 matches the RH-predicted √x growth almost exactly. The
ratio |δ_n|/log²(p_n) is *not* bounded — it grew from < 1 at n ≤ 1000 to 6.65
at n = 10⁶ and would continue to ≈ √(p)/log²(p) → ∞.

## Verdict — CLOSED (failure mode: equivalence/circularity)

The proposal conflated two different bounds:
- Cramér's conjecture bounds **prime gaps** g_n = p_{n+1} − p_n by O(log²p_n).
- The R⁻¹(n) approximation error is **Ω(√p_n / log p_n)** under RH — it is
  the oscillatory contribution of the Riemann zeros, not a prime-gap quantity.

So the window required to bracket p_n around R⁻¹(n) is Ω(√p_n), not polylog.
A polylog-window sieve does not localise p_n given only R⁻¹(n).

The "anchor counting" bottleneck is therefore the original problem: any
fast π evaluation at any anchor would already give nth-prime in O(T_π log n).

## Failure mode
Equivalence (E): the proposal reduces to "compute π(x) on a polylog-anchor
set", which is the original problem. The (min,+) machinery does not change
this because the sieve work is dominated by anchor cost, not within-window
cost.

## Useful empirical constants
- For n ≤ 10⁶, **|δ_n| ≤ 0.59 √(p_n)** (better than RH's 1/(8π) bound for
  these n).
- The 99-th-percentile ratio |δ_n|/log²(p_n) ≈ 4 within tested range — would
  matter if we could solve anchor counting separately.

These constants are consistent with values reported by prior sessions.

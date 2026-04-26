# P4 — Iterated Newton with progressive zero-budget

**Script:** `session63fresh_newton_zerobudget.py`

## What was tested

Newton iteration to solve π(x) = n converges quadratically *if* the
function evaluation is accurate enough. Hypothesis: a *progressive*
zero budget — K_k = 2^(k+1) zeros at iteration k — gives sufficient
accuracy because Newton tolerates O(error²) noise in early steps,
and total zero work is geometric (≤ 2 · K_final).

## Method

Initial guess x₀ = R⁻¹(n) (Newton on R alone, no zeros). At iteration k:
1. Use K_k = 2^(k+1) zeros to compute π_K(x) = R(x) − Σ 2 Re Ei(ρ log x).
2. Compute π_K'(x) likewise (analytic derivative).
3. Newton step: x_{k+1} = x_k − (π_K(x_k) − n) / π_K'(x_k).
4. Track residual = x_k − p_n_truth.

Tested n ∈ {100, 1000, 5000, 10000} with K up to 256.

## Result

| n | p(n) truth | x₀ residual | iterations | K_total | final residual |
|---|-----------|-------------|-------------|---------|----------------|
| 100 | 541 | −4.52 | 2 | 6 | **+0.016** |
| 1000 | 7919 | +3.57 | 8 | 510 | −0.95 (oscillates) |
| 5000 | 48611 | −56.0 | 8 | 510 | **+15.99 (DRIFTS)** |
| 10000 | 104729 | +38.8 | 8 | 510 | +18.85 (oscillates) |

For n = 100, the method works beautifully: starting from R⁻¹(100) ≈
536.5 (residual −4.5), 2 Newton iterations with 6 total zeros locks
onto x = 541.016 — within ±0.5 of the true 541.

For n ≥ 1000, the iterates **oscillate** in a band of size 5–25 around
the true p(n) and never close in. For n = 5000, the iteration even
diverges *away* from the better starting point.

## Why it fails (root cause)

The Newton update assumes pi_K(x) is a smooth approximation that gets
*better* as x → p(n). But pi_K(x) has *high-frequency oscillations* of
amplitude O(√x / log x) (sum-of-zeros tail) that don't disappear as x
moves; they merely change shape. With K = 2 zeros you've smoothed
those oscillations to ~constant, but the constant is wrong; with
K = 256 you've added high-frequency oscillations whose troughs and
peaks are O(1) — the same scale as inter-prime gaps.

Quantitatively: at x ≈ 10000, the prime gap is ~ log(10000) ≈ 9.
So |π_K(x) − π(x)| stays at order 1 for all reasonable K, and
Newton's update direction (pi_K(x) − n) is essentially noise of
amplitude 1 in the *target*. Newton then jumps log(x) ≈ 9 in x per
spurious unit of π-noise. Result: the iterate bounces in a band of
size ~ 10–30, which is what we observe.

The premise that "residual halves per Newton iteration" is wrong:
quadratic convergence requires the function value's error to decay
to zero with x, not stay at O(1).

## Verdict

**CLOSED.** Failure mode: **Information loss (I)** — the explicit-
formula truncation has irreducible noise of order 1 at any K, which
Newton amplifies to order log(x) in the iterate. For n ≥ 1000, the
iterate cannot localize p(n) to ±0.5 by zero summation alone.

The method *does* succeed for n ≤ ~500 where R⁻¹(n) is already within
~5 of p(n) and 2-zero correction does the rest. This bounds the regime
where the method is useful at sub-O(p(n)^{2/3}) cost.

## What would change the verdict

If we could reduce |π_K(x) − π(x)| below the prime-gap scale (i.e.,
below log(p(n))) using only polylog zeros, this method would work.
That requires sub-Riemann-Siegel zero accuracy or a different
expansion. None is currently known.

## One-line summary

Progressive zero-budget Newton converges in 2 iterations for n=100
(error 0.016 with 6 zeros total) but fails for n ≥ 1000 because the
zero-tail noise in π_K(x) is ~1, the same scale as prime gaps —
Newton then oscillates with amplitude log(x) in the iterate.

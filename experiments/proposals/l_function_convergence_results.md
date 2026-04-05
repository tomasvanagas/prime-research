# L-function Zero Sum Convergence vs Riemann Zeta Zero Sum -- Results

**Script:** l_function_convergence.py

## What Was Tested
Does the explicit formula for pi(x; q, a) (primes in arithmetic progressions) converge faster with L-function zeros than pi(x) converges with zeta zeros? If so, CRT-based reconstruction might be cheaper per modulus.

## Key Findings
- Computed E(x) = pi(x) - Li(x) and E(x; q, a) = pi(x; q, a) - Li(x)/phi(q) for various q and x up to 10^6.
- The roughness (variation) of E(x; q, a) is comparable to E(x) when normalized by 1/phi(q). No smoothing advantage.
- L-function zeros for chi mod q are distributed similarly to zeta zeros (same density ~ T/(2*pi) * log(T/(2*pi))). No sparsity gain.
- The explicit formula for pi(x; q, a) involves 1/phi(q) * sum_chi chi(a) * sum_rho terms. Each character contributes its own zero sum.
- For fixed q, the total number of zeros needed across all characters is phi(q) * O(sqrt(x)) -- more work, not less.
- Conclusion: L-function zero sums for arithmetic progressions do NOT converge faster; the per-modulus cost is at least O(sqrt(x)).

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- L-function zero sums have the same density and convergence rate as zeta zero sums; no per-modulus savings.

## One-Line Summary
L-function vs zeta convergence: same zero density, same O(sqrt(x)) cost per modulus; no advantage for arithmetic progression counting.

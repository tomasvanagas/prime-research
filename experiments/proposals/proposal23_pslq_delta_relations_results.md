# Proposal 23: PSLQ/LLL Integer Relation Discovery for delta(n) — Results

## Idea
Use PSLQ/LLL integer relation algorithms to search for algebraic relationships between delta(n) and known mathematical constants (pi, log values, Euler gamma, zeta values, etc.). If delta(n) has a closed form in terms of computable quantities, we get O(polylog).

## Key Results

### Polynomial fit of delta(n)
All polynomial bases in log-quantities give R² ≈ 0.001 or less:
- log(n), log(log(n)): R² = 0.0002
- sqrt(p)/log(p): R² = 0.0000
- 1, log(n), n^{1/3}/log(n): R² = 0.0004

**delta(n) is essentially uncorrelated with all smooth functions of n**.

### PSLQ integer relations
PSLQ found relations at specific n values (error < 10^{-6}), but they are **different for each n**:
- n=100: -5δ ≈ -14·sqrt(p)/ln(p)² + 10·π - 5·π/ln(n) + ln(2)·ln(n)
- n=500: δ ≈ -11 + 4·π + 12·π/ln(n) - γ·ln(n)
- n=1000: -6δ ≈ 4 + 4·ln(n) - 4·ln(n)·ln(ln(n)) - 3/ln(n)

Different relations at different n = **overfitting**, not a genuine formula.

### Normalized delta statistics
- delta_hat(n) = delta(n) · log(p(n)) / sqrt(p(n))
- Mean: -0.357, Std: 2.157, Range: [-7.57, 6.50]
- **Passes normality test** (p=0.19 > 0.05) — consistent with Gaussian distribution
- Autocorrelation at lag 1: 0.96 (smooth), but decays to 0.39 by lag 50

## Verdict: CLOSED
The normalized delta is approximately Gaussian with slowly decaying autocorrelation. PSLQ finds spurious point-wise relations that don't generalize. No universal algebraic formula for delta(n) exists in the tested constant space. The Gaussian distribution further confirms that delta encodes random-like information from zeta zeros.

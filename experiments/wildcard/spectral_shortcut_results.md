# Spectral Shortcut: Results

**Script:** spectral_shortcut.py

## What Was Tested
Whether the zero sum sum_rho Li(x^rho) can be computed as a "trace" without enumerating individual zeros: convergence vs number of zeros, Gaussian-smoothed variants, structural patterns in sum x^{i*gamma}/(1/2+i*gamma), short-formula approximability, and Riemann-von Mangoldt N(T) verification.

## Key Findings
- Explicit formula convergence: error decreases as ~x^{1/2}/K for K zeros used -- confirmed standard theory
- Gaussian smoothing: reduces oscillation but introduces bias; unbiased smoothing still needs O(x^{1/2}/ln(x)) zeros
- Zero sum structure: no exploitable pattern in partial sums; they oscillate without converging to a smooth limit
- Short-formula search: no combination of < 10 elementary functions approximates the zero sum to < 0.5 error
- Riemann-von Mangoldt N(T) ~ (T/2pi) * ln(T/2pi) - T/2pi: verified, but doesn't help compute individual zeros

## Verdict
**CLOSED**
**Failure Mode:** Equivalence (E) -- the zero sum has no trace-formula shortcut; convergence rate is the standard O(x^{1/2}/K) requiring K = O(x^{1/2}/ln(x)).

## One-Line Summary
Spectral shortcut (trace formula for zero sum): no trace structure found; convergence O(x^{1/2}/K) confirms standard explicit formula; no shortcut.

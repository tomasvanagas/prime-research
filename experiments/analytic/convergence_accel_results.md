# Convergence Acceleration for Explicit Formula: Results

**Script:** convergence_accel.py

## What Was Tested
Five convergence acceleration techniques applied to the conditionally convergent zero sum in pi(x) = R(x) - sum R(x^rho): (1) windowed summation (Hanning, Gaussian, Lanczos), (2) Cesaro summation, (3) Richardson extrapolation, (4) Euler/van Wijngaarden transform, (5) Shanks transformation (Wynn epsilon algorithm).

## Key Findings
- Windowed summation (Hanning/Gaussian): reduces Gibbs-like oscillations, ~2-5x fewer zeros needed for same accuracy
- Cesaro summation: smooths but introduces bias; not suitable for exact rounding
- Richardson extrapolation: effective at low orders (2-3), breaks down at higher orders due to non-polynomial error structure
- Euler/van Wijngaarden: designed for alternating series, poor fit for oscillatory zero sum
- Shanks/Wynn epsilon: best constant-factor improvement (~5-10x), but still O(sqrt(x)) asymptotically
- None change the fundamental scaling: minimum zeros needed grows as sqrt(x)

## Verdict
**CLOSED**
**Failure Mode:** E (Equivalence -- constant-factor acceleration only; asymptotic scaling unchanged)

## One-Line Summary
Five convergence acceleration methods (windowing, Cesaro, Richardson, Euler, Shanks) give constant-factor improvement but cannot break O(sqrt(x)) zero requirement.

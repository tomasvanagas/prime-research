# Differential Equation Search: Results

**Script:** diffeq_search.py

## What Was Tested

Whether the smoothed remainder f(x) = pi(x) - R(x) satisfies any ODE or integral equation. Five sections: (1) Gaussian-smoothed f(x) and derivatives at multiple scales, (2) linear ODE search via SVD, (3) multiplicative (Euler-type) ODE search, (4) non-linear ODE search, (5) Volterra integral equation test.

## Key Findings

- Gaussian smoothing at sigma=5,10,20,50 produces clean derivatives but no ODE structure emerges
- Linear ODE search (SVD): no combination of f, f', f'', f''' with polynomial coefficients gives residual below 1%
- Euler-type ODE (x*f' + ... form): same negative result
- Non-linear ODE search: no improvement over linear; residuals remain large
- Volterra integral equation: no kernel found that maps f to a simpler function
- f(x) does NOT satisfy any simple ODE or integral equation, consistent with encoding oscillatory zeta-zero contributions not annihilated by any fixed-order DE

## Verdict

**CLOSED** -- Failure Mode: Information Loss (I). f(x) is not annihilated by any fixed-order linear/nonlinear differential or integral equation; the oscillatory content is too rich.

## One-Line Summary

No ODE or integral equation found for f(x) = pi(x) - R(x) -- the oscillatory zeta-zero content cannot be captured by any fixed-order differential operator.

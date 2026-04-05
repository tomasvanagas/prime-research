# Inversion Search: Results

**Date:** 2026-04-04 (Session context)
**Script:** inversion_search.py

## What Was Tested
Exact formula for the nth prime via inverting the Riemann R function. Six parts: (1) li^{-1}(n) as approximation, (2) correction functions g(n), (3) R^{-1}(n) inversion, (4-5) fixed-point iteration p_{k+1} = R^{-1}(n + sum_rho R(p_k^rho)), (6) analysis of the zero-sum correction.

## Key Findings
- li^{-1}(n) overestimates systematically; error grows as sqrt(p(n))*ln(p(n)) consistent with RH
- R^{-1}(n) is much better than li^{-1}(n) but still not always roundable to correct prime
- Fixed-point iteration converges with 500 zeta zeros: exact by rounding at k=0 or k=1 for most tested n up to 10000
- The correction sum_rho R(x^rho) is oscillatory, grows as ~sqrt(x)/ln(x), cannot be replaced by elementary function
- The iteration is a direct consequence of inverting pi(x) = R(x) - sum_rho R(x^rho) -- not a new discovery per se
- Computing the zero sum requires O(sqrt(x)) zeta zeros for exact results -- same complexity barrier

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss) -- the correction requires O(sqrt(x)) zeta zeros, encoding incompressible oscillatory information

## One-Line Summary
R^{-1} inversion with fixed-point iteration works but requires O(sqrt(x)) zeta zeros for exactness -- no polylog shortcut.

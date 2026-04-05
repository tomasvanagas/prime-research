# Pattern Mining: Results

**Script:** pattern_mining.py

## What Was Tested
Mining hidden patterns in the prime sequence to find a direct formula p(n) = f(n). Includes: residual analysis after known approximation, residual scaling (best alpha fitting), spectral analysis (FFT of residuals), autocorrelation of residuals, higher-order difference analysis, and systematic search for algebraic/transcendental relations via curve fitting on 100K primes.

## Key Findings
- Residual r(n) = p(n) - approx(n) scales as ~n^(1/3) on average, with std growing as ~n^(1/2)
- Best scaling exponent alpha ~ 0.45 for r(n)/n^alpha having constant variance
- FFT of residuals shows broad spectrum with no dominant frequencies (spectral flatness ~0.9)
- Autocorrelation of gaps drops to ~0 by lag 2 -- gaps are essentially uncorrelated
- Finite differences of all orders grow without bound -- p(n) is not polynomial
- No algebraic relation found between p(n) and simple functions of n beyond known asymptotics

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- residuals after smooth approximation are spectrally flat and uncorrelated; no hidden pattern exists)

## One-Line Summary
Systematic pattern mining on 100K primes finds no exploitable structure in residuals; spectral flatness ~0.9 confirms random-like behavior.

# AlphaEvolve-Inspired Formula Evolution: Results

**Script:** formula_evolution.py

## What Was Tested
Genetic programming to evolve symbolic expressions for p(n) using building blocks: n, log, exp, sqrt, floor, ceil, Lambert W, li, R, R_inv, and arithmetic operators. Strategy: p(n) = known_approx(n) + evolved_correction(n). Population of expression trees mutated/crossed over generations, evaluated against first 10000 primes.

## Key Findings
- Best evolved formulas rediscover known asymptotic expansions (n*ln(n) + n*ln(ln(n)) - n + corrections)
- Correction terms evolved by GP have train RMSE ~2-5 but test RMSE ~10-50 (severe overfitting)
- No formula achieves >10% exact match rate on held-out test data
- The evolved "corrections" are effectively curve-fitting noise in delta(n), which does not generalize
- Increasing population size and generations increases train accuracy but worsens generalization ratio

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- delta(n) is information-theoretically incompressible; no finite symbolic expression can capture it)

## One-Line Summary
Evolutionary formula search (AlphaEvolve-style GP) rediscovers known asymptotics; correction terms overfit and do not generalize.

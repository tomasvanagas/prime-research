# Symbolic Regression / Genetic Programming: Results

**Script:** symbolic_regression.py

## What Was Tested
Symbolic regression via custom genetic programming search and parametric formula optimization (differential evolution). Trained on n=1..5000, tested on n=5001..10000. Parametric families: n*(ln(n) + ln(ln(n)) + a + b/ln(n) + c*ln(ln(n))/ln(n)^2 + ...) with up to 8 free parameters. GP evolves expression trees with arithmetic, log, exp, sqrt, floor.

## Key Findings
- Best parametric formula achieves ~0.003% relative error but only ~2-5% exact match rate on test
- Differential evolution finds optimal constants that match known asymptotic coefficients
- GP expression trees converge to equivalent of known expansions plus noise-fitting terms
- Train exact match: up to ~8%; test exact match: ~3-5% (severe overfitting of correction)
- Generalization ratio (test RMSE / train RMSE) consistently 2-5x across all GP runs
- No evolved formula beats R_inv(n) in exact match rate on test data

## Verdict
**CLOSED**
**Failure Mode:** I (Information Loss -- parametric optimization and GP both converge to known asymptotics; residual correction overfits)

## One-Line Summary
Parametric optimization and GP symbolic regression converge to known asymptotic formulas; correction terms overfit with 2-5x generalization gap.

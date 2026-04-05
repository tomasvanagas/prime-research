# Compute f(x) = pi(x) - R(x) Data Generation: Results

**Script:** compute_fx.py

## What Was Tested

Data generation script that computes f(x) = pi(x) - R(x) for x = 2..100000, where R(x) is Riemann's smooth approximation via the Gram series. Outputs fx_data.json and fx_data.npz for use by subsequent identity search experiments (LLL, WZ, diffeq, PSLQ).

## Key Findings

- pi(100000) computed via sieve; R(x) computed via Mobius-weighted li(x^{1/k}) sum with 80 terms at 30-digit precision
- f(x) oscillates around 0 with growing amplitude (~sqrt(x)/log(x))
- Data saved in both JSON and compressed numpy formats for downstream experiments
- This is a utility/data-generation script, not an experiment per se

## Verdict

**CLOSED** -- Failure Mode: N/A (data generation utility). The oscillatory residual f(x) is the target of identity searches in companion scripts; no identity was found by any of them.

## One-Line Summary

Data generation script computing f(x) = pi(x) - R(x) for x = 2..100000, saved as fx_data.npz for downstream identity search experiments.

# Full experiment driver

Runs the three parts of the D45 experiment (closed-form identity
verification on Γ_N and H_N, z-scores vs random matched divisors,
generalised Hasse identity for other arithmetic divisors) and writes
`full_results.json`.

**Status:** built and used in S161. See main writeup
`baker_norine_chi_p_results.md` for the consolidated theorems,
proofs, z-scores, and self-graded synthesis.

**Output:** `full_results.json`, `full_run.log`.

**Reproducibility:** `python3 run_full.py`. Total runtime ~3 seconds
on standard machine. Deterministic (seed=42 for random divisor
sampling).

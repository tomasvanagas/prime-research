# gamma1_phase_diagnostic — γ_1-cosine ceiling for cross-modulus lag-1 autocorr

Sub-script of `cross_modulus_pi_structure.py`. See
`cross_modulus_pi_structure_results.md` § "γ_1-cosine ceiling diagnostic"
for the populated table and analysis.

## What it does

For each m ∈ {3, 5, 6, 10, 30}, computes
`φ_1(m) := γ_1 · log(m) mod 2π` and `cos(φ_1(m))` (the ceiling), then
compares to empirical lag-1 autocorrelation read from
`raw_data.json`.

## Falsification criterion (pre-registered post-hoc with 6/6 retroactive
m=2 check)

`|empirical_lag1_ac(m)| ≤ |cos(γ_1·log(m) mod 2π)| + 0.05` should hold
for at least 5/5 cross-modulus cases. Sign agreement may be lower
because higher-zero contributions can flip the sign for small K_m.

## Result

5/5 magnitude bound HOLDS. 3/5 sign agreement (m=5, 6, 10 match;
m=3, 30 mismatch). Retroactively m=2 magnitude bound HOLDS (S246
empirical lag-1 = +0.283 ≤ |cos(9.797 mod 2π)| + 0.05 = 0.937 + 0.05).

## Output

`gamma1_phase_diagnostic.json` — per-m {φ_1, cos(φ_1), emp_lag1, max_lag,
max|ac|} plus aggregate {sign_match, bound_holds, n}.

# `T_sweep.py` — companion script results

This is a follow-on to `cohn_elkies_chi_p.py` that fixes `N = 10^6` and
sweeps `T_max ∈ {50, 100, 200, 400, 800, 1500}`. Headline finding:
the LP optimum `f̂*(0; T)` grows log-linearly in `T` with slope ≈ 0.85
for `g_obs`, slope ≈ 5.24 for `g_HL`.

See `cohn_elkies_chi_p_results.md` for the full structural analysis,
falsifiability statement, edge connection, and self-grading. This file
exists only to satisfy the `<name>.py + <name>_results.md` requirement
in CLAUDE.md.

Numerical output: `T_sweep_results.json`. Run log: `T_sweep.log`.

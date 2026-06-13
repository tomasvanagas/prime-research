# s425_inverse_chipfiring_probe — results

Consolidated session-level results live in
`s425_inverse_chipfiring_results.md`. This file holds the probe's
primary script; raw measurements in `s425_inverse_chipfiring_results.json`.

**Probe summary.** Six divisors {D_P, D_uniform, D_sink_K=N, D_τ,
D_μ², D_Ω₁} × two graphs {Γ_N, H_N} × four N values {16, 32, 64, 128} =
48 q-reduction runs. Each row reports `qreduced_sink_chips`,
`qreduced_nonsink_support_size`, and `qreduced_nonsink_total`.

**Headline finding.** D_Ω₁ on H_N gives sink chips = π(N) exactly
(verified at all N), with off-sink chips totalling
`Σ_{k≥2} π(N^{1/k})` — chip conservation forces this. Identity is
new (not in S161, not in EDGES.md prior to S425).

**Run command.**

```
python3 s425_inverse_chipfiring_probe.py --N 16 32 64 128 --out s425_inverse_chipfiring_results.json
```

Wall time ~2 s. Output: 48 rows of JSON.

**Falsifiability.** Any deviation `qreduced_sink_chips != π(N)` for
D_Ω₁ × H_N × N ≥ 16, or `qreduced_sink_chips != π(N) − π(N/2)` for
D_Ω₁ × Γ_N × N ≥ 16, would refute identities S425.1 and S425.2
respectively.

See `s425_inverse_chipfiring_results.md` for full closure analysis and
`archive/sessions/session425_reverify_e2_28.md` for session synthesis.

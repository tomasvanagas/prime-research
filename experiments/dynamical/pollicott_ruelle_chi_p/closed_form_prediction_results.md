# closed_form_prediction.py — Rayleigh-quotient analytical formula for `λ_0^h`

Implements the closed form `T_n = (1/log²2) ∫₀¹ dx/[(1+x)(x+n)(x+n+1)]`
via partial fractions:

  `T_n = (1/log²2)·[ln 2/(n(n-1)) − ln((n+1)/n)/(n-1) + ln((n+2)/(n+1))/n]`

for `n ≥ 2`, with `T_1 = (1/log²2)·[−ln 2 + 1/2 + ln(3/2)]`.

Defines `a_n = T_n / ‖g‖²` (where `‖g‖² = (2 log²2)^{-1}`), with
`Σ_{n≥1} a_n = 1`. Predicts `λ_0^h ≈ Σ_n h(n) a_n`.

**Numerical agreement at `n_max = 400`:**

| `h` | predicted | measured (M=120, n=400) | rel. error |
|-----|-----------|-------------------------|------------|
| `1` (unweighted) | Σ a_n = 0.99655 | 0.99640 | +0.015% |
| `χ_P` (primes) | Σ_p a_p = **0.36187** | **0.35961** | **+0.629%** |
| `Λ` (von Mangoldt) | 0.5206 | 0.4968 | +4.78% |
| `λ` (Liouville, signed) | 0.1749 | 0.0902 | (signed; Rayleigh fails) |

**`χ_P` prediction holds to better than 1%.** This is the headline
analytical result: the χ_P-weighted Pollicott-Ruelle leading resonance
admits an explicit closed form `Σ_{p prime} a_p` with `a_n` an
arithmetic-content-independent Gauss-Kuzmin Rayleigh-quotient
coefficient sequence.

Asymptotic: `a_n ~ 2 log 2 / n²` for large `n`, so `Σ_p a_p ~
2 log 2 · P(2)` where `P(2) = Σ_p 1/p² ≈ 0.4523` is the prime zeta at
`s = 2`. The leading-order asymptotic recovers `2 log 2 · 0.4523 ≈ 0.627`,
but small-prime corrections (`a_2 = 0.170`, `a_3 = 0.092`, `a_5 = 0.040`)
matter quantitatively.

**For Λ**: prediction is +4.8% off, larger than for χ_P. Likely cause:
the Λ-weighted operator carries small contributions from prime-power
indices `p^k` with `k ≥ 2` that have weights `log p` rather than 0
(unlike χ_P which is 0 at composites including prime powers). The
leading right eigenvector deviation from `g` is correspondingly
larger.

**For λ (signed)**: Rayleigh-on-`g` predicts 0.175, measured 0.090
— factor-of-2 discrepancy. The reason: signed cancellation makes the
LEFT eigenvector deviate substantially from constant. Rayleigh
quotient `R(g) = ⟨g, L_λ g⟩/⟨g, g⟩` is only `≈ λ_0` when `g` is close
to BOTH right AND left eigenvectors. For unsigned positive weights
(χ_P, Λ), the left eigenvector is approximately constant, so the
formula works to sub-5%; for signed weights, it fails.

For full discussion + numerics see `pollicott_ruelle_chi_p_results.md`.

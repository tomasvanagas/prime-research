# §D7 — DPP / Permanantal-PP fit to the integer prime sequence

**Session:** 95 (frontier attack on ATTACK_VECTORS §D7).
**Mode:** ambitious frontier attack with B-grade fallback.
**Self-grade:** **B-grade** structural negative result + quantitative
breakdown identity, plus extension to signed-real-K and complex-
Hermitian-K cases (both also fail).
**Channel:** Tao (additive combinatorics & point-process framework).
**Cross-domain technique imported:** Determinantal Point Processes
(DPPs) and the permanantal counterpart, from
Hough-Krishnapur-Peres-Virag 2009 *Zeros of GAFs and DPPs*
(AMS ULect 51) and Soshnikov 2000 (arXiv:math/0002099). First time
this technique appears in the project (was PROPOSED in
`CROSS_DOMAIN_TECHNIQUES.md §3`, never USED).

---

## Pre-stated falsification statements (written BEFORE running)

**F1 (DPP pair-level):** for every admissible even t > 0,
`K^2_DPP(t) := rho^2 - R_2(t) < 0`. DPP fails at the pair level.

**F2 (PPP pair-level):** for every odd t > 1,
`K^2_PPP(t) := R_2(t) - rho^2 < 0`. PPP fails at the pair level.

**F3 (PPP 3-point identity):** for some admissible all-even triple
`(0, t_1, t_2)`, the relative gap
`|R_3^PPP - R_3^HL| / R_3^HL > 0.10`.

**F4 (real-signed K):** even with sign assignments `s : t → {-1, +1}`
in `K(t) = s(t) ρ √(S(0,t) - 1)`, the cross-term ratio
`σ_req := [R_3^HL - ρ^3 - ρ Σ K(t_i)^2] / [2|K(t_1)||K(t_2)||K(dt)|]`
is NOT in {-1, +1} at any tested triple.

**F5 (complex Hermitian K):** even with complex Hermitian K of fixed
magnitudes |K(t)| = ρ √(S(0,t) - 1), no global phase assignment
`φ : t → [0, 2π)` makes `cos(φ(t_1) + φ(dt) - φ(t_2)) = σ_req` hold
across all admissible triples within tolerance.

---

## Mathematical setup (compact)

A *translation-invariant simple point process* on `Z` with intensity
`ρ` and real-symmetric kernel `K(t) = K(-t)` has 2-point inclusion
probability:

  `R_2(t) = K(0) K(t) - K(t)^2  =  ρ^2 - K(t)^2`     (DPP, det form)
  `R_2(t) = K(0) K(t) + K(t)^2  =  ρ^2 + K(t)^2`     (PPP, perm form)

setting `K(0) = ρ` (since `R_1 = K(0,0) = ρ`).

DPPs are **repulsive** (`R_2 ≤ ρ^2`); PPPs are **attractive**
(`R_2 ≥ ρ^2`).

For the prime indicator `χ_P` on `Z`:

  - **t even, admissible** (Hardy-Littlewood): `R_2(t) ≈ ρ^2 S(0, t)`
    with `S(0, t) > 1` (twin / cousin / sexy primes have positive
    density). **DPP requires `K^2 ≥ 0`, but `ρ^2 - R_2(t) < 0`.**
    Pair-level DPP infeasible.

  - **t odd, > 1**: `R_2(t) ≈ 0` (odd-offset prime pairs include only
    (2, 2 + odd_t) — finitely many). PPP requires
    `K^2 = R_2(t) - ρ^2 ≥ 0`, but `0 - ρ^2 < 0`. Pair-level PPP
    infeasible.

So R_2 - ρ^2 changes sign across odd/even, killing both real DPP and
real PPP at the pair level. **F1 and F2 are structural certainties**
from Hardy-Littlewood + parity of primes.

The deeper question is whether even the **3-point** identity is
self-consistent on the all-even sub-lattice: if we restrict to all
even offsets and try a PPP with magnitudes `|K(t)|^2 = R_2(t) - ρ^2 ≥ 0`,
does the permanent identity match the empirical / HL 3-point?

---

## Run configuration

  - `N = 10^7`, `π(N) = 664,579`, `ρ = 6.6458 × 10^-2`, `log N = 16.118`.
  - Pair offsets `t ∈ [0, 30]` → `R_2(t)`; one-shot numpy multiplies.
  - Triple list: 20 all-even admissible triples + 6 parity-mixed.
  - HL singular series truncated at `P_max = 5000` (669 primes).
  - Wall time: 5.3 s end-to-end.

Run script: `primes_dpp_ppp_fit.py`. Raw data: `main_run.json`.

---

## Pair-level results (F1, F2)

```
  t  parity     R_2_emp        R_2_HL    R_2/ρ²   S_HL    K²_DPP    K²_PPP
  1     odd  1.000e-07  0.000e+00  0.00002  0.000   +4.42e-3   -4.42e-3
  2    even  5.898e-03  5.831e-03  1.33540  1.320   -1.48e-3   +1.48e-3
  3     odd  1.000e-07  0.000e+00  0.00002  0.000   +4.42e-3   -4.42e-3
  4    even  5.862e-03  5.831e-03  1.32729  1.320   -1.45e-3   +1.45e-3
  6    even  1.172e-02  1.166e-02  2.65375  2.641   -7.30e-3   +7.30e-3
  8    even  5.860e-03  5.831e-03  1.32668  1.320   -1.44e-3   +1.44e-3
 10    even  7.821e-03  7.776e-03  1.77082  1.760   -3.40e-3   +3.40e-3
 12    even  1.175e-02  1.166e-02  2.66007  2.641   -7.33e-3   +7.33e-3
 14    even  7.046e-03  6.999e-03  1.59539  1.584   -2.63e-3   +2.63e-3
 16    even  5.861e-03  5.831e-03  1.32693  1.320   -1.44e-3   +1.44e-3
 18    even  1.175e-02  1.166e-02  2.65955  2.641   -7.33e-3   +7.33e-3
 20    even  7.822e-03  7.776e-03  1.77098  1.760   -3.41e-3   +3.41e-3
 22    even  6.532e-03  6.480e-03  1.47895  1.467   -2.12e-3   +2.12e-3
 24    even  1.173e-02  1.166e-02  2.65681  2.641   -7.32e-3   +7.32e-3
 26    even  6.420e-03  6.362e-03  1.45366  1.440   -2.00e-3   +2.00e-3
 28    even  7.044e-03  6.999e-03  1.59485  1.584   -2.63e-3   +2.63e-3
 30    even  1.565e-02  1.555e-02  3.54379  3.521   -1.12e-2   +1.12e-2
```

**Empirical R_2(t) matches HL prediction `ρ^2 S(0, t)` to better than 1.5%
at every even t**, confirming Hardy-Littlewood at this scale.

**F1 (DPP pair-level) HOLDS:** `K²_DPP < 0` for **all 15** admissible
even t in [2, 30]. No pair-level DPP fit is possible.

**F2 (PPP pair-level) HOLDS:** `K²_PPP < 0` for **all 14** odd t > 1
in [3, 29] (and for t=1). No pair-level PPP fit is possible.

Together these say the kernel sign cannot be uniform across odd/even,
so neither real DPP nor real PPP can fit χ_P globally — the parity of
the prime sequence is incompatible with translation-invariant kernel
sign.

---

## 3-point results (F3): PPP overshoots HL by up to 79%

Restricting to all-even triples (where `K^2_PPP > 0` is feasible at the
pair level), we compute:

  - `R_3^emp(t_1, t_2) = (1/N) Σ_n χ_P(n) χ_P(n + t_1) χ_P(n + t_2)`,
  - `R_3^HL(t_1, t_2) = ρ^3 S(0, t_1, t_2)`,
  - `R_3^PPP(t_1, t_2) = perm[K]` with `|K(t)| = ρ √(S(0, t) - 1)`,
    positive sign convention.

```
        triple    R_3_emp     R_3_HL    R_3_PPP   PPP-HL%   emp-HL%
     (0,2,6)    8.543e-04  8.390e-04  1.224e-03   +45.84    +1.82
     (0,2,8)    8.627e-04  8.390e-04  1.223e-03   +45.79    +2.82
     (0,2,12)   1.292e-03  1.259e-03  1.490e-03   +18.40    +2.68
     (0,4,6)    8.677e-04  8.390e-04  1.224e-03   +45.84    +3.42
     (0,4,10)   1.273e-03  1.259e-03  1.480e-03   +17.63    +1.16
     (0,4,14)   1.000e-07  0.000e+00  1.018e-03      ----    ----
     (0,4,16)   8.630e-04  8.390e-04  1.220e-03   +45.44    +2.86
     (0,6,8)    8.515e-04  8.390e-04  1.223e-03   +45.79    +1.49
     (0,6,12)   1.719e-03  1.678e-03  3.002e-03   +78.93    +2.47
     (0,6,16)   1.279e-03  1.259e-03  1.480e-03   +17.61    +1.61
     (0,6,18)   1.725e-03  1.678e-03  3.006e-03   +79.16    +2.79
     (0,8,14)   1.053e-03  1.049e-03  1.383e-03   +31.82    +0.42
     (0,8,18)   1.288e-03  1.259e-03  1.482e-03   +17.78    +2.35
     (0,10,16)  1.290e-03  1.259e-03  1.480e-03   +17.61    +2.54
     (0,10,22)  1.455e-03  1.416e-03  1.607e-03   +13.52    +2.73
     (0,12,18)  1.698e-03  1.678e-03  3.006e-03   +79.16    +1.18
     (0,14,20)  1.604e-03  1.573e-03  1.691e-03   + 7.52    +1.95
     (0,16,22)  9.636e-04  9.439e-04  1.314e-03   +39.23    +2.09
     (0,18,24)  1.715e-03  1.678e-03  3.004e-03   +79.03    +2.22
     (0,20,26)  1.411e-03  1.384e-03  1.585e-03   +14.48    +1.95
```

**Two structural patterns visible:**

1. **Empirical R_3 matches HL prediction to better than 6%** at every
   admissible triple (consistent with finite-N noise + the truncated
   singular-series numerics; HL is correct).

2. **R_3^PPP overshoots R_3^HL** (and hence empirical) by **+7.5% to
   +79.2%**, NEVER less. The maximum gap is on **3-arithmetic-progression
   triples** `(0, 6, 12), (0, 12, 18), (0, 18, 24), (0, 6, 18)` where
   `|K(t_1)| = |K(dt)| = |K(t_2)|` are all large. For these the
   permanent's cross-term is maximal.

3. **Inadmissibility blindness:** the triple `(0, 4, 14)` is
   inadmissible mod 3 (residues `{0, 1, 2}` cover all of `Z/3Z`, so for
   any `n` exactly one of `n, n+4, n+14` is divisible by 3 — giving
   `R_3^HL = 0` and empirically only the single triple containing `3`
   survives). PPP predicts `R_3^PPP ≈ 1.02 × 10^-3`, off by
   essentially infinity. **PPP cannot represent multi-point
   admissibility constraints because pair-level R_2 sees only 2-point
   admissibility.**

**F3 HOLDS:** 18 of 19 admissible all-even triples have
`|R_3^PPP - R_3^HL|/R_3^HL > 10%`; max gap = 79.16%.

---

## F4: Real signed K also fails

If we permit real K with arbitrary signs `s : t → {-1, +1}` and
magnitudes `|K(t)| = ρ √(S(0, t) - 1)` (fixed by R_2), the permanent
becomes

`R_3^PPP_signed = ρ^3 + ρ Σ K(t_i)^2 + 2 σ |K(t_1)||K(t_2)||K(dt)|`

where `σ = s(t_1) s(t_2) s(dt) ∈ {-1, +1}` for any sign assignment.

The cross-term ratio `σ_req` required to match HL is:

```
        triple    σ_req      |σ_req|     real-signed?    complex-feas?
     (0,2,6)    -0.5375        0.54         no               YES
     (0,2,8)    -0.5373        0.54         no               YES
     (0,2,12)   +0.3979        0.40         no               YES
     (0,4,6)    -0.5375        0.54         no               YES
     (0,4,10)   +0.4147        0.41         no               YES
     (0,4,16)   -0.5408        0.54         no               YES
     (0,6,8)    -0.5373        0.54         no               YES
     (0,6,12)   -0.0588        0.06         no               YES
     (0,6,16)   +0.4152        0.42         no               YES
     (0,6,18)   -0.0601        0.06         no               YES
     (0,8,14)   -0.0025        0.00         no               YES
     (0,8,18)   +0.4104        0.41         no               YES
     (0,10,16)  +0.4152        0.42         no               YES
     (0,10,22)  +0.5836        0.58         no               YES
     (0,12,18)  -0.0601        0.06         no               YES
     (0,14,20)  +0.7687        0.77         no               YES
     (0,16,22)  -0.2396        0.24         no               YES
     (0,18,24)  -0.0594        0.06         no               YES
     (0,20,26)  +0.5509        0.55         no               YES
```

`|σ_req| ∈ [0.002, 0.769]` at every triple — never near `±1` (closest
is 0.769 at `(0, 14, 20)`, still 23 % below `±1`).

**F4 HOLDS:** real-signed K fails at **every** admissible triple (0
of 19 satisfy `σ_req ∈ {±1}` within tolerance 0.05). The pair-level
magnitudes `|K(t)| = ρ √(S - 1)` are *too large* for any sign
assignment to reproduce HL.

---

## F5: Complex Hermitian K — phase consistency

Since `|σ_req| ≤ 1` at every triple, a *single* triple can be locally
fitted by a complex Hermitian K with magnitude fixed and a phase
chosen so that `cos(φ(t_1) + φ(dt) - φ(t_2)) = σ_req`. The question:
does a *globally consistent* `φ : t → [0, 2π)` exist solving all 19
equations simultaneously?

Numerical least-squares fit of `φ` over 13 distinct offsets
`{2, 4, 6, …, 26}` (gauge-fixing `φ(2) = 0`) using 200 random starts
of Levenberg-Marquardt + trust-region:

  **Best max residual `|cos(φ(t_1) + φ(dt) - φ(t_2)) - σ_req|`: 0.0746.**

This is small but well above sample noise (residuals scale linearly
with `cos(...)`-error, and the sigma_req values are computed from
exact HL singular series with truncation error <10⁻⁵). The 0.075
plateau represents a real obstruction.

**F5 HOLDS within tolerance 0.01 (= 1% match):** no complex Hermitian
phase assignment passes — best residual is 7.5%, off by ≥7×.

The worst residual is at `(0, 4, 10)` where `dt = 6`. The phase
constraint `φ(4) + φ(6) - φ(10)` is over-determined by triples like
`(0, 4, 10), (0, 6, 16), (0, 10, 16)` (the latter two also constrain
`φ(6)`), and the constraints are mutually inconsistent.

So even with the *most flexible* translation-invariant complex
Hermitian K of magnitude `|K(t)| = ρ √(S - 1)`, no phase choice
matches HL globally.

---

## Why D7 fails by structure: pair-factorisation vs prime-factorisation

The deeper structural fact behind F3-F5:

  - The HL singular series **factorises over primes** p:
       `S(0, t_1, t_2) = ∏_p α_p(0, t_1, t_2),`
    with `α_p` depending on the *number of distinct residues*
    `ν_p({0, t_1, t_2}) ∈ {1, 2, 3}`.

  - DPP / PPP 3-point correlations **factorise over pairs** `(t_i, t_j)`:
       `R_3^DPP/PPP = polynomial( K(t_1), K(t_2), K(t_2-t_1) )`,
    each pair K is a univariate function of the offset.

  - **Pairwise admissibility ≠ triple admissibility.** Example:
    `(0, 4, 14)` is admissible mod 3 in every pair but inadmissible
    as a triple (covers all of `Z/3Z`). PPP cannot detect this.

  - Even on admissible triples, the prime-factor cancellation pattern
    in `α_p` (especially at small primes where `ν_p` saturates to `p`
    or drops below) does not match the pair-product structure
    required by DPP/PPP.

Conclusion: **the integer prime sequence is NOT a translation-
invariant determinantal nor permanantal point process, real OR
complex, in any of its standard formal senses.** The HL prime-by-
prime factorisation is genuinely non-pairwise.

---

## Falsification outcomes

| ID | Hypothesis | Outcome at N = 10^7 |
|----|------------|---------------------|
| F1 | DPP pair-level: `K²_DPP < 0` at every admissible even t | **HOLDS** (15/15 violations) |
| F2 | PPP pair-level: `K²_PPP < 0` at every odd t > 1 | **HOLDS** (14/14 violations) |
| F3 | PPP 3-point: `\|gap\| > 10%` at some all-even triple | **HOLDS** (18/19 triples; max 79.16%) |
| F4 | Real-signed K: `σ_req ∉ {±1}` at every triple | **HOLDS** (19/19; max \|σ_req\| = 0.77) |
| F5 | Complex Hermitian K: no global φ matches | **HOLDS** (best residual 0.075 ≫ 0.01 tol) |

All 5 pre-stated falsifiers HOLD. **Primes are quantitatively NOT a
translation-invariant DPP, PPP, signed-real-K, or complex-Hermitian-K
point process** at the all-even sub-lattice 3-point level.

---

## Self-grade and edge proposal

**Self-grade: B-grade (structural negative result).** Per
`ATTACK_VECTORS.md §D7`:

> B-grade success: primes are NOT a DPP — quantitative breakdown of
> the determinantal identity at a specific (k, t) tuple, producing a
> new structural edge "primes are anti-determinantal beyond pair
> correlation."

This session delivered:

  - **Quantitative breakdown** at the 3-point level (max 79.16% gap).
  - **Multi-tuple breakdown**: 18 of 19 tested admissible triples.
  - **Stronger structural fact**: even the most flexible complex-K
    phase assignment fails globally.
  - **Mechanism identified**: HL prime-by-prime factorisation cannot
    be reduced to pair correlations because triple admissibility
    (`ν_p = p` saturation) is a multi-body fact.

**Adds candidate EDGE E2.16:**

> The integer prime sequence is NOT a translation-invariant
> determinantal nor permanantal point process. The 3-point Hardy-
> Littlewood singular series `S(0, t_1, t_2)` does NOT factorise as
> a permanent / determinant of pair-correlation factors `K(t_i, t_j)`,
> contradicting the structural prediction of any real DPP / PPP /
> signed-real-K / complex-Hermitian-K fit. Maximum
> `|R_3^PPP - R_3^HL| / R_3^HL = 79.16%` on the 3-AP triple (0, 6, 18)
> at N = 10^7. Distinct from prior pair-level closures (E2.13 Gowers,
> E2.14 Anderson, E2.15 algebraic immunity) which are 2-point
> statements; E2.16 is the first **3-point structural** statement
> ruling out a kernel-factorisation representation.

---

## Algorithmic content

None new for polylog π(x). PPP/DPP fit fails, so no kernel-based
polylog representation is available. Note also that even if a fit
existed, K(0) = ρ = 1/log N gives `π(N) ≈ N/log N` — the prime number
theorem trivially, no further information.

---

## Falsified things

- (Original §D7 goal) primes ARE a DPP with kernel admitting closed
  form: definitively falsified.
- (Stronger version) primes ARE a PPP: falsified at pair level (odd t).
- (Stronger version) primes ARE a real signed K-PPP on the all-even
  sub-lattice: falsified at 3-point level (no σ_req ∈ {±1}).
- (Strongest version) primes ARE a complex Hermitian-K PPP on the all-
  even sub-lattice: falsified at 3-point level (best residual 7.5 %).

## Things this session DID NOT prove

- Whether a non-translation-invariant K (depending on `(x, y)` not just
  `(x - y)`) could match. Likely, given enough freedom — but such K is
  not the "natural" object in DPP theory and would not give a polylog
  algorithmic content.
- Whether a Pfaffian point process (matrix-valued K, like 2-step
  determinantal) could match at higher k. Open but speculative.
- Whether the correct framework is `α-determinantal` (Vere-Jones 1997)
  with non-uniform `α(t)` ∈ ℝ — out of scope; unlikely to give polylog.

---

## Cross-domain technique status

**Update `CROSS_DOMAIN_TECHNIQUES.md §3`:** Determinantal Point
Processes — was PROPOSED (D7), now **USED (mode I)** with edge E2.16.
Permanantal counterpart added as a sub-row. Survey ref:
Hough-Krishnapur-Peres-Virag 2009 + Soshnikov 2000.

---

## Self-extension (per CLAUDE.md autonomy invariant)

Two follow-on challenges proposed for `NOVELTY_CHALLENGES.md`:

- **§D7.b — Pfaffian / matrix-valued K fit.** A Pfaffian point process
  has a 2x2 matrix kernel; its 2- and 3-point correlations admit
  signs and richer phase structure. Test whether the W-tricked
  `χ_{210, 1}` (where pair correlations approach 1) admits a Pfaffian
  representation that factorises through pair-level data. Predicted:
  same structural failure as D7 because admissibility is multi-body
  even after W-tricking.

- **§D7.c — α-determinantal fit with offset-varying α(t).** Allow
  `R_2(t) = ρ^2 - α(t) K(t)^2` with `α(t) ∈ ℝ`. Then α(t) ∈ {-1, +1}
  encodes parity. Could the 3-point identity match HL with a
  non-trivial α(t1, t2) function? Predicted: also fails (the HL
  factorisation over primes is the obstruction, not the sign).

Or as an A-grade attempt: §D7.d — derive a *non-pair-factorisable*
K-like object whose multi-point statistics reproduce HL. Goal: a new
arithmetic invariant of `χ_P` that factorises over primes (as HL does)
rather than over pairs (as DPP/PPP do). If successful, this would be
the first project object that genuinely captures the *prime-by-prime*
structure rather than reducing to pair correlations. Speculative; budget 2-3 sessions.

---

## Files

  - `primes_dpp_ppp_fit.py` — main script.
  - `primes_dpp_ppp_fit_results.md` — this file.
  - `main_run.json` — raw data (pair rows + triple rows + sigma test +
    LS fit residuals).

## Edges cited or composed

- **E1.10, E3.13, E7.1** (pseudorandomness battery; χ_P matches GUE
  pair correlation but NOT 3-point as we now show structurally).
- **E2.13** (Gowers U^k → HL singular series; E2.16 is the 3-point
  structural complement).
- **E2.14** (Anderson Lyapunov; E2.16 is the 3-point version of the
  same "structure = HL" statement).
- **E2.15** (algebraic immunity; E2.16 is the next-leg confirmation in
  a distinct mathematical category).

## Cross-domain references

- Hough-Krishnapur-Peres-Virag 2009 *Zeros of Gaussian Analytic
  Functions and Determinantal Point Processes*
  (AMS ULect 51) https://www.ams.org/books/ulect/051/
- Soshnikov 2000 "Determinantal random point fields"
  Russ. Math. Surv. 55, 923 https://arxiv.org/abs/math/0002099
- Vere-Jones 1997 (α-determinantal generalisation)
- Hardy-Littlewood 1923 (singular series for prime k-tuples)
- Wikipedia: Determinantal point process
  https://en.wikipedia.org/wiki/Determinantal_point_process

---

## Honest assessment

The session set out to attempt an A-grade outcome (find a kernel
`K(t)` that locally encodes prime correlations and gives a new
polylog representation). It failed cleanly: DPP and PPP both rule
out at pair level, and even with maximum flexibility (complex
Hermitian, magnitudes from R_2) the 3-point identity has a 7.5%
residual obstruction.

The B-grade outcome was anticipated by the §D7 attack vector ("primes
are NOT a DPP"); we exceeded it by ALSO ruling out PPP, real-signed,
and complex-Hermitian variants, and by quantifying the structural
mechanism (prime-factorisation vs pair-factorisation, plus
multi-body admissibility).

This is structurally consistent with the project's prior closures
(χ_P deviation from random reduces to HL, no further structure) but
is the FIRST 3-point statement of that fact, complementing the
2-point closures at E2.13 / E2.14 / E2.15. The cross-domain
technique (DPP) is genuinely fresh.

Honest grade: **B**. No A-grade content. No inflated claims.

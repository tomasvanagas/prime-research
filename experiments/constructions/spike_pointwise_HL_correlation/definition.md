# T_Q autocorrelation as truncated Hardy–Littlewood singular series

## Object

For positive integer N, integer cutoff Q ≥ 1, and integer shift h ≥ 0,
define the **shift-h autocorrelation of the pointwise spike approximator**

```
    R_h(Q, N)  :=  (1/(N - h)) · Σ_{n=1}^{N - h}  T_Q(n) · T_Q(n + h),
```

where `T_Q(n) = (π(N)/N) · Σ_{q sqf ≤ Q} (μ(q)/φ(q)) · c_q(n)` is the
S191 / C9 pointwise approximator (`experiments/constructions/
ramanujan_spike_pointwise/`). The **connected** autocorrelation
subtracts the squared mean (the q = 1 background mode):

```
    R_h^conn(Q, N)  :=  R_h(Q, N)  −  (π(N) / N)^2.
```

## Identity used (proven from definitions, verified empirically)

By the Ramanujan-Fourier orthogonality

```
    (1/N) Σ_{n=1}^{N} c_{q1}(n) c_{q2}(n + h)
        =  c_q(h)             if q1 = q2 = q,
        =  0                   if (q1, q2) coprime, both ≥ 2.
```

(the second relation holds when q1 q2 | N; finite-N correction below).
Plugging the diagonal into the Ramanujan-Fourier expansion of T_Q,

```
    T_Q(n)  =  (π(N)/N) · Σ_{q sqf ≤ Q}  (μ(q)/φ(q)) · c_q(n),
```

the diagonal q1 = q2 contribution to the average is

```
    R_h^diag(Q, N)
        =  (π(N)/N)^2  ·  Σ_{q sqf ≤ Q}  (μ(q)/φ(q))^2 · c_q(h)
        =  (π(N)/N)^2  ·  S_Q(h).
```

The function

```
    S_Q(h)  :=  Σ_{q sqf ≤ Q}  (μ²(q) / φ²(q)) · c_q(h)
```

is the **conductor-Q truncation of the Hardy–Littlewood twin-prime
singular series**, which has the classical Ramanujan-Fourier expansion

```
    S(h)  =  Σ_{q ≥ 1}  (μ²(q) / φ²(q)) · c_q(h)
```

(see Hardy–Littlewood 1923; Wintner 1944). The infinite-Q limit is

```
    S(h) =
        2 C_2 · ∏_{p | h, p ≥ 3} (p − 1) / (p − 2)     for even h ≠ 0,
        0                                              for odd h,
        +∞                                              for h = 0,
```

where `C_2 = ∏_{p ≥ 3} (1 − 1/(p − 1)²) ≈ 0.660 161 815 8…`.

## Why this object exists (project-internal derivation)

S191 / C9 built `T_Q(n)` as a closed-form pointwise spike approximator
of `chi_P` matching the S168 / S169 L²-energy spike content. The
**single-point** L² content was established. The natural next moment is
the **two-point** content:

```
    < chi_P(n) · chi_P(n + h) >_n  ~  S(h) · (π(N) / N)²
```

(Hardy–Littlewood twin-prime conjecture) for even h ≠ 0. Since T_Q is a
truncation of `chi_P`'s Ramanujan-Fourier projection onto `⊕_{q sqf ≤
Q} V_q^prim`, its two-point function should reproduce the
**Q-truncated** singular series, identically, up to the controlled
cross-conductor leakage and the Gallagher S168 remainder.

This composes:

* **E2.1** — MPS bond-dim identity (S168 spike subspace lives at
  squarefree conductors);
* **E2.13** — Gowers `U^k` of chi_P matches HL singular series (k = 2
  cube version, which is the four-point cube; this construction is the
  **2-point** companion);
* **E1.6** — parity bisection (q = 2 term `μ(2)/φ(2) c_2(h) = c_2(h) =
  +1` for even h, `−1` for odd h: the parity pattern is recovered as
  the leading non-constant term of `S_Q(h)`);
* **E2.2** — Liouville / parity identity (the q = 2 contribution is
  exactly the Liouville parity sign at shift h);
* **C9 / S191** — pointwise spike approximator (the object being squared).

## Predicted moments (closed-form)

For h ≥ 1:

```
    R_h^conn(Q, N)  ≈  (π(N)/N)² · S_Q(h)             [main prediction]
```

with `S_Q(h)` the truncated singular series above. Cross-conductor
leakage (from squarefree q1 ≠ q2 with shared prime factor) and the
finite-N remainder contribute a `O((π(N)/N)² · Q · log Q / N)` error.

For h = 0:

```
    R_0(Q, N)  =  ⟨T_Q²⟩  =  (π(N)/N)² + (π(N)²/N) · Σ_{q sqf ≤ Q, q ≥ 2} 1/φ(q),
```

which is the energy formula recovered from S190 / S191. **This
provides a built-in self-consistency check**: the C9 L²-content
prediction must match the h = 0 column of this construction.

For h = 2 and Q → ∞:

```
    S(2) = 2 C_2 ≈ 1.32032 36313…
```

so `R_2^conn(Q → √N, large N) ≈ (π(N)/N)² · 1.3203…`.

For odd h ≥ 1: `S(h) = 0`. The truncated `S_Q(h)` may oscillate around
zero at small Q (q = 1 contributes +1, q = 2 contributes −1, finite-Q
truncation determines the sign), but should approach 0 as Q grows. The
prime correlation `π_h := ⟨χ_P(n) χ_P(n+h)⟩` is also 0 in the limit,
so this is a non-trivial agreement at h odd.

## Falsification criteria (pre-stated, to be checked AFTER running)

1. **(F1) Identity at moderate Q.** For h ∈ {2, 4, 6, 8, 10, 12, 30,
   210} and Q ∈ {30, primorial(4) = 210, primorial(5) = 2310,
   round(N^{0.185})}, at d ∈ {16, 18, 20}, the ratio

   ```
       R_h^conn(Q, N) / [(π(N)/N)² · S_Q(h)]
   ```

   must lie in `[0.85, 1.15]` (10–15% relative band, allowing for
   cross-conductor leakage and finite-N).

2. **(F2) Convergence to full HL.** As Q → √N, the empirical
   `R_h^conn / (π(N)/N)²` must converge to the textbook full singular
   series `S(h) = 2 C_2 ∏_{p|h, p≥3} (p−1)/(p−2)` for even h ∈ {2, 4,
   6, 8, 10, 12, 30}. Relative agreement < 5 % at d = 18, 20, Q = √N.

3. **(F3) Odd-h vanishing.** For h ∈ {1, 3, 5, 7, 9}, both
   `R_h^conn / (π(N)/N)²` and the textbook `S(h) = 0` must agree to
   `|·| < 0.1` (an order of magnitude below the typical even-h
   `S(h) ~ 1.3 - 2.5`).

4. **(F4) h = 0 self-consistency.** `R_0(Q, N) − (π(N)/N)²` must equal
   `(π(N)²/N) · Σ_{q sqf ≤ Q, q ≥ 2} 1/φ(q)` to within 5 % at d = 18,
   20.

5. **(F5) Prime baseline.** The empirical prime-correlation
   `π_h(N) := (1/(N − h)) #{n ≤ N − h : n, n + h prime}` must match the
   HL prediction `(π(N)/N)² · S(h)` to within 5–10 % at d = 20, h ∈
   {2, 4, 6, 8, 10, 12, 30, 210} (HL is known to hold empirically at
   this scale; this is a sanity check, not a novel claim).

A construction that fails F1 but passes F2 is **partial** (the
identity holds asymptotically but cross-conductor leakage at finite Q
breaks the strict equality). A construction that fails both is
duplicate-PLUS of E2.13.

## Edges composed (cite by ID)

* **E2.1** (MPS bond-dim identity) — squarefree-conductor spike
  subspace.
* **E2.2** (Liouville / parity identity) — q = 2 term recovers parity
  sign.
* **E2.13** (Gowers `U^k` of chi_P matches HL singular series) — the
  family this construction extends to the **2-point** case.
* **E1.6** (parity bisection) — `c_2(h)` distinguishes even / odd h.
* **C9 / S191** (pointwise spike approximator) — the squared object.

## What is novel relative to E2.13 + S191

* E2.13 is the **4-point cube** (Gowers `U^2`) singular-series match.
  This construction is the **2-point shift** singular-series match. The
  algebraic structure (Ramanujan-Fourier, `μ²/φ²`) is identical, but
  the 2-point version is `O(N)` to evaluate (vs `Θ(N² log N)` for `U^2`).
* S191 is a **single-point** L² statement. This construction is the
  **two-point** generalisation. The two-point identity carries strictly
  more information (the diagonal recovers S191's L² content as the
  h = 0 case).
* The closed-form pointwise function `T_Q(n)` ⇒ truncated HL singular
  series via simple shift-and-multiply averaging is, to my knowledge,
  not previously catalogued in the project.

## Algorithmic content

* **Cost** for the autocorrelation curve `{R_h(Q, N) : h ∈ H}`:
  `O(N · |H| + Q · ω · N)` (the latter for the T_Q array).
* This is the second-moment fingerprint of T_Q. The full Hardy–
  Littlewood signature is recovered with `Q ≪ √N` (sub-`√N` cost).
* **No polylog opening expected.** Same as S191: the singular-series
  recovery is structural, not algorithmic.

## Reproducing

```
cd experiments/constructions/spike_pointwise_HL_correlation
python3 tq_correlation.py            # main sweep, ~5 min at d=20
python3 sanity_singular_series.py    # standalone S(h) numerical check
```

# Critique — Session 63 (2026-04-26)

Critique of the four S63-fresh proposals (`archive/ephemeral/proposals_session.md`)
and their experimental verdicts (`experiments/proposals/session63fresh_*_results.md`).

The S63 session, instructed not to consult CLOSED_PATHS.md while drafting, ran each
of its four proposals and self-reported all four as CLOSED. My job here is to verify
those verdicts against the 537+-entry closed-paths catalogue and decide whether any
genuinely novel component survives.

**Bottom line:** all four self-verdicts hold up. P1 is a **legitimate new closure**
of a sub-question not previously tested at this granularity (D-finiteness specifically
on δ(n), not on π(n) or 1_P). P2 and P4 are **duplicate-plus** of prior closures
(line 693 and line 685 respectively). P3 is **CONFIRMED CLOSED (circular)** with
one empirical refinement worth recording. Nothing escalates to `novel/` or
`OPEN_PROBLEMS.md`.

---

## Verdict summary

| # | Proposal | Self-verdict | Critic verdict | Mode | Closest prior |
|---|----------|--------------|----------------|------|----------------|
| P1 | D-finite recurrence on δ(n) | CLOSED (I) | **CONFIRMED-NOVEL TEST, CLOSED** | I | extends lines 576, 577, 680 (D-finite on π and π/n GF) |
| P2 | Mollifier-corrected explicit formula | CLOSED (E) | **DUPLICATE-PLUS, CLOSED** | E | line 693 (Hermite/Gaussian/Riesz mollification) |
| P3 | RMT local-moment predictor for Δ(x) | CLOSED-circular (C) | **CONFIRMED, CLOSED** with side-finding worth keeping | C | lines 353, 656 |
| P4 | Newton with progressive 2^k zero-budget | CLOSED (I) | **DUPLICATE-PLUS, CLOSED** | I | line 685 (R⁻¹ fixed-point with zero correction) |

---

## P1 — D-finite (Apéry-style) recurrence hunt for δ(n)

### Has the exact approach been tried?

**Adjacent but distinct.** Three near-misses in CLOSED_PATHS.md:

- **Line 576** (S23): "Holonomic (D-finite) recurrence for π(n) — FAIL (I), order ≤ 20,
  polynomial degree ≤ 8". Target was π(n).
- **Line 577** (S23): "Prime indicator holonomic — FAIL (I), same ratio to random as
  π(n)." Target was 1_P(n).
- **Line 680** (S43): "J-fraction of π(n)/n generating function (D-finiteness test) —
  FAIL (E)". Hankel-determinant route on π(n) generating function.
- **Line 703** (S49): PSLQ + zeta-zero/log basis on δ(n). Algebraic, not D-finite.

So the **specific test** "is δ(n) = p(n) − R⁻¹(n) D-finite at low (L, d)?" has not
appeared before. P1 is a **legitimate new closure**, not a duplicate, even though
the broader holonomic-on-prime-counting question is closed.

### Failure-mode analysis

Mode **I**. Held-out skill 0.5–1.3 across (L,d) ∈ [1..4]×[1..4] is at the noise floor.
The train-side rank ratio shrinking with d is **purely a column-conditioning artifact**
(rows in n^k δ(n+j) span 12 orders of magnitude). The proposer caught this and switched
to held-out prediction, which is the right test. Methodology is sound.

### Complexity claim

Correct: if a recurrence of order L and polynomial-degree d existed, Bostan–Salvy–Schost
binary-splitting shift would evaluate at index n in O((L+d)·M(log n)). The complexity
analysis is fine; the issue is empirical absence of such a recurrence.

### Specific obstacle

δ(n) inherits GUE-random oscillations from ζ-zeros (`novel/pseudorandomness_of_pi.md`,
21+ measures). Generic prediction: δ-style sequences with GUE-distributed Fourier phases
are not P-recursive at any low order, consistent with line 576 closing at d ≤ 8 for π
and line 703 closing PSLQ on δ.

### Verdict — **CONFIRMED CLOSED, mode I.** Add to CLOSED_PATHS as a new entry.

### What's worth keeping

The proposer's **methodology lesson** is genuinely useful: when running null-space
searches on data with multiplicative column-scale spread, validate via held-out
**prediction**, not via singular-value ratios on the training matrix — column conditioning
can fake rank deficiency at 1e-13. Worth a one-line record for any future PSLQ /
D-finite / regression search in this project.

---

## P2 — Mollifier-corrected explicit formula

### Has the exact approach been tried?

**Yes — covered structurally by line 693 (S48):**

> "Hermite/Gaussian/Riesz mollification of truncated explicit formula for ψ(x) — FAIL (E).
> Tested Gaussian kernel exp(-(γ-log x)²/2σ²) for σ ∈ {0.5, 1, 2} and Riesz mean
> (1-|γ|/T)² for T ∈ {50, 200, 800}, K up to 800 zeros. Linear-functional argument:
> any kernel w(γ) with sum|w| ≥ cN(T) inherits Ω(√x · √(log T / log x)) tail bound;
> σ → ∞ recovers unmollified. To beat √x, need NONLINEAR operation on zeros, not
> re-weighting."

P2's Selberg-style **Dirichlet-polynomial** mollifier M(s) = Σ_{n≤Y} a_n n^{-s} is a
*different shape* of kernel (multiplicative-arithmetic structure) but is captured by the
same general-kernel tail-bound argument. The proposer's trade-off observation
("Y growth pays back what K kills") is the manifestation of that tail bound on the
specific kernel shape.

Adjacent prior closures: lines 31, 32 (Beurling–Selberg uncertainty), 36 (sieve-style
mollification), 519, 685.

### Methodology issue (worth flagging)

**P2's experiment does NOT actually run the explicit formula for ζ·M.** The proposer
explicitly notes this in the script docstring (lines 37–49 of
`session63fresh_mollifier_pi.py`): the test re-weighs each ζ-zero contribution by
M(ρ_j)/M(1) without including the contribution from M's own zeros/poles or the cross
terms in the von-Mangoldt-type coefficients of ζ'/ζ + M'/M. So the empirical
"mollifier worse than sharp" result is **not** decisive evidence about the actual
ζ·M explicit formula — only that the naive heuristic weighting doesn't help.

The **decisive closure** is line 693's theoretical argument (kernel-w bound), which
applies regardless of which kernel shape one picks. The empirical run is corroborating
evidence at best.

### Verdict — **DUPLICATE-PLUS (refines line 693), CLOSED, mode E.**

Worth a CLOSED_PATHS entry because the *Dirichlet-polynomial* mollifier shape is not
listed by name, but with explicit citation of line 693 as the parent closure.

---

## P3 — RMT local-moment predictor for Δ(x)

### Has the exact approach been tried?

**Adjacent but not identical.** Two prior entries cover the GUE/RMT route:

- **Line 353** (S8/S10): "RMT/GUE approximation — FAIL (I). Correct statistics, wrong values."
- **Line 656** (S32): "Hilbert-Pólya trace / GUE model / functional equation — FAIL (E,I).
  GUE: correct statistics, wrong values (std ≈ 20-60). δ(x) incompressible."

P3's **specific** twist — using a *local empirical window* of Δ(x±k) for k=1..H and a
weighted-mean estimator to predict Δ(X) — has not been tried this way. The empirical
result (RMSE 0.33 with H=100) is stronger than line 353 alone would predict because it
exploits **local smoothness of Δ on small scales**, which is a real phenomenon.

### The circularity argument

Verified. Δ(x) = π(x) − li(x), so building the window samples requires π(x) at every
window point. Two strategies:

1. **Global sieve up to X:** O(X log log X) — strictly worse than Meissel–Lehmer's
   O(X^{2/3}).
2. **Anchor + Miller–Rabin sweep:** anchor π(x_0) costs O(x_0^{2/3}) once, then
   incremental sweep through (x_0, X+H] is O(H · log² X). For a *single* π(X) query
   this is dominated by the anchor.

The proposer's analysis is correct. **The only path that would convert this into a
polylog algorithm is a polylog anchor — which is the original problem.** So C is the
right mode.

### Side-finding worth recording

The empirical observation — **std(Δ) over a 200-window is < 2 even at X = 5×10^4** —
is a quantitative refinement of the GUE prediction "Δ varies slowly on scales below √x".
Specifically, weighted-mean predictor RMSE 0.33 over windows of size 200 says:

> Δ(X) is determined to within ±0.5 by Δ(x) on a window of 200 points around X,
> for X up to 5×10^4.

This supports (does not prove) the conjecture that **the entropy of Δ at scale X is
concentrated at frequencies ≥ 1/√X**, not at unit frequency. It's not strong enough for
a `novel/` entry on its own (the underlying "GUE → smooth on small scales" claim is in
the literature), but worth a CLOSED_PATHS line. **Caveat:** would need a Cramér-random
control to confirm this is a *prime-specific* phenomenon and not a generic property of
any zero-mean drift process at scale √x. Flagged for follow-up next session.

### Verdict — **CONFIRMED CLOSED, mode C.** Useful side-finding.

---

## P4 — Newton with progressive 2^k zero-budget

### Has the exact approach been tried?

**Yes, structurally — line 685 (S44):**

> "Fixed-point iteration p_{k+1} = R⁻¹(n + sum_ρ R(p_k^ρ)) — FAIL (I+E). Iteration
> is non-contracting and non-monotone. n=10: k=2 hits err=-0.14 then k=3 drifts back
> to -1.89. n=1000: k=1 hits err=-0.16, k=2/k=3 drift to -0.6. ... Map f(x) =
> R⁻¹(n + S(x)) has |f'| oscillating ~ |S'(x)|, GUE-random; no contraction."

P4's variant uses **Newton-on-π** rather than R⁻¹-fixed-point, and adds **progressive
zero budget** K_k = 2^(k+1). Mathematically these are different maps but the **failure
mechanism is identical**: π_K(x) carries O(1) noise from the GUE-random tail, Newton
amplifies it by 1/π'(x) ≈ log(x), and for n ≥ 1000 this noise (~log x ≈ 9) is at the
scale of prime gaps. So Newton oscillates instead of converging.

Adjacent prior closures: line 36 (heuristic candidate gen + AKS), line 337
(self-correcting R⁻¹), line 657 (convergence acceleration of zero sum).

### The "geometric K_k" claim

Even granting K_total = O(K_final), the final K must satisfy
|π_K(p(n)) − n| < 0.5 / log(p(n)) ≈ 0.05 (so that Newton's per-iter step is < 0.5).
The empirical zero-tail behaves as |π_K(x) − π(x)| ~ √x / (√K · log x · √|γ_K|).
Plugging: x = 10^4, K = 256 ⇒ tail ≈ 100 / (16 · 9 · 6) ≈ 0.12 — at the threshold.
x = 10^6, K = 256 ⇒ tail ≈ 1000 / (16 · 13 · 6) ≈ 0.8 — exceeds threshold.

So the claimed budget K_total = O(log² n) **does not** suffice asymptotically; one
needs K = Ω(x^{1/2-ε} / polylog), and Newton-on-π_K cannot localize p(n) to ±0.5
within polylog total work. This is the same √x barrier as direct explicit-formula
evaluation. The proposer's empirical observation (oscillates at n ≥ 1000 in a band of
size 5–25) is consistent with this analysis.

### Verdict — **DUPLICATE-PLUS (refines line 685 via different fixed-point map), CLOSED, mode I.**

### Quantitative obstruction worth keeping

The proposer's identification of the failure mode — π_K-noise × Newton's 1/log x
amplification ≈ prime gap scale at x ≥ 10^3 — is a clean quantitative statement that
complements line 685's "|f'| oscillates GUE-randomly" framing.

---

## Cross-cutting observations

1. **All four proposals self-closed correctly.** The proposer's experimental discipline
   was high: each test was hardened against a specific failure mode (P1: column
   normalization + held-out prediction; P3: control by varying H and X; P4: progressive
   K with explicit residual tracking). No reruns needed.

2. **None of the four warrants `novel/`.** Each closure either duplicates or refines an
   existing closure.

3. **None of the four reopens an OPEN_PROBLEMS.md direction.** OPEN_PROBLEMS still has
   only Circuit Complexity of π(x) as the sole genuinely viable direction (with FOCUS-1
   sub-attacks 1 and 3 still un-built per line 714).

4. **The strongest by-product is P3's empirical local smoothness of Δ.** This *could*
   eventually belong in `novel/pseudorandomness_of_pi.md` as a **complementary** measure
   (local smoothness below √x scale is the *flip side* of pseudorandomness at scale √x).
   Recommend adding one paragraph there, but only **after** running a Cramér-random
   baseline. Skipping for this critique session because that requires new code; flagging
   for next session.

5. **Methodology lesson worth pinning:** P1's column-normalization vs held-out validation
   distinction is general — applies to any future PSLQ / null-space / regression search.
   Adding a note to the project conventions (or a one-liner in CLAUDE.md) might save
   future sessions from the same trap.

---

## Action items

The proposer already wrote four `_results.md` files. The remaining work for this
critique:

1. **Add four entries to `status/CLOSED_PATHS.md`** (S63 batch). Drafted below.
2. **Optionally add P3's local-smoothness measure** to `novel/pseudorandomness_of_pi.md`
   *after* running a Cramér baseline. Deferred to next session.
3. **Update `status/SESSION_INSIGHTS.md`** with S63 entry.
4. **Append session synthesis** at `archive/sessions/session63_proposals_critique.md`.

### Drafted CLOSED_PATHS.md entries (insert after line 714)

```
| D-finite (Apéry-style) recurrence for δ(n) = p(n) − R⁻¹(n) at (L,d) ≤ (4,4) (S63 fresh) | FAIL | I | Held-out prediction skill 0.5–1.3 (≈ noise) for all (L,d) tested at N=400, 200-bit precision. Train-side singular-value drop with d is purely column-conditioning artefact (rows span 12 OOM in n^k); held-out validation is the correct test, training rank-ratio is artifact. Distinct from line 576 (D-finite on π itself) and line 680 (J-fraction on π/n GF) — first explicit test on the analytic-detrending residual δ. δ inherits GUE-random oscillation; not P-recursive at low order. Methodology lesson: validate null-space relations via prediction, not training SVD ratios. See experiments/proposals/session63fresh_dfinite_delta.py | 63 |
| Selberg-style Dirichlet-polynomial mollifier weight on truncated explicit formula (S63 fresh) | FAIL | E | DUPLICATE-PLUS of line 693 (Hermite/Gaussian/Riesz mollification, S48). Tested length-Y mollifier M(s)=Σ_{n≤Y} a_n n^{-s} with a_n LSQ-tuned to zero |M(ρ_j)| for j ≤ K, sweep (K,Y) ∈ {(5,30),(10,50),(20,100),(40,200)}, T ≤ 1000. Mollified estimator strictly WORSE than sharp at every T (sharp at x=10000, T=50: err 0.44; best mollified: 1.72). CAVEAT: experiment uses naive M(ρ)/M(1) per-zero weighting without ζ·M's M-zero/pole contributions; decisive closure is line 693's general-kernel argument (any kernel w with Σ|w| ≥ cN(T) inherits Ω(√x √(log T/log x)) tail bound). Y-vs-K duality: any K zeros killed cost length-Y compensation, no net gain. To beat √x need NONLINEAR ops on zeros. See experiments/proposals/session63fresh_mollifier_pi.py | 63 |
| RMT local-moment Bayesian predictor for Δ(x) = π(x) − li(x) over H=100 window (S63 fresh) | FAIL | C | Weighted-mean predictor (1/|x−X|) on Δ over [X−H, X+H] window predicts Δ(X) at RMSE 0.33 < 0.5 across X ∈ {1000, 5000, 10000, 50000} — mathematically passes ±0.5 criterion. Algorithm is CIRCULAR: building Δ(x) for x in window requires π(x) at every window point. Two strategies both fail: global sieve to X is O(X log log X) (worse than Meissel-Lehmer); anchor + Miller-Rabin sweep needs anchor π(x_0) at O(x_0^{2/3}) which IS the original problem. Side-finding (worth keeping for novel/pseudorandomness_of_pi.md after Cramér baseline): std(Δ) over 200-window < 2 even at X=5e4, supports "Δ entropy at scale X concentrated at freq ≥ 1/√X". Adjacent: lines 353 (RMT statistics-only), 656 (HP/GUE), 36 (anchor + AKS sieve). See experiments/proposals/session63fresh_rmt_moments.py | 63 |
| Newton-on-π with progressive 2^k zero-budget for inverting π(x)=n (S63 fresh) | FAIL | I | DUPLICATE-PLUS of line 685 (R⁻¹ fixed-point with zero correction, S44). Different map (Newton on π vs fixed-point on R⁻¹) but identical failure mechanism: π_K(x) carries O(1) GUE-random tail noise, Newton amplifies by 1/π'(x) ≈ log(x), so iterate band-oscillates at scale ~log x for n ≥ 1000. Empirical: n=100 converges in 2 steps with K=6 (err 0.016); n=1000/5000/10000 oscillate in band 5–25, n=5000 even drifts away from initial guess. Asymptotic obstruction: required K satisfies √x/(√K·log x·√γ_K) < 0.05/log x ⇒ K ≥ x^{1/2-ε}/polylog, identical to direct explicit-formula evaluation. Geometric K_k=2^k summing to K_total=O(K_final) does not break the √x barrier. See experiments/proposals/session63fresh_newton_zerobudget.py | 63 |
```

---

## Critic's bottom line

All four proposer self-verdicts hold up. Two are duplicates-plus of existing closures
(P2 → line 693; P4 → line 685), one is a new specific test that confirms a broader
pattern (P1 — first explicit D-finite test on δ rather than π), one is genuinely
circular but yields a small empirical refinement to the local-smoothness picture (P3).

The session adds four well-documented closures, one methodology lesson on
column-conditioning vs held-out validation, and a small empirical handle on Δ-smoothness
at scales below √x. **Nothing reopens any closed direction; nothing reaches the
`novel/` bar.** OPEN_PROBLEMS.md remains as it was: only circuit complexity of π(x)
is genuinely open.
